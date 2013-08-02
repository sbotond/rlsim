/*
* Copyright (C) 2013 EMBL - European Bioinformatics Institute
*
* This program is free software: you can redistribute it
* and/or modify it under the terms of the GNU General
* Public License as published by the Free Software
* Foundation, either version 3 of the License, or (at your
* option) any later version.
*
* This program is distributed in the hope that it will be
* useful, but WITHOUT ANY WARRANTY; without even the
* implied warranty of MERCHANTABILITY or FITNESS FOR A
* PARTICULAR PURPOSE. See the GNU General Public License
* for more details.
*
* Neither the institution name nor the name rlsim
* can be used to endorse or promote products derived from
* this software without prior written permission. For
* written permission, please contact <sbotond@ebi.ac.uk>.

* Products derived from this software may not be called
* rlsim nor may rlsim appear in their
* names without prior written permission of the developers.
* You should have received a copy of the GNU General Public
* License along with this program. If not, see
* <http://www.gnu.org/licenses/>.
 */

package main

import (
	"os"
	"runtime"
	"runtime/pprof"
	"time"
)

var L *Log
var Rg Rander

func main() {

	// Intialize logging on stderr:
	L = NewLog(1)

	// Parse command line arguments:
	args := new(CmdArgs)
	args.Parse()

	// Re-initialize logging:
	var logLevel int
	if args.Verbose {
		logLevel = 1
	}
	L = NewLog(logLevel)

	// Start up CPU profiling:
	if args.ProfFile != "" {
		f, err := os.Create(args.ProfFile)
		if err != nil {
			L.Fatalf("Could not create file \"%s\" for profiling output: %s", args.ProfFile, err.Error())
		}
		pprof.StartCPUProfile(f)
		defer f.Close()
		defer pprof.StopCPUProfile()
	}

	// Set the number of processors to use:
	runtime.GOMAXPROCS(int(args.MaxProcs))
	L.PrintfV("Starting up rlsim using maximum %d cores.", args.MaxProcs)

	// Initialize global random number generator:
	var seed int64
	if args.InitSeed != 0 {
		seed = args.InitSeed
	} else {
		time := time.Now().UTC()
		seed = time.Unix()
	}
	L.PrintfV("Initial random seed: %d", seed)
	Rg = NewRandGen(seed)

	// Report about raw parameter file:
	if args.RawParamsFile != "" {
		L.PrintfV("Using raw parameter file: %s", args.RawParamsFile)
	}

	// Initialize reporter:
	var report Reporter
	report = NewReport(args.ReportFile)

	// Initialize input:
	var input Inputer
	input = NewInput(args.InputFiles)

	// Initialize target:
	var target Targeter
	if args.RawLenProbs != nil {
		target = NewRawTarget(args.ReqFrags, args.RawLenProbs, Rg)
	} else {
		target = NewTarget(args.ReqFrags, args.TargetMix, Rg)
	}

	// Report target lengths:
	target.ReportTargetLengths(report)

	// Initialize fragmentor:
	var fragmentor Fragmentor
	fragmentor = NewFragmentor(args.FragMethod, args.PrimingTemp, args.KmerLength, target, args.FragLossProb)

	// Initialize fragment statistics:
	var stats FragStater
	stats = NewFragStats()

	// Initialize thermocycler:
	var cycler Thermocycler
	cycler = NewTechne(args.NrCycles, args.FixedEff, args.GcEffParam, args.RawGcEffs, args.MinRawGcEff, args.LenEffParam, target)
	// Report efficiency functions:
	cycler.ReportEffFunctions(target, report)

	// Initialize sampler
	var sampler Sampler
	sampler = NewLenSampler(args.StrandBias)

	//Initialize pool:
	var pool Pooler
	pool = NewPool(args.GobDir)

	// Deal with the PCR seed:
	var pcrRand Rander
	if args.PcrSeed != 0 {
		pcrRand = NewRandGen(args.PcrSeed)
		L.PrintfV("PCR random seed: %d", args.PcrSeed)
	} else {
		pcrRand = Rg
	}
	// We still need the old global generator wehen simulating fragmentation.

	// Initialize Transcripts
	pool.InitTranscripts(input, target, fragmentor, cycler, stats, args.GCFreq, args.PolyAParam, args.ExprMul, Rg, pcrRand)

	// Deal with sampling seed:
	if args.SamplingSeed != 0 {
		Rg = NewRandGen(args.SamplingSeed)
		L.PrintfV("Sampling random seed: %d", args.SamplingSeed)
	} else if args.PcrSeed != 0 {
		// If we had a PCR seed than replace the global random number generator:
		Rg = pcrRand
	}

	// Sample fragments:
	sampler.SampleFragments(pool, target, stats, Rg)

	// Report fragment statistics:
	stats.ReportFragStats(report)

	// Cleanup pool and transcripts:
	pool.Cleanup()

	// Write out report in JSON format:
	report.WriteJSON()
}
