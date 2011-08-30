/*
 *  Copyright (C) 2011 by Botond Sipos, European Bioinformatics Institute
 *  sbotond@ebi.ac.uk
 *
 *  This file is part of the rlsim software for simulating RNA-seq
 *  library preparation with PCR biases and size selection.
 *
 *  rlsim is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  rlsim is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with rlsim.  If not, see <http://www.gnu.org/licenses/>.
 */

package main

import (
	"runtime"
	"runtime/pprof"
	"time"
	"os"
)

var L *Log
var Rg Rander

func main() {

	// Intialize logging on stdout:
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
			L.Fatalf("Could not create file \"%s\" for profiling output: %s", args.ProfFile, err.String())
		}
		pprof.StartCPUProfile(f)
		defer f.Close()
		defer pprof.StopCPUProfile()
	}

	// Initialize global random number generator:
	time := time.UTC()
	Rg = NewRandGen(time.Seconds())

	// Set the number of processors to use:
	runtime.GOMAXPROCS(int(args.MaxProcs))
	L.PrintfV("Starting up rlsim is using %d cores.", args.MaxProcs)

	// Initialize reporter:
	var report Reporter
	report = NewReport(args.ReportFile)

	// Initialize input:
	var input Inputer
	input = NewInput(args.InputFiles)

	// Initialize target:
	var target Targeter
	target = NewTarget(args.ReqFrags, args.TargetMix, Rg)

	// Report target lengths:
	target.ReportTargetLengths(report)

	// Initialize fragmentor:
	var fragmentor Fragmentor
	fragmentor = NewFragmentor(args.FragMethod, args.PrimingBias, args.KmerLength, target)

	// Initialize fragment statistics:
	var stats FragStater
	stats = NewFragStats()

	// Initialize thermocycler:
	var cycler Thermocycler
	cycler = NewTechne(args.NrCycles, args.FixedEff, args.GcEffA, args.GcEffB, args.LenEffPar)
	// Report efficiency functions:
	cycler.ReportEffFunctions(target, report)

	// Initialize sampler
	var sampler Sampler
	sampler = NewLenSampler(args.StrandBias)

	//Initialize pool:
	var pool Pooler
	pool = NewPool(args.GobDir)

	// Initialize Transcripts
	pool.InitTranscripts(input, target, fragmentor, cycler, stats, args.GCFreq, args.PolyAMean, args.PolyAMax, Rg)

	// Sample fragments:
	sampler.SampleFragments(pool, target, stats, Rg)

	// Report fragment statistics:
	stats.ReportFragStats(report)

	// Cleanup pool and transcripts:
	pool.Cleanup()
}
