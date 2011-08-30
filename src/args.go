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
	"flag"
	"fmt"
	"strings"
	"os"
)

// Mixture component:
type MixComp struct {
	Mean uint64
	Sd   uint64
	Low  uint64
	High uint64
}

func (comp *MixComp) String() (s string) {
	s = fmt.Sprintf("(%d, %d, %d, %d)", comp.Mean, comp.Sd, comp.Low, comp.High)
	return
}

type TargetMix struct {
	Components map[*MixComp]float64
}

func (mix *TargetMix) String() (s string) {
	s += "["
	for k, v := range mix.Components {
		s += " " + fmt.Sprintf("%.2f", v) + ":" + k.String() + " "
	}
	s += " ]"
	return
}

type FragMethod struct {
	Name  string
	Param int
}

func (fm *FragMethod) String() (s string) {
	s = fmt.Sprintf("[ name: %s, param: %d ]", fm.Name, fm.Param)
	return
}

// Struct to hold command line arguments:
type CmdArgs struct {
	ReqFrags    int64
	NrCycles    int64
	StrandBias  float64
	PrimingBias float64
	FixedEff    float64
	GcEffA      float64
	GcEffB      float64
	LenEffPar   float64
	ReportFile  string
	Verbose     bool
	TargetMix   *TargetMix
	FragMethod  *FragMethod
	InputFiles  []string
	MaxProcs    int64
	ProfFile    string
	KmerLength  uint32
	GobDir      string
	GCFreq      int
	PolyAMax    int
	PolyAMean   int
}

// Parse command line arguments using the flag package.
func (a *CmdArgs) Parse() {
	var targMix string
	var fragMethod string
	var help, version, gob bool
	var kmerLenght int
	var randtest bool

	// Default target mixture:
	const mix_default = " 0.9:(450,50,100,600) + 0.1:(400,2,100,600) "

	// Process simple command line parameters:
	flag.Int64Var(&a.ReqFrags, "n", 0, "Number of requested fragments.")
	flag.StringVar(&targMix, "d", mix_default, "Target mixture.")
	flag.IntVar(&a.PolyAMean, "a", 150, "Mean length of poly-A tail.")
	flag.IntVar(&a.PolyAMax, "amax", 300, "Maximum length length of poly-A tail.")
	flag.Int64Var(&a.NrCycles, "c", 11, "Number of PCR cycles.")
	flag.Float64Var(&a.StrandBias, "b", 0.5, "Strand bias.")
	flag.Float64Var(&a.PrimingBias, "p", 0.0, "Priming bias intentsity parameter.")
	flag.IntVar(&kmerLenght, "k", 6, "Primer length.")
	flag.Float64Var(&a.FixedEff, "e", 0.0, "Fixed per-cyle PCR efficiency.")
	flag.Float64Var(&a.GcEffA, "gca", 8.0, "GC efficiency parameter: a")
	flag.Float64Var(&a.GcEffB, "gcb", 0.8, "GC efficiency parameter: b")
	flag.Float64Var(&a.LenEffPar, "l", 0.001, "Length efficinecy parameter.")
	flag.StringVar(&a.ReportFile, "r", "rlsim_report.tab", "Report file.")
	flag.Int64Var(&a.MaxProcs, "t", 2, "Number of cores to use.")
	flag.StringVar(&fragMethod, "f", "after_prim", "Fragmentation method.")
	flag.BoolVar(&gob, "g", false, "Store fragments on disk.")
	flag.StringVar(&a.GobDir, "gobdir", "", "Directory to store gob files.")
	flag.IntVar(&a.GCFreq, "gcfreq", 0, "Force garbage collection after processing <gcfreq> transcripts.")
	flag.StringVar(&a.ProfFile, "prof", "", "Write out CPU profiling information.")
	flag.BoolVar(&help, "h", false, "Print out help message.")
	flag.BoolVar(&version, "V", false, "Print out version.")
	flag.BoolVar(&a.Verbose, "v", true, "Toggle verbose mode.")
	flag.BoolVar(&randtest, "randt", false, "Generate random numbers for testing.")

	// Redefine usage:
	flag.Usage = func() {
		fmt.Printf("Simulate RNA-seq library preparation with PCR biases and size selection (version: %s).\n\n", VERSION)
		fmt.Printf(`Usage:
        rlsim -n requested fragments [optional arguments] [optional input file] 

Optional arguments:
                argument                    type    default  
        -d      target mixture              string  [check source]
        -f      fragmentation method        string  "after_prim"
        -a      mean poly-A tail length     int     150
        -amax   max poly-A tail lengths     int     300
        -b      strand bias                 float   0.5
        -c      PCR cycles                  int     11
        -p      priming bias parameter      float   0.0
        -k      primer length               int     6
        -e      fixed PCR efficiency        float   0.0
        -gca    GC efficiency parameter: a  float   8.0
        -gcb    GC efficiency parameter: b  float   0.8
        -l      length efficiency parameter float   0.001
        -r      report file                 string  "rlsim_report.tab"
        -t      number of cores to use      int     2
        -g      store fragments on disk     bool    false
        -gobdir fragment directory          string  "rlsim_gob_$PID"
        -v      toggle verbose mode         bool    false
        -h      print usage and exit        bool    false
        -V      print version and exit      bool    false
        -prof   write CPU profiling info    string  "prof.out" 
        -gcfreq trigger garbage collection  int     0
        -randt  generate RNG test files     bool    false

Example:
        rlsim -n 2000000 transcripts.fa
        cat transcripts.fa | rlsim -n 2000000
`)
		os.Exit(0)
	}

	flag.Parse()
	// Print usage:
	if help {
		flag.Usage()
		os.Exit(0)
	}
	// Print version:
	if version {
		fmt.Printf("%s\n", VERSION)
		os.Exit(0)
	}
	// Check primer length:
	if kmerLenght < 1 {
		flag.Usage()
		L.Fatal("Invalid primer length!")
	}
	a.KmerLength = uint32(kmerLenght)

	// Check poly-A tail parameters:
	if a.PolyAMax < 0 {
		L.Fatalf("The maximum length of poly-A tails must be non-negative!")
	}
	if a.PolyAMean < 0 {
		L.Fatalf("The mean length of poly-A tails must be non-negative!")
	}

	// Generate random numbers for testing:
	if randtest {
		RandTest()
		os.Exit(0)
	}

	// Check the number of requested fragments:
	if a.ReqFrags < 1 {
		flag.Usage()
		L.Fatal("No fragments requested, exiting!")
	}

	// Parse target mixture string:
	a.TargetMix = parseTargetMixString(targMix)

	// Parse fragmentation method string:
	a.FragMethod = parseFragMethodString(fragMethod)

	// Set gob directory
	if gob && a.GobDir == "" {
		a.GobDir = "rlsim_gob_" + fmt.Sprintf("%d", os.Getpid())
	}

	// Set input files
	a.InputFiles = flag.Args()

}

func parseTargetMixString(s string) *TargetMix {
	// The target mixture string format is:
	// proportion:(mean, sd, low, high) + ...

	// Construct mix structure, initialize map:
	mix := &TargetMix{}
	mix.Components = make(map[*MixComp]float64)

	// Split mixture string:
	com_strs := strings.Split(s, "+")

	for _, str := range com_strs {
		str = strings.TrimSpace(str)
		if len(str) == 0 {
			L.Fatal("Empty target mixture string!")
		}
		comp, weight := parseTargetMixComponent(str)
		mix.Components[comp] = weight
	}

	return mix
}

func parseTargetMixComponent(s string) (*MixComp, float64) {
	// The target mixture string format is:
	// proportion:(mean, sd, low, high) + ...
	mixComp := &MixComp{}
	var weight float64

	spl1 := strings.Split(s, ":")

	if len(spl1) != 2 {
		L.Fatal("Invalid mixture component: " + s)
	}

	// Parse componenet weight:
	nr, err := fmt.Sscanf(spl1[0], "%e", &weight)
	if err != nil || nr == 0 {
		L.Fatal("Invalid mixture weight: " + spl1[0])
	}

	// Check component string:
	tmp := spl1[1]
	if tmp[0] != '(' || tmp[len(tmp)-1] != ')' {
		L.Fatal("Malformed mixture component string: " + tmp)
	}

	// Trim parantheses:
	tmp = strings.Trim(tmp, "()")
	spl2 := strings.Split(tmp, ",")

	if len(spl2) != 4 {
		L.Fatalf("Invalid mixture string: \"%s\"", string(tmp))
	}

	// Parse mean:
	nr, err = fmt.Sscanf(spl2[0], "%d", &mixComp.Mean)
	if err != nil || nr == 0 {
		L.Fatalf("Invalid mean \"%s\" in mixture: \"%s\"", spl2[0], tmp)
	}

	// Parse sd:
	nr, err = fmt.Sscanf(spl2[1], "%d", &mixComp.Sd)
	if err != nil || nr == 0 {
		L.Fatalf("Invalid sd \"%s\" in mixture: \"%s\"", spl2[1], tmp)
	}

	// Parse lower bound:
	nr, err = fmt.Sscanf(spl2[2], "%d", &mixComp.Low)
	if err != nil || nr == 0 {
		L.Fatalf("Invalid lower bound \"%s\" in mixture: \"%s\"", spl2[2], tmp)
	}

	// Parse upper bound:
	nr, err = fmt.Sscanf(spl2[3], "%d", &mixComp.High)
	if err != nil || nr == 0 {
		L.Fatalf("Invalid upper bound \"%s\" in mixture: \"%s\"", spl2[3], tmp)
	}

	return mixComp, weight
}

func parseFragMethodString(s string) *FragMethod {
	fm := new(FragMethod)

	spl := strings.Split(s, ":")
	err_str := "Malformed fragmentation method string: "

	if len(spl) > 2 || len(spl) < 1 {
		L.Fatal(err_str + s)
	}

	switch spl[0] {
	case "pre_prim", "after_prim", "prim_jump":
		fm.Name = spl[0]
	default:
		L.Fatal(err_str + s)
	}

	if len(spl) == 2 {
		_, err := fmt.Sscanf(spl[1], "%d", &fm.Param)
		if err != nil {
			L.Fatal(err_str + s)
		}
	}

	// pre_prim method needs a parameter no matter what!
	if fm.Name == "pre_prim" && fm.Param == 0 {
		fm.Param = 2000
	}

	return fm
}
