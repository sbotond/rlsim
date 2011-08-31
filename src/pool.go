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
)

type Pooler interface {
	InitTranscripts(input Inputer, tg Targeter, fg Fragmentor, tc Thermocycler, st FragStater, GCFreq int, polyAParam *TargetMix, exprMul float64, rand Rander, pcrRand Rander)
	AddTranscript(tr Transcripter)
	GetGobDir() string
	RegisterFragments(tr Transcripter, length uint32, count uint64)
	GetTranscripts() []Transcripter
	GetNrTranscripts() int
	Flatten()
	SampleTranscript(length uint32, rand Rander) (Transcripter, bool)
	String() string
	JettisonLenTrCounts(l uint32)
	Cleanup()
}

type TrCountMap map[Transcripter]uint64
type LenTrCountMap map[uint32]TrCountMap

type TrCountStruct struct {
	tr    []Transcripter
	count []uint64
}

type Pool struct {
	Transcripts       *[]Transcripter
	LenTrCountMap     LenTrCountMap
	LenTrCountStructs map[uint32]*TrCountStruct
	GobDir            string
}

// Pool constructor:
func NewPool(gobDir string) *Pool {
	p := new(Pool)
	tmp := make([]Transcripter, 0)
	p.Transcripts = &tmp
	p.LenTrCountMap = make(LenTrCountMap)
	p.LenTrCountStructs = make(map[uint32]*TrCountStruct)
	p.GobDir = gobDir
	if p.GobDir != "" {
		err := os.Mkdir(p.GobDir, 0700)
		if err != nil {
			L.Fatalf("Cannot create gob directory \"%s\": %s", p.GobDir, err.Error())
		}
	}
	return p
}

// Initialize transcripts from input
func (p Pool) InitTranscripts(input Inputer, tg Targeter, fg Fragmentor, tc Thermocycler, st FragStater, GCFreq int, polyAParam *TargetMix, exprMul float64, rand Rander, pcrRand Rander) {
	_, polyAmax := getGlobalMinMax(polyAParam)
	L.PrintfV("Poly(A) tail distribution components:")
	for _, l := range StringToSlice(polyAParam.String()) {
		L.PrintfV("%s\n", l)
	}
	c := input.GetTranscriptChan(p.GobDir, int(polyAmax), st, exprMul)

	if p.GobDir != "" {
		L.PrintfV("Fragments will be cached to %s.", p.GobDir)
	}

	L.PrintfV("Fragmenting transcripts and amplifying fragments:\n")
	var trCount uint64

	// Collect garbage before starting fragmentation:
	runtime.GC()

	for tr := range c {
		L.PrintfV("\t%d\t%s\n", trCount, tr.String())
		// Skip if not expressed:
		if tr.GetExprLevel() == 0 {
			trCount++
			continue
		}
		// Fragment transcript:
		tr.Fragment(tg, fg, polyAParam, int(polyAmax), st, rand) // uses initial seed
		// Flatten fragment table:
		tr.Flatten()
		// Amplify fragments:
		tr.Pcr(p, tc, st, pcrRand) // uses initial or PCR seed
		// Add transcript to pool:
		p.AddTranscript(tr)
		// Store fragments:
		tr.Gob()
		// Explicitly trigger GC:
		if (GCFreq > 0) && ((int(trCount) % GCFreq) == 0) {
			runtime.GC()
		}
		// Increase transcript counter:
		trCount++
	}

	// Jettison primer cache:
	fg.JettisonPrimerCache()
	L.PrintfV("Initialized %d transcripts.\n", p.GetNrTranscripts())
	// Flatten pool:
	p.Flatten()
	// Explicitly trigger GC:
	if GCFreq > 0 {
		runtime.GC()
	}
}

func (p Pool) GetGobDir() string {
	return p.GobDir
}

func (p Pool) GetTranscripts() []Transcripter {
	return *p.Transcripts
}

func (p Pool) GetNrTranscripts() int {
	return len(*p.Transcripts)
}

// Add a transcript to pool:
func (p Pool) AddTranscript(tr Transcripter) {
	tmp := p.Transcripts
	*tmp = append(*p.Transcripts, tr)
}

// Register fragments into LenTrCountMap:
func (p Pool) RegisterFragments(tr Transcripter, length uint32, count uint64) {
	trc, ok := p.LenTrCountMap[length]
	// New length:
	if !ok {
		trc = make(TrCountMap)
		trc[tr] = count
		p.LenTrCountMap[length] = trc
	} else {
		trc[tr] = count
	}

}

// Linearize LenTrCountMap:
func (p Pool) Flatten() {
	m := p.LenTrCountStructs
	// Iterate over lengths:
	for length, trcm := range p.LenTrCountMap {
		// Linearize transcript/count map into slices:
		size := len(trcm)
		trS := make([]Transcripter, size)
		countS := make([]uint64, size)
		i := 0
		for tr, count := range trcm {
			trS[i] = tr
			countS[i] = count
			i++
		}
		m[length] = &TrCountStruct{trS, countS}
		// Delete transcript/count map:
		delete(p.LenTrCountMap, length)
	}

}

func (p Pool) SampleTranscript(length uint32, rand Rander) (Transcripter, bool) {
	tmp, ok := p.LenTrCountStructs[length]
	if !ok {
		return nil, false
	}
	trS := tmp.tr
	countS := tmp.count
	index, oks := rand.SampleIndexUint64(countS)
	if !oks {
		return nil, false
	}
	countS[index]--
	return trS[index], true
}

func (p Pool) JettisonLenTrCounts(l uint32) {
	delete(p.LenTrCountStructs, l)
}

func (p Pool) Cleanup() {
	if p.GobDir == "" {
		return
	}
	for _, tr := range *p.Transcripts {
		tr.Cleanup()
	}
	err := os.RemoveAll(p.GobDir)
	if err != nil {
		L.Fatalf("Could not remove gob directory \"%s\": %s", p.GobDir, err.Error())
	}
}

func (p Pool) String() string {
	s := "[ "
	for _, tr := range *p.Transcripts {
		s += tr.String() + " "
	}
	s += "]"
	return s
}
