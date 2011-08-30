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
	"os"
	"runtime"
)

type Pooler interface {
	InitTranscripts(input Inputer, tg Targeter, fg Fragmentor, tc Thermocycler, st FragStater, GCFreq int, polyAmean int, polyAmax int, rand Rander)
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
			L.Fatalf("Cannot create god directory \"%s\": %s", p.GobDir, err.String())
		}
	}
	return p
}

// Initialize transcripts from input
func (p Pool) InitTranscripts(input Inputer, tg Targeter, fg Fragmentor, tc Thermocycler, st FragStater, GCFreq int, polyAmean int, polyAmax int, rand Rander) {
	c := input.GetTranscriptChan(p.GobDir, polyAmax)

	if p.GobDir != "" {
		L.Printf("Fragments will be cached to %s.", p.GobDir)
	}

	L.PrintfV("Fragmenting transcripts and amplifying fragments:\n")
	var trCount uint64
	for tr := range c {
		L.PrintfV("\t%d\t%s\n", trCount, tr.String())
		// Fragment transcript:
		tr.Fragment(tg, fg, polyAmean, polyAmax, st, rand)
		// Flatten fragment table:
		tr.Flatten()
		// Amplify fragments:
		tr.Pcr(p, tc, st, rand)
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
		p.LenTrCountMap[length] = nil, false
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
	p.LenTrCountStructs[l] = nil, false
}

func (p Pool) Cleanup() {
	if p.GobDir == "" {
		return
	}
	for _, tr := range *p.Transcripts {
		tr.Cleanup()
	}
	err := os.Remove(p.GobDir)
	if err != nil {
		L.Fatalf("Could not remove gob directory \"%s\": %s", p.GobDir, err.String())
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
