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
	"fmt"
	"path"
	"gob"
	"os"
	"strings"
)

type Transcripter interface {
	GetName() string
	GetSeq() string
	SimulatePolyA(polyAmean int, polyAmax int, st FragStater, rand Rander) int
	GetRevSeq() string
	GetExprLevel() uint64
	GetLen() uint32
	RegisterFragment(length uint32, start uint32, end uint32)
	SampleFragment(length uint32, rand Rander) Fragment
	Flatten()
	Fragment(tg Targeter, fg Fragmentor, polyAmean int, polyAmax int, st FragStater, rand Rander)
	Pcr(p Pooler, tc Thermocycler, st FragStater, rand Rander)
	GetFragStructs() *map[uint32]StartEndCountStruct
	String() string
	Cleanup()
	Gob()
	Ungob()
	JettisonFragStructs()
}

type EndCountMap map[uint32]uint64
type StartEndCountMap map[uint32]EndCountMap
type StartEndCountStruct struct {
	Start []uint32
	End   []uint32
	Count []uint64
}

type Transcript struct {
	id   uint64
	name string
	seq  struct {
		plus  string
		minus string
	}
	len         uint32
	exprLevel   uint64
	gobFile     string
	FragMap     *map[uint32]StartEndCountMap
	FragStructs *map[uint32]StartEndCountStruct
}

var maxTranscriptId uint64

func init() {
	maxTranscriptId = 0
}

func NewTranscript(name string, seq string, level uint64, polyAmax int, gobDir string) (tr *Transcript) {
	tr = new(Transcript)
	tr.id = maxTranscriptId
	maxTranscriptId++
	tr.name = name
	tr.seq.plus = seq + strings.Repeat("A", polyAmax)
	tr.seq.minus = RevCompDNA(tr.seq.plus)
	tr.len = uint32(len(tr.seq.plus)) // Length with maximal poly-A tail!
	tr.exprLevel = level
	valFragMap := make(map[uint32]StartEndCountMap, 0)
	tr.FragMap = &valFragMap
	valFragStructs := make(map[uint32]StartEndCountStruct, 0)
	tr.FragStructs = &valFragStructs
	if gobDir != "" {
		tr.gobFile = path.Join(gobDir, fmt.Sprintf("%d_%s.gob", tr.id, tr.name))
	}
	return
}

func (tr Transcript) GetName() string {
	return tr.name
}

func (tr Transcript) GetSeq() string {
	return tr.seq.plus
}

func (tr Transcript) GetRevSeq() string {
	return tr.seq.minus
}

func (tr Transcript) GetLen() uint32 {
	return tr.len
}

func (tr Transcript) GetExprLevel() uint64 {
	return tr.exprLevel
}

func (tr Transcript) Fragment(tg Targeter, fg Fragmentor, polyAmean int, polyAmax int, st FragStater, rand Rander) {
	level := tr.GetExprLevel()
	// Calculate binding profile:
	bindProf := fg.GetBindingProfile(tr)
	var i uint64
	for ; i < level; i++ {
		// Simulate poly-A tail:
		polyAend := tr.SimulatePolyA(polyAmean, polyAmax, st, rand)
		// Fragment transcript:
		fg.Fragment(tr, bindProf, polyAend, st, rand)
	}
}

func (tr Transcript) Pcr(p Pooler, tc Thermocycler, st FragStater, rand Rander) {
	tc.Pcr(&tr, p, st, rand)
}

func (tr Transcript) RegisterFragment(length uint32, start uint32, end uint32) {
	startMap, oks := (*tr.FragMap)[length]
	if !oks {
		startMap = make(StartEndCountMap, 0)
		(*tr.FragMap)[length] = startMap
	}

	endMap, oke := startMap[start]
	if !oke {
		endMap = make(EndCountMap, 0)
		startMap[start] = endMap
	}

	endMap[end]++
}

func (tr Transcript) SampleFragment(length uint32, rand Rander) Fragment {
	s, okt := (*tr.FragStructs)[length]
	if !okt {
		L.Fatalf("Cannot sample fragment: missing length %s!", length)
	}

	index, oks := rand.SampleIndexUint64(s.Count)
	if !oks {
		L.Fatal("Transcript %s is out of fragments! Simulation is inconsistent!", tr.String())
	}
	// Update counts:
	s.Count[index] -= 1

	fg := NewFrag(tr, s.Start[index], s.End[index])
	return fg
}

func (tr Transcript) Flatten() {
	m := *tr.FragStructs
	// Iterate over lengths:
	for length, secm := range *tr.FragMap {

		if length == 0 {
			L.Fatal("Serious trouble: fragment length is zero!")
		}

		t := StartEndCountStruct{}
		t.Start = make([]uint32, 0)
		t.End = make([]uint32, 0)
		t.Count = make([]uint64, 0)

		// Iterate over starts:
		for start, ecm := range secm {
			for end, count := range ecm {
				t.Start = append(t.Start, start)
				t.End = append(t.End, end)
				t.Count = append(t.Count, count)
			}
		}
		m[length] = t
	}
	// Discard FragMap
	*tr.FragMap = nil
}

func (tr Transcript) String() string {
	s := fmt.Sprintf("%s|%d", tr.GetName(), tr.GetExprLevel())
	return s
}

func (tr Transcript) GetFragStructs() *map[uint32]StartEndCountStruct {
	return tr.FragStructs
}

func (tr Transcript) Cleanup() {
	if tr.gobFile == "" {
		return
	}
	err := os.Remove(tr.gobFile)
	if err != nil {
		L.Fatalf("Could not remove gob file \"%s\": %s", tr.gobFile, err.String())
	}
}

func (tr Transcript) Gob() {
	if tr.gobFile == "" {
		return
	}
	f, err := os.Create(tr.gobFile)
	if err != nil {
		L.Fatalf("Could not create gob file \"%s\": %s", tr.gobFile, err.String())
	}
	defer f.Close()
	encoder := gob.NewEncoder(f)
	err = encoder.Encode(*tr.FragStructs)
	if err != nil {
		L.Fatalf("Failed to encode fragemnts for transcript %s: %s", tr.name, err.String())
	}
	// Discard fragment structure:
	*tr.FragStructs = nil
}

func (tr Transcript) Ungob() {
	if tr.gobFile == "" {
		return
	}
	if *tr.FragStructs != nil {
		L.Fatalf("Ungob tries to replace data for transcript %s", tr.name)
	}
	f, err := os.Open(tr.gobFile)
	if err != nil {
		L.Fatalf("Could not open gob file \"%s\": %s", tr.gobFile, err.String())
	}
	defer f.Close()
	val := make(map[uint32]StartEndCountStruct)
	decoder := gob.NewDecoder(f)
	err = decoder.Decode(&val)
	if err != nil {
		L.Fatalf("Failed to decode fragments for transcript %s!", tr.name)
	}
	*tr.FragStructs = val
}

func (tr Transcript) JettisonFragStructs() {
	*tr.FragStructs = nil
}

func (tr Transcript) SimulatePolyA(polyAmean int, polyAmax int, st FragStater, rand Rander) int {
	polyAlength := rand.TruncExpInt(polyAmean, polyAmax)
	st.UpdatePolyALen(uint32(polyAlength))
	// If polyAlength === polyAmax => tr.len is returned
	return (int(tr.len) - polyAmax + polyAlength)
}
