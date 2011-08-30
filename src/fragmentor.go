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

import "sort"

type Fragmentor interface {
	Fragment(tr Transcripter, bindProf BindingProfile, polyAend int, st FragStater, rand Rander)
	GetBindingProfile(tr Transcripter) BindingProfile
	JettisonPrimerCache()
}

func NewFragmentor(method *FragMethod, primingBias float64, kmerLength uint32, tg Targeter) Fragmentor {
	var fg Fragmentor
	L.PrintfV("Fragmentation method: \"%s\" with parameter %d.", method.Name, method.Param)
	switch method.Name {
	case "after_prim":
		fg = NewFragAfterPrim(uint32(method.Param), primingBias, kmerLength, tg)
	case "pre_prim":
		fg = NewFragPrePrim(uint32(method.Param), primingBias, kmerLength, tg)
	case "prim_jump":
		fg = NewFragPrimJump(uint32(method.Param), primingBias, kmerLength, tg)
	default:
		L.Fatal("Invalid fragmentation method.")
	}
	return fg
}

// Fragmentation method: after_prim
type FragAfterPrim struct {
	primingBias float64
	target      Targeter
	primer      Primer
	fragParam   uint32
}

func NewFragAfterPrim(fragParam uint32, pbias float64, kmerLength uint32, tg Targeter) *FragAfterPrim {
	fr := new(FragAfterPrim)
	fr.primingBias = pbias
	fr.target = tg
	fr.primer = NewNNthermo(pbias, kmerLength)
	fr.fragParam = uint32(fragParam)
	return fr
}

func (fr FragAfterPrim) GetBindingProfile(tr Transcripter) BindingProfile {
	return fr.primer.GetBindingProfile(tr)
}

func (fr FragAfterPrim) Fragment(tr Transcripter, bindProf BindingProfile, polyAend int, st FragStater, rand Rander) {
	// Sample mixture component:
	d := fr.target.SampleMixComp(rand)

	// Define the length based on polyAend:
	length := polyAend // Shift end to end of poly-A tail
	// Get the boundaries:
	low := uint32(fr.target.GetLow())
	high := uint32(fr.target.GetHigh())
	// Calculate Poisson rate:
	var rate float64
	if fr.fragParam == 0 {
		rate = (float64(length) / float64(d.Mean)) / 2.0 // Double mean.
	} else {
		rate = float64(length) / float64(fr.fragParam)
	}

	// Sample the number of breakpoints:
	nrBreaks := Rg.Poisson(rate) - 1
	if nrBreaks <= 0 {
		nrBreaks = 1
	}

	// 0 and polyAend are predefined breakpoints:
	breaks := make(sort.IntSlice, nrBreaks+2)
	breaks[0] = 0
	breaks[1] = int(length)

	// Sample breakpoints
	for i := int32(0); i < nrBreaks; i++ {
		breaks[i+2] = int(rand.Int63n(int64(length)))
	}

	breaks.Sort()
	// Iterate over fragments:
FRAG:
	for i := 0; i < len(breaks)-2; i++ {
		start := uint32(breaks[i])
		end := uint32(breaks[i+1])
		// filter out too short fragments:
		if end-start < low {
			continue FRAG
		}

		// Simulate priming:
		new_start := fr.primer.SimulatePriming(bindProf, start, end, rand)
		size := (end - new_start)

		//Filter out too long or too short fragments:
		if size < low || size > high {
			continue FRAG
		}

		// Update transcript tables:
		tr.RegisterFragment(size, new_start, end)
		// Update fragment statistics:
		st.UpdateAfterFrag(size)
	}
}

func (fg FragAfterPrim) JettisonPrimerCache() {
	fg.primer.JettisonCache()
}

// Fragmentation method: pre_prim
type FragPrePrim struct {
	primingBias float64
	target      Targeter
	primer      Primer
	fragParam   float64
}

func NewFragPrePrim(fragParam uint32, pbias float64, kmerLength uint32, tg Targeter) *FragPrePrim {
	fr := new(FragPrePrim)
	fr.primingBias = pbias
	fr.target = tg
	fr.primer = NewNNthermo(pbias, kmerLength)
	fr.fragParam = (1.0 / float64(fragParam))
	return fr
}

func (fr FragPrePrim) GetBindingProfile(tr Transcripter) BindingProfile {
	return fr.primer.GetBindingProfile(tr)
}

func (fr FragPrePrim) Fragment(tr Transcripter, bindProf BindingProfile, polyAend int, st FragStater, rand Rander) {
	// Sample mixture component:
	d := fr.target.SampleMixComp(rand)

	// Define length based on polyAend:
	length := polyAend
	// Get the boundaries:
	low := uint32(fr.target.GetLow())
	high := uint32(fr.target.GetHigh())

	// Pick new start by simulating elongation:
	elong := int(rand.ExpFloat64(fr.fragParam))
	if elong > int(length) {
		elong = int(length)
	}
	elong = int(length) - elong

	// Caluclate Poisson rate for fragmentation:
	rate := (float64(length) / float64(d.Mean))

	// Sample the number of breakpoints:
	nrBreaks := Rg.Poisson(rate) - 1
	if nrBreaks <= 0 {
		nrBreaks = 1
	}

	// elong and tr.Len are predefined breakpoints:
	breaks := make(sort.IntSlice, nrBreaks+2)
	breaks[0] = elong
	breaks[1] = int(length)

	// Sample breakpoints
	for i := int32(0); i < nrBreaks; i++ {
		breaks[i+2] = elong + int(rand.Int63n(int64(length)-int64(elong)))
	}

	breaks.Sort()
	// Iterate over fragments:
FRAG:
	for i := 0; i < len(breaks)-2; i++ {
		start := uint32(breaks[i])
		end := uint32(breaks[i+1])
		size := end - start

		//Filter out too long or too short fragments:
		if size < low || size > high {
			continue FRAG
		}

		// Update transcript tables:
		tr.RegisterFragment(size, start, end)
		// Update fragment statistics:
		st.UpdateAfterFrag(size)
	}
}

func (fg FragPrePrim) JettisonPrimerCache() {
	fg.primer.JettisonCache()
}

// Fragmentation method: prim_jump
type FragPrimJump struct {
	primingBias float64
	target      Targeter
	primer      Primer
	fragParam   float64
	kmerLength  uint32
}

func NewFragPrimJump(fragParam uint32, pbias float64, kmerLength uint32, tg Targeter) *FragPrimJump {
	fr := new(FragPrimJump)
	fr.primingBias = pbias
	fr.target = tg
	fr.primer = NewNNthermo(pbias, kmerLength)
	if fragParam != 0 {
		fr.fragParam = (1.0 / float64(fragParam))
	}
	fr.kmerLength = kmerLength
	return fr
}

func (fr FragPrimJump) GetBindingProfile(tr Transcripter) BindingProfile {
	return fr.primer.GetBindingProfile(tr)
}

func (fr FragPrimJump) Fragment(tr Transcripter, bindProf BindingProfile, polyAend int, st FragStater, rand Rander) {
	var start, end, new_end uint32

	// Define final base based on polyAend:
	ltmp := int(tr.GetLen()) - len(bindProf)
	final := uint32(polyAend - ltmp)

	// Get the boundaries:
	low := uint32(fr.target.GetLow())
	high := uint32(fr.target.GetHigh())

FRAG:
	for end < final {
		// Simulate priming:
		start = fr.primer.SimulatePriming(bindProf, end, final, rand)

		// Sample fragment length:
		if fr.fragParam != 0.0 {
			new_end = start + uint32(rand.ExpFloat64(fr.fragParam))
		} else {
			new_end = start + fr.target.SampleMixLen(rand)
		}

		// Check end overflow:
		if new_end > final {
			new_end = final
		}
		// Update end:
		end = new_end
		size := (end - start)

		//Filter out too long or too short fragments:
		if size < low || size > high {
			continue FRAG
		}

		// Update transcript tables:
		tr.RegisterFragment(size, start, end)
		// Update fragment statistics:
		st.UpdateAfterFrag(size)

	}

}

func (fg FragPrimJump) JettisonPrimerCache() {
	fg.primer.JettisonCache()
}
