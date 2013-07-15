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

import "sort"

type Fragmentor interface {
	Fragment(tr Transcripter, bindProfs *BindingProfiles, polyAend int, st FragStater, rand Rander)
	GetBindingProfiles(tr Transcripter) *BindingProfiles
	JettisonPrimerCache()
}

func NewFragmentor(method *FragMethod, primingTemp float64, kmerLength uint32, tg Targeter, fragLossProb float64) Fragmentor {
	var fg Fragmentor
	fragFilter := NewFragFilter(fragLossProb)
	L.PrintfV("Fragmentation method: \"%s\" with parameter %d.", method.Name, method.Param)
	switch method.Name {
	case "after_prim":
		fg = NewFragAfterPrim(uint32(method.Param), primingTemp, kmerLength, tg, fragFilter, true, false)
	case "after_prim_double":
		fg = NewFragAfterPrim(uint32(method.Param), primingTemp, kmerLength, tg, fragFilter, true, true)
	case "after_noprim":
		fg = NewFragAfterPrim(uint32(method.Param), primingTemp, kmerLength, tg, fragFilter, false, false)
	case "after_noprim_double":
		fg = NewFragAfterPrim(uint32(method.Param), primingTemp, kmerLength, tg, fragFilter, false, true)
	case "pre_prim":
		fg = NewFragPrePrim(uint32(method.Param), primingTemp, kmerLength, tg, fragFilter)
	case "prim_jump":
		fg = NewFragPrimJump(uint32(method.Param), primingTemp, kmerLength, tg, fragFilter)
	default:
		L.Fatal("Invalid fragmentation method.")
	}
	return fg
}

// Fragmentation method: after_prim
type FragAfterPrim struct {
	primingTemp float64
	simPriming  bool
	doublePrime bool
	target      Targeter
	primer      Primer
	fragParam   uint32
	fragFilter  FragFilter
}

func NewFragAfterPrim(fragParam uint32, pbias float64, kmerLength uint32, tg Targeter, fragFilter FragFilter, simPriming bool, doublePrime bool) *FragAfterPrim {
	fr := new(FragAfterPrim)
	fr.simPriming = simPriming
	fr.doublePrime = doublePrime
	fr.primingTemp = pbias
	fr.target = tg
	fr.primer = NewNNthermo(pbias, kmerLength)
	fr.fragParam = uint32(fragParam)
	fr.fragFilter = fragFilter
	return fr
}

func (fr FragAfterPrim) GetBindingProfiles(tr Transcripter) *BindingProfiles {
	if !fr.simPriming {
		return nil
	}
	if fr.doublePrime {
		return fr.primer.GetBindingProfiles(tr, true)
	}
	return fr.primer.GetBindingProfiles(tr, false)
}

func (fr FragAfterPrim) Fragment(tr Transcripter, bindProfs *BindingProfiles, polyAend int, st FragStater, rand Rander) {
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
		rate_scaler := 2.0
		rate = (float64(length) / float64(d.Location)) / rate_scaler // Multiply produced mean fragment length.
	} else {
		rate = float64(length) / float64(fr.fragParam)
	}

	// Sample the number of breakpoints:
	nrBreaks := rand.Poisson(rate) - 1
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
		var new_start uint32
		if fr.simPriming {
			// Use binding energies to select fragment start:
			new_start = fr.primer.SimulatePriming(bindProfs.Forward, start, end, rand)
			// Select end by priming simulation as well:
			if fr.doublePrime {
				final := uint32(tr.GetLen() - 1)                   // needed for coordinate transformation.
				s, e := (final - end + 1), (final - new_start + 1) // "reverse complement" coordinates.
				ns := fr.primer.SimulatePriming(bindProfs.Reverse, s, e, rand)
				end = final - ns + 1 // "reverse complement" coordinates.
			}
		} else {
			// Do not use binding energies:
			new_start = start + uint32(rand.Int31n(int32(end-start)))
			if fr.doublePrime {
				end = new_start + uint32(rand.Int31n(int32(end-new_start)))
			}
		}
		size := (end - new_start)

		//Filter out too long or too short fragments:
		if size < low || size > high {
			continue FRAG
		}

		// Simulate fragment loss:
		if fr.fragFilter.Filter(rand) == false {
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
	primingTemp float64
	target      Targeter
	primer      Primer
	fragParam   float64
	fragFilter  FragFilter
}

func NewFragPrePrim(fragParam uint32, pbias float64, kmerLength uint32, tg Targeter, fragFilter FragFilter) *FragPrePrim {
	fr := new(FragPrePrim)
	fr.primingTemp = pbias
	fr.target = tg
	fr.primer = NewNNthermo(pbias, kmerLength)
	fr.fragParam = (1.0 / float64(fragParam))
	fr.fragFilter = fragFilter
	return fr
}

func (fr FragPrePrim) GetBindingProfiles(tr Transcripter) *BindingProfiles {
	return nil
}

func (fr FragPrePrim) Fragment(tr Transcripter, bindProfs *BindingProfiles, polyAend int, st FragStater, rand Rander) {
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
	rate := (float64(length) / float64(d.Location))

	// Sample the number of breakpoints:
	nrBreaks := rand.Poisson(rate) - 1
	if nrBreaks <= 0 {
		nrBreaks = 1
	}

	// elong and tr.Len are predefined breakpoints:
	breaks := make(sort.IntSlice, nrBreaks+2)
	breaks[0] = elong
	breaks[1] = int(length)

	// Sample breakpoints
	for i := int32(0); i < nrBreaks; i++ {
		len_tmp := int64(length) - int64(elong)
		// Return if priming event is too close to start:
		if len_tmp < 1 {
			return
		}
		breaks[i+2] = elong + int(rand.Int63n(len_tmp))
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

		// Simulate fragment loss:
		if fr.fragFilter.Filter(rand) == false {
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
	primingTemp float64
	target      Targeter
	primer      Primer
	fragParam   float64
	kmerLength  uint32
	fragFilter  FragFilter
}

func NewFragPrimJump(fragParam uint32, pbias float64, kmerLength uint32, tg Targeter, fragFilter FragFilter) *FragPrimJump {
	fr := new(FragPrimJump)
	fr.primingTemp = pbias
	fr.target = tg
	fr.primer = NewNNthermo(pbias, kmerLength)
	if fragParam != 0 {
		fr.fragParam = (1.0 / float64(fragParam))
	}
	fr.kmerLength = kmerLength
	fr.fragFilter = fragFilter
	return fr
}

func (fr FragPrimJump) GetBindingProfiles(tr Transcripter) *BindingProfiles {
	return fr.primer.GetBindingProfiles(tr, false)
}

func (fr FragPrimJump) Fragment(tr Transcripter, bindProfs *BindingProfiles, polyAend int, st FragStater, rand Rander) {
	var start, end, new_end uint32

	// Define final base based on polyAend:
	ltmp := int(tr.GetLen()) - len(bindProfs.Forward)
	final := uint32(polyAend - ltmp)

	// Get the boundaries:
	low := uint32(fr.target.GetLow())
	high := uint32(fr.target.GetHigh())

FRAG:
	for end < final {
		// Simulate priming:
		start = fr.primer.SimulatePriming(bindProfs.Forward, end, final, rand)

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

		// Simulate fragment loss:
		if fr.fragFilter.Filter(rand) == false {
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
