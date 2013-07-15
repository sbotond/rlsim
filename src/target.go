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

type Targeter interface {
	SampleMixComp(rand Rander) *MixComp
	NextLenCount() (length uint32, count uint64, ok bool)
	ReportTargetLengths(rep Reporter)
	GetReqFrags() int64
	GetLow() uint64
	GetHigh() uint64
	SampleMixLen(rand Rander) uint32
}

type LenCountStruct struct {
	Length []uint32
	Count  []uint64
}

type Target struct {
	ReqFrags      int64
	TargetMix     *TargetMix
	TargetLengths *LenCountStruct
	cursor        *int
	Low           uint64
	High          uint64
}

func NewTarget(ReqFrags int64, TargetMix *TargetMix, rand Rander) Target {
	tg := Target{}
	tg.cursor = new(int)
	tg.TargetLengths = new(LenCountStruct)
	tg.ReqFrags = ReqFrags
	L.PrintfV("Number of requested fragments: %d\n", ReqFrags)
	tg.TargetMix = TargetMix
	L.PrintfV("Target fragment length distribution components:")
	for _, l := range StringToSlice(tg.TargetMix.String()) {
		L.PrintfV("%s\n", l)
	}
	tg.SampleLengths(rand)
	L.PrintfV("Finished sampling target lengths.")
	tg.Low, tg.High = getGlobalMinMax(TargetMix)
	return tg
}

func getGlobalMinMax(m *TargetMix) (min, max uint64) {
	i := 0
	for comp, _ := range m.Components {
		if i == 0 {
			min = comp.Low
		}
		if max < comp.High {
			max = comp.High
		}
		if min > comp.Low {
			min = comp.Low
		}
		i++
	}
	return
}

func (tg Target) GetLow() uint64 {
	return tg.Low
}

func (tg Target) GetHigh() uint64 {
	return tg.High
}

func (tg Target) SampleMixComp(rand Rander) *MixComp {
	mix := tg.TargetMix
	return mix.SampleMixComp(rand)
}

func (tg Target) SampleLengths(rand Rander) {
	// Sample lengths:
	tmp := make(map[uint32]uint64)
	var i int64
	for i = 0; i < tg.ReqFrags; i++ {
		comp := tg.SampleMixComp(rand)
		l := comp.SampleLength(rand)
		tmp[l]++
	}
	// Flatten map:
	lcs := tg.TargetLengths
	size := len(tmp)
	lcs.Length = make([]uint32, size)
	lcs.Count = make([]uint64, size)
	i = 0
	for length, count := range tmp {
		lcs.Length[i] = length
		lcs.Count[i] = count
		delete(tmp, length)
		i++
	}
}

func (tg Target) NextLenCount() (length uint32, count uint64, ok bool) {
	i := *tg.cursor
	if *tg.cursor < len(tg.TargetLengths.Length) {
		*tg.cursor++
		return tg.TargetLengths.Length[i], tg.TargetLengths.Count[i], true
	}
	return 0, 0, false
}

func (tg Target) ReportTargetLengths(rep Reporter) {
	rep.ReportSliceInt32t64(tg.TargetLengths.Length, tg.TargetLengths.Count, "Length", "Count", "Target lenghts", "bar")
}

func (tg Target) GetReqFrags() int64 {
	return tg.ReqFrags
}

func (tg Target) SampleMixLen(rand Rander) uint32 {
	comp := tg.TargetMix.SampleMixComp(rand)
	return comp.SampleLength(rand)
}
