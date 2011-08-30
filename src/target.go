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
	p := make([]float64, len(mix.Components))
	comps := make([]*MixComp, len(mix.Components))

	i := 0
	for k, v := range mix.Components {
		p[i] = v
		comps[i] = k
		i++
	}

	index, _ := rand.SampleIndexFloat64(p)
	return comps[index]
}

func (tg Target) SampleLengths(rand Rander) {
	// Sample lengths:
	tmp := make(map[uint32]uint64)
	var i int64
	for i = 0; i < tg.ReqFrags; i++ {
		comp := tg.SampleMixComp(rand)
		l := rand.TruncNormUint64(comp.Mean, comp.Sd, comp.Low, comp.High)
		tmp[uint32(l)]++
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
		tmp[length] = 0, false
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
	comp := tg.SampleMixComp(rand)
	return uint32(rand.TruncNormUint64(comp.Mean, comp.Sd, comp.Low, comp.High))
}
