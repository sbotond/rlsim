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
	"math"
)

type RawTarget struct {
	ReqFrags      int64
	LenProbs      *LenProbStruct
	TargetLengths *LenCountStruct
	cursor        *int
	Low           uint64
	High          uint64
	PseudoMixComp *MixComp // Pseudo mixture component.
}

func NewRawTarget(reqFrags int64, lenProbs *LenProbStruct, rand Rander) RawTarget {
	t := RawTarget{}
	t.cursor = new(int)
	t.ReqFrags = reqFrags
	t.LenProbs = lenProbs
	t.TargetLengths = new(LenCountStruct)
	L.PrintfV("Number of requested fragments: %d\n", reqFrags)
	t.SampleLengths(rand)
	L.PrintfV("Finished sampling target lengths from raw size distribution.")
	t.Low, t.High = MinMax(t.TargetLengths.Length)
	t.PseudoMixComp = t.calcMixComp()
	return t
}

func MinMax(d []uint32) (uint64, uint64) {
	if len(d) == 0 {
		return 0, 0
	}
	min := d[0]
	max := d[0]
	for i := 0; i < len(d); i++ {
		if d[i] < min {
			min = d[i]
		}
		if d[i] > max {
			max = d[i]
		}
	}
	return uint64(min), uint64(max)
}

func (tg RawTarget) calcMixComp() *MixComp {
	//Calculate mean fragment length:
	var mean float64
	var m2 float64
	l := tg.TargetLengths.Length
	c := tg.TargetLengths.Count
	t := float64(tg.ReqFrags)
	for i := 0; i < len(l); i++ {
		mean += (float64(c[i]) / t) * float64(l[i])
		m2 += (float64(c[i]) / t) * math.Pow(float64(l[i]), 2)
	}
	sd := math.Sqrt(m2 - math.Pow(mean, 2))
	comp := &MixComp{Type: "raw", Location: mean, Scale: sd, Low: tg.Low, High: tg.High}
	return comp
}

func (tg RawTarget) GetLow() uint64 {
	return tg.Low
}

func (tg RawTarget) GetHigh() uint64 {
	return tg.High
}

func (tg RawTarget) SampleMixComp(rand Rander) *MixComp {
	return tg.PseudoMixComp
}

func (tg RawTarget) SampleLengths(rand Rander) {
	// Sample lengths:
	tmp := make(map[uint32]uint64)
	var i int64
	l := tg.LenProbs.Length
	p := tg.LenProbs.Prob
	for i = 0; i < tg.ReqFrags; i++ {
		pos, ok := rand.SampleIndexFloat64(p)
		if !ok {
			L.Fatal("Error when sampling from raw fragment size distribution!")
		}
		tmp[l[pos]]++
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

func (tg RawTarget) NextLenCount() (length uint32, count uint64, ok bool) {
	i := *tg.cursor
	if *tg.cursor < len(tg.TargetLengths.Length) {
		*tg.cursor++
		return tg.TargetLengths.Length[i], tg.TargetLengths.Count[i], true
	}
	return 0, 0, false
}

func (tg RawTarget) ReportTargetLengths(rep Reporter) {
	rep.ReportSliceInt32t64(tg.TargetLengths.Length, tg.TargetLengths.Count, "Length", "Count", "Target lenghts", "bar")
}

func (tg RawTarget) GetReqFrags() int64 {
	return tg.ReqFrags
}

func (tg RawTarget) SampleMixLen(rand Rander) uint32 {
	pos, ok := rand.SampleIndexFloat64(tg.LenProbs.Prob)
	if !ok {
		L.Fatal("Error when sampling from raw fragment size distribution!")
	}
	return tg.LenProbs.Length[pos]
}
