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

import "math"

type Thermocycler interface {
	Pcr(tr Transcripter, p Pooler, st FragStater, rand Rander)
	ReportEffFunctions(tg Targeter, rep Reporter)
}

type Techne struct {
	NrCycles  int64
	FixedEff  float64
	GcEffA    float64
	GcEffB    float64
	LenEffPar float64
}

func NewTechne(NrCycles int64, FixedEff float64, GcEffA float64, GcEffB float64, LenEffPar float64) *Techne {
	tn := new(Techne)
	tn.NrCycles = NrCycles
	tn.FixedEff = FixedEff
	tn.GcEffA = GcEffA
	tn.GcEffB = GcEffB
	tn.LenEffPar = LenEffPar
	return tn
}

func (tn Techne) Pcr(tr Transcripter, p Pooler, st FragStater, rand Rander) {
	frags := tr.GetFragStructs()
	// Iterate over lengths:
	for length, sec := range *frags {
		// Total fragments with current length:
		total := uint64(0)
		// Calculate length efficiency:
		lengthEff := tn.CalcLengthEff(length)
		size := len(sec.Count)
		// Iterate over fragments:
		for i := 0; i < size; i++ {
			// Amplify fragment:
			ampliCount := tn.AmplifyFragment(tr, sec.Start[i], sec.End[i], sec.Count[i], lengthEff, rand)
			// Update fragment count:
			sec.Count[i] = ampliCount
			total += ampliCount
		}
		// Register into pool:
		p.RegisterFragments(tr, length, total)
		// Update AfterPcr stats:
		st.UpdateAfterPcr(length, total)
	}
}

func (tn Techne) AmplifyFragment(tr Transcripter, start uint32, end uint32, icount uint64, lengthE float64, rand Rander) uint64 {
	e := tn.FixedEff
	var oldIcount uint64
	if e == 0.0 {
		e = lengthE * tn.CalcGcEff(tr, start, end)
	}
	for i := int64(0); i < tn.NrCycles; i++ {
		oldIcount = icount
		icount += rand.Binomial(icount, e)
		if icount < oldIcount {
			L.Fatal("Integer overflow detected when amplifying fragments!")
		}
	}
	return icount
}

func (tn Techne) CalcLengthEff(l uint32) (e float64) {
	// eff := l^-tn.LenEffPar
	return math.Pow(float64(l), -tn.LenEffPar)
}

func (tn Techne) CalcGcEff(tr Transcripter, start uint32, end uint32) (e float64) {
	a := tn.GcEffA
	b := tn.GcEffB
	seq := tr.GetSeq()[start:end]
	// Calculate GC content:
	var gc float64
	for i := 0; i < len(seq); i++ {
		if seq[i] == 'G' || seq[i] == 'C' {
			gc++
		}
	}
	gc = gc / float64(len(seq))
	// eff := b + (1 - b) * [ (1-gc^a)^a ]
	eff := b + (1-b)*math.Pow((1-math.Pow(gc, a)), a)
	return eff
}

func (tn Techne) CalcGcEffFixed(gc float64) (e float64) {
	a := tn.GcEffA
	b := tn.GcEffB
	// eff := b + (1 - b) * [ (1-gc^a)^a ]
	eff := b + (1-b)*math.Pow((1-math.Pow(gc, a)), a)
	return eff
}

func (tn Techne) ReportEffFunctions(tg Targeter, rep Reporter) {
	// Report GC efficiency function:
	var x [201]float64
	var y [201]float64
	var gc float64
	gc = 0.0
	for i := 0; i <= 200; i++ {
		x[i] = gc
		y[i] = tn.CalcGcEffFixed(gc / 100.0)
		gc = gc + 0.5
	}
	rep.ReportSliceFloat64f64(x[:], y[:], "GC content (%)", "Efficiency", "GC efficiency function", "line")

	// Report length efficiency function:
	low := uint32(tg.GetLow())
	high := uint32(tg.GetHigh())
	xl := make([]uint32, high-low+1)
	yl := make([]float64, high-low+1)
	for i := uint32(0); i <= (high - low); i++ {
		xl[i] = low + i
		yl[i] = tn.CalcLengthEff(low + i)
	}
	rep.ReportSliceInt32f64(xl, yl, "Length", "Efficiency", "Length efficiency function", "line")
}
