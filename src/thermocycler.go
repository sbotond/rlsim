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

import "math"

type Thermocycler interface {
	Pcr(tr Transcripter, p Pooler, st FragStater, rand Rander)
	ReportEffFunctions(tg Targeter, rep Reporter)
}

type Techne struct {
	NrCycles    int64
	FixedEff    float64
	GcEffParam  *EffParam
	LenEffParam *EffParam
	Target      Targeter
	LenScalers  *LenScalers
	RawGcEffs   []float64
	hasRawEffs  bool
}

type LenScalers struct {
	A float64
	B float64
}

func NewTechne(NrCycles int64, FixedEff float64, gcEffParam *EffParam, rawGcEffs []float64, minRawGcEff float64, lenEffParam *EffParam, tg Targeter) *Techne {
	tn := new(Techne)
	tn.NrCycles = NrCycles
	L.PrintfV("Number of PCR cycles: %d", NrCycles)
	tn.FixedEff = FixedEff
	if FixedEff != 0.0 {
		// Got fixed amplification efficiency:
		L.PrintfV("Using fixed amplifcation efficiency: %g", FixedEff)
	} else {

		if rawGcEffs != nil {
			// Got raw gc efficiencies:
			L.PrintfV("Using raw GC efficiencies with a minimum efficiency: %g", minRawGcEff)
			tn.RawGcEffs = tn.ProcessRawGcEffs(rawGcEffs, minRawGcEff)
			tn.hasRawEffs = true
		} else {
			L.PrintfV("GC dependent efficiency parameters: (%g,%g,%g)", gcEffParam.Shape, gcEffParam.Min, gcEffParam.Max)
			L.PrintfV("Length dependent efficiency parameters: (%g,%g,%g)", lenEffParam.Shape, lenEffParam.Min, lenEffParam.Max)
			tn.GcEffParam = gcEffParam
		}

		tn.LenEffParam = lenEffParam
		// Precalculate scaling factors:
		tn.LenScalers = CalcLenScalers(tg, lenEffParam)
	}
	tn.Target = tg
	return tn
}

func CalcLenScalers(tg Targeter, p *EffParam) *LenScalers {
	s := new(LenScalers)
	alpha := math.Pow(float64(tg.GetLow()), -p.Shape)
	beta := math.Pow(float64(tg.GetHigh()), -p.Shape)

	s.A = (p.Max - p.Min) / (alpha - beta)
	s.B = p.Min - s.A*beta

	return s
}

func (tn Techne) ProcessRawGcEffs(effs []float64, minEff float64) []float64 {
	for i := 0; i < len(effs); i++ {
		if effs[i] > 1.0 {
			effs[i] = 1.0
		}
		if effs[i] < minEff {
			effs[i] = minEff
		}
	}
	return effs
}

func (tn Techne) Pcr(tr Transcripter, p Pooler, st FragStater, rand Rander) {
	frags := tr.GetFragStructs()
	// Iterate over lengths:
	for length, sec := range *frags {
		// Total fragments with current length:
		total := uint64(0)
		// Calculate length efficiency:
		var lengthEff float64
		if tn.FixedEff == 0.0 {
			lengthEff = tn.CalcLengthEff(length)
		}
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
	shape := tn.LenEffParam.Shape
	if shape == 0.0 {
		return 1.0
	}
	sc := tn.LenScalers
	return sc.A*math.Pow(float64(l), -shape) + sc.B
}

func (tn Techne) CalcGcEff(tr Transcripter, start uint32, end uint32) (e float64) {
	seq := tr.GetSeq()[start:end]
	// Calculate GC content:
	var gc float64
	for i := 0; i < len(seq); i++ {
		if seq[i] == 'G' || seq[i] == 'C' {
			gc++
		}
	}
	gc = gc / float64(len(seq))

	if tn.hasRawEffs {
		return tn.RawGcEffs[int(gc*100.0)]
	}

	min := tn.GcEffParam.Min
	max := tn.GcEffParam.Max
	shape := tn.GcEffParam.Shape

	eff := min + (max-min)*math.Pow((1-math.Pow(gc, shape)), shape)

	return eff
}

func (tn Techne) CalcGcEffFixed(gc float64) (e float64) {
	if tn.hasRawEffs {
		return tn.RawGcEffs[int(gc*100)]
	}
	min := tn.GcEffParam.Min
	max := tn.GcEffParam.Max
	shape := tn.GcEffParam.Shape

	eff := min + (max-min)*math.Pow((1-math.Pow(gc, shape)), shape)

	return eff
}

func (tn Techne) ReportEffFunctions(tg Targeter, rep Reporter) {
	if tn.FixedEff > 0.0 {
		return
	}
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
