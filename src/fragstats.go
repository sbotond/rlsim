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

type FragStater interface {
	UpdateTrLengths(length uint32, count uint64)
	UpdateExprLevels(level uint32)
	UpdateAfterFrag(length uint32)
	UpdateAfterPcr(length uint32, count uint64)
	UpdateAfterSampling(length uint32, count uint64)
	UpdateSampled(length uint32, count uint64)
	UpdateNrSampled(count uint64)
	UpdateMissing(length uint32, count uint64)
	UpdatePolyALen(length uint32)
	ReportFragStats(rep Reporter)
	LogSamplingRatio(sampled uint64)
}

type FragStats struct {
	ExprLevels    LenCountMap
	TrLengths     LenCountMap
	AfterFrag     LenCountMap
	AfterPcr      LenCountMap
	AfterSampling LenCountMap
	Sampled       LenCountMap
	Missing       LenCountMap
	PolyALen      LenCountMap
	TotalFrags    *uint64
	NrSampled     *uint64
}

func NewFragStats() (st *FragStats) {
	st = new(FragStats)
	st.ExprLevels = make(LenCountMap)
	st.TrLengths = make(LenCountMap)
	st.AfterFrag = make(LenCountMap)
	st.AfterPcr = make(LenCountMap)
	st.AfterSampling = make(LenCountMap)
	st.Sampled = make(LenCountMap)
	st.Missing = make(LenCountMap)
	st.PolyALen = make(LenCountMap)
	st.TotalFrags = new(uint64)
	st.NrSampled = new(uint64)
	return
}

func (st FragStats) UpdateTrLengths(length uint32, count uint64) {
	st.TrLengths[length] += count
}
func (st FragStats) UpdateExprLevels(level uint32) {
	st.ExprLevels[level]++
}

func (st FragStats) UpdateAfterFrag(length uint32) {
	st.AfterFrag[length]++
}

func (st FragStats) UpdateAfterPcr(length uint32, count uint64) {
	st.AfterPcr[length] += count
	st.AfterSampling[length] += count

	// Safely update the total fragment count:
	var newTotal uint64
	newTotal = *st.TotalFrags + count
	if newTotal < *st.TotalFrags {
		L.Fatalf("Integer overflow when updating total fragment count!")
	}
	*st.TotalFrags = newTotal
}

func (st FragStats) UpdateAfterSampling(length uint32, count uint64) {
	if count > st.AfterSampling[length] {
		L.Panic("BIG trouble", count, st.AfterSampling[length])
	}
	st.AfterSampling[length] -= count
}

func (st FragStats) UpdateSampled(length uint32, count uint64) {
	st.Sampled[length] += count
}

func (st FragStats) UpdateNrSampled(count uint64) {
	*st.NrSampled = count
}

func (st FragStats) UpdateMissing(length uint32, count uint64) {
	st.Missing[length] += count
}

func (st FragStats) UpdatePolyALen(length uint32) {
	st.PolyALen[length]++
}

func (st FragStats) ReportFragStats(r Reporter) {
	r.ReportMapInt32t64(st.AfterFrag, "Length", "Count", "Fragdist after fragmentation", "bar")
	r.ReportMapInt32t64(st.AfterPcr, "Length", "Count", "Fragdist after PCR", "bar")
	r.ReportMapInt32t64(st.AfterSampling, "Length", "Count", "Fragdist after sampling", "bar")
	r.ReportMapInt32t64(st.Missing, "Length", "Count", "Missing fragments", "bar")
	r.ReportMapInt32t64(st.Sampled, "Length", "Count", "Sampled fragments", "bar")
	r.ReportMapInt32t64(st.PolyALen, "Length", "Count", "Poly-A tail lengths", "bar")
	r.ReportMapInt32t64(st.TrLengths, "Length", "Count", "Transcript lengths", "hist")
	r.ReportMapInt32t64(st.ExprLevels, "Expression level", "Count", "Expression levels", "hist")
	st.ReportSamplingRatio(r)
}

func (st FragStats) LogSamplingRatio(sampled uint64) {
	L.PrintfV("Total number of fragments in the pool: %g", float64(*st.TotalFrags))
	L.PrintfV("Sampling ratio: %g", float64(sampled)/float64(*st.TotalFrags))
}

func (st FragStats) ReportSamplingRatio(r Reporter) {
	x := [...]float64{0}
	y := [...]float64{float64(*st.NrSampled) / float64(*st.TotalFrags)}
	r.ReportSliceFloat64f64(x[:], y[:], "", "Sampling ratio", "Sampling ratio", "bar")
}
