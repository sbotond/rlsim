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

type FragStater interface {
	UpdateAfterFrag(length uint32)
	UpdateAfterPcr(length uint32, count uint64)
	UpdateAfterSampling(length uint32, count uint64)
	UpdateSampled(length uint32, count uint64)
	UpdateMissing(length uint32, count uint64)
	UpdatePolyALen(length uint32)
	ReportFragStats(rep Reporter)
}

type FragStats struct {
	AfterFrag     LenCountMap
	AfterPcr      LenCountMap
	AfterSampling LenCountMap
	Sampled       LenCountMap
	Missing       LenCountMap
	PolyALen      LenCountMap
}

func NewFragStats() (st *FragStats) {
	st = new(FragStats)
	st.AfterFrag = make(LenCountMap)
	st.AfterPcr = make(LenCountMap)
	st.AfterSampling = make(LenCountMap)
	st.Sampled = make(LenCountMap)
	st.Missing = make(LenCountMap)
	st.PolyALen = make(LenCountMap)
	return
}

func (st FragStats) UpdateAfterFrag(length uint32) {
	st.AfterFrag[length]++
}

func (st FragStats) UpdateAfterPcr(length uint32, count uint64) {
	st.AfterPcr[length] += count
	st.AfterSampling[length] += count
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
}
