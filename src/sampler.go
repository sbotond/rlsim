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
	"runtime"
	"fmt"
)

type LenCountMap map[uint32]uint64
type TrLenCountMap map[Transcripter]LenCountMap

type Sampler interface {
	SampleFragments(p Pooler, tg Targeter, st FragStater, rand Rander)
	String() string
}

type LenSampler struct {
	StrandBias float64
}

type Request struct {
	Tr  Transcripter
	Len uint32
	Req uint64
}

func NewRequest(tr Transcripter, length uint32, count uint64) *Request {
	req := &Request{tr, length, count}
	return req
}

func NewLenSampler(StrandBias float64) (sl *LenSampler) {
	sl = new(LenSampler)
	sl.StrandBias = StrandBias
	return
}

func (sl LenSampler) SampleTranscripts(p Pooler, tg Targeter, st FragStater, rand Rander) TrLenCountMap {
	m := make(TrLenCountMap)

LEN:
	for {
		length, count, ok := tg.NextLenCount()
		if !ok {
			break LEN
		}

		// Sample transcripts:
		var i uint64 = 0
	TR:
		for i < count {
			tr, okt := p.SampleTranscript(length, rand)
			// Break if we are out of fragments:
			if !okt {
				break TR
			}
			i++

			// Update map:
			lcm, okm := m[tr]
			if !okm {
				lcm = make(LenCountMap)
				m[tr] = lcm
			}
			lcm[length]++
		}
		// Update fragment statistics:
		st.UpdateSampled(length, i)
		st.UpdateAfterSampling(length, i)
		missing := int64(count) - int64(i)
		if missing < 0 {
			L.Panic("Negative value for missing fragments! ", missing)
		}
		st.UpdateMissing(length, uint64(missing))
		// Jettison lengths from pool:
		p.JettisonLenTrCounts(length)
	} // Length iter.
	return m
}

func (sl LenSampler) SampleFragments(p Pooler, tg Targeter, st FragStater, rand Rander) {
	trLenCount := sl.SampleTranscripts(p, tg, st, rand)

	// Calculate the number of sampler gorutines:
	procs := runtime.GOMAXPROCS(0) - 1
	if procs < 1 {
		procs = 1
	}

	// Split random number generator:
	randS := make([]Rander, procs)
	for i := 0; i < procs; i++ {
		randS[i] = rand.Split()
	}

	// Iterate over transcripts:
	var fragCount uint64
	for tr, LenCount := range trLenCount {
		// Ungob fragments:
		tr.Ungob()

		// Make request channel:
		reqChan := make(chan *Request, 500)
		// Make fragment chanel slice:
		fragChans := make([](chan Fragment), procs)

		// Launch length sampler gorutines:
		for i := 0; i < procs; i++ {
			fragChans[i] = sl.GetFragChan(reqChan, randS[i])
		}

		// Send requests:
		go sl.FillReqChan(reqChan, tr, LenCount)

		// Receive and count fragments:
		fragCount += sl.ReceiveFragments(fragChans)

		// Jettison fragment structures:
		tr.JettisonFragStructs()

	}

	L.PrintfV("Sampled %d fragments.\n", fragCount)
	L.PrintfV("Missing fragments: %d\n", uint64(tg.GetReqFrags())-fragCount)
}

func (sl LenSampler) FillReqChan(reqChan chan *Request, tr Transcripter, LenCount LenCountMap) {

	for l, count := range LenCount {
		reqChan <- NewRequest(tr, l, count)
	}

	close(reqChan)
	return
}

func (sl LenSampler) GetFragChan(reqChan chan *Request, rand Rander) (c chan Fragment) {
	// Fork random number generator:
	fragChan := make(chan Fragment, 500)

	go func() {
		for req := range reqChan {
			for i := uint64(0); i < req.Req; i++ {
				frag := req.Tr.SampleFragment(req.Len, rand)
				// Sample fragment strand:
				frag = sl.SampleFragStrand(frag, rand)
				fragChan <- frag
			}
		}
		close(fragChan)
	}()

	return fragChan

	return
}

func (sl LenSampler) ReceiveFragments(fragChans [](chan Fragment)) uint64 {
	nrChans := len(fragChans)
	var fragCount uint64
EVER:
	for {

		closedChans := 0
	CHANS:
		for i := 0; i < nrChans; i++ {
			// Receive fragment:
			frag, Open := <-fragChans[i]
			// Skip if chanel is closed:
			if !Open {
				closedChans++
				continue CHANS
			}
			// Set fragment id:
			frag = frag.SetId(fragCount)
			// Print out fragment:
			fmt.Printf("%s\n", frag.String())
			fragCount++

		} // CHANS
		// Break outer loop if all channels are closed:
		if closedChans == nrChans {
			break EVER
		}

	} // EVER
	return fragCount
}

func (sl LenSampler) SampleFragStrand(f Fragment, rand Rander) Fragment {
	u := rand.Float64()
	if u < sl.StrandBias {
		f = f.SetStrand("-")
	} else {
		f = f.SetStrand("+")
	}
	return f
}

func (sl LenSampler) String() string {
	return "[ parallel length sampler ]"
}
