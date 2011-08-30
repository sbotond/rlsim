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

type Primer interface {
	GetBindingProfile(tr Transcripter) BindingProfile
	SimulatePriming(profile BindingProfile, start uint32, end uint32, rand Rander) uint32
	JettisonCache()
	GetKmerLength() uint32
}

type NNthermo struct {
	nn_dG       map[string]float64
	kmerCache   *map[string]float64
	kmerLength  uint32
	initiation  float64
	termAT      float64
	symmetry    float64
	primingBias float64
}

type BindingProfile []float64

func NewNNthermo(primingBias float64, kmerLength uint32) *NNthermo {
	nn := new(NNthermo)
	// Doublet binding enerigies:
	nn.nn_dG = map[string]float64{
		"AA": -1.0,
		"TT": -1.0,
		"AG": -1.28,
		"CT": -1.28,
		"AC": -1.44,
		"GT": -1.44,
		"GA": -1.30,
		"TC": -1.30,
		"GG": -1.84,
		"CC": -1.84,
		"TG": -1.45,
		"CA": -1.45,
		"CG": -2.17,
		"GC": -2.24,
		"AT": -0.88,
		"TA": -0.58,
	}

	// Initialize hexamer cache:
	hexTmp := make(map[string]float64, 1296)
	nn.kmerCache = &hexTmp

	// Intilize other parameters:
	nn.initiation = +1.96
	nn.termAT = +0.05
	nn.symmetry = +0.43
	nn.primingBias = primingBias
	nn.kmerLength = kmerLength
	return nn
}

func (nn NNthermo) KmerAffinity(kmer string) float64 {
	// Check in k-mer cache:
	ce, ok := (*nn.kmerCache)[kmer]
	if ok {
		return ce
	}

	// Initialize:
	kl := len(kmer)
	e := -nn.initiation

	// Sum up doublet energies:
	for i := 0; i < kl-2; i++ {
		doublet := kmer[i : i+2]
		e -= nn.nn_dG[doublet]
	}
	// Terminal A/T penalties:
	if kmer[0] == 'A' || kmer[0] == 'T' {
		e -= nn.termAT
	}
	if kmer[kl-1] == 'A' || kmer[kl-1] == 'T' {
		e -= nn.termAT
	}

	// Symmetry correction:
	sym := true
	for j := 0; j <= int(kl/2); j++ {
		if kmer[j] != kmer[kl-j-1] {
			sym = false
			break
		}
	}
	if sym {
		e -= nn.symmetry
	}

	// Store in cache:
	(*nn.kmerCache)[kmer] = e
	return e
}

func (nn NNthermo) GetBindingProfile(tr Transcripter) BindingProfile {
	seq := tr.GetSeq()
	profLen := tr.GetLen() - nn.kmerLength
	p := make(BindingProfile, profLen)

	for i := uint32(0); i < profLen; i++ {
		end := i + nn.kmerLength
		p[i] = nn.KmerAffinity(seq[i:end]) + nn.primingBias
	}
	return p
}

func (nn NNthermo) SimulatePriming(profile BindingProfile, start uint32, end uint32, rand Rander) uint32 {
	profLen := uint32(len(profile))
	if end > profLen {
		end = profLen
	}
	if start >= end {
		return 0
	}
	p := profile[start:end]
	ind, ok := rand.SampleIndexFloat64(p)
	if !ok {
		L.Fatalf("Error when simulating priming!\n")
	}
	return (uint32(ind) + start)
}

func (nn NNthermo) JettisonCache() {
	*nn.kmerCache = nil
}
func (nn NNthermo) GetKmerLength() uint32 {
	return nn.kmerLength
}
