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

type Primer interface {
	GetBindingProfile(tr Transcripter, reverse bool) BindingProfile
	GetBindingProfiles(tr Transcripter, rev bool) *BindingProfiles
	SimulatePriming(profile BindingProfile, start uint32, end uint32, rand Rander) uint32
	JettisonCache()
	GetKmerLength() uint32
}

type NNthermo struct {
	nn_dS        map[string]float64
	nn_dH        map[string]float64
	kmerCache    *map[string]float64
	kmerLength   uint32
	initiation_H float64
	initiation_S float64
	termAT_H     float64
	termAT_S     float64
	symmetry_S   float64
	T            float64
}

type BindingProfile []float64

type BindingProfiles struct {
	Forward BindingProfile
	Reverse BindingProfile
}

func NewNNthermo(primingTemp float64, kmerLength uint32) *NNthermo {
	nn := new(NNthermo)

	// Doublet binding enerigies from:
	// SantaLucia and Hicks(2004) The thermodynamics of DNA structural motifs - Annu Rev.Biophys. Biomol. Struct. 33:415-440
	nn.nn_dH = map[string]float64{
		"AA": -7.6,
		"TT": -7.6,
		"AG": -7.8,
		"CT": -7.8,
		"AC": -8.4,
		"GT": -8.4,
		"GA": -8.2,
		"TC": -8.2,
		"GG": -8.0,
		"CC": -8.0,
		"TG": -8.5,
		"CA": -8.5,
		"CG": -10.6,
		"GC": -9.8,
		"AT": -7.2,
		"TA": -7.2,
	}

	nn.nn_dS = map[string]float64{
		"AA": -21.3,
		"TT": -21.3,
		"AG": -21.0,
		"CT": -21.0,
		"AC": -22.4,
		"GT": -22.4,
		"GA": -22.2,
		"TC": -22.2,
		"GG": -19.9,
		"CC": -19.9,
		"TG": -22.7,
		"CA": -22.7,
		"CG": -27.2,
		"GC": -24.4,
		"AT": -20.4,
		"TA": -21.3,
	}

	// Initialize kmer cache:
	hexTmp := make(map[string]float64, 1296)
	nn.kmerCache = &hexTmp

	// Initialize other parameters:
	nn.initiation_H = +0.2
	nn.initiation_S = -5.7

	// Terminal AT penalty:
	nn.termAT_H = +2.2
	nn.termAT_S = +6.9

	// Symmetry correction for entropy:
	nn.symmetry_S = -1.4

	nn.T = primingTemp
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
	dH := nn.initiation_H
	dS := nn.initiation_S

	// Sum up doublet energies:
	for i := 0; i <= kl-2; i++ {
		doublet := kmer[i : i+2]
		dH += nn.nn_dH[doublet]
		dS += nn.nn_dS[doublet]
	}
	// Terminal A/T penalties:
	if kmer[0] == 'A' || kmer[0] == 'T' {
		dH += nn.termAT_H
		dS += nn.termAT_S
	}
	if kmer[kl-1] == 'A' || kmer[kl-1] == 'T' {
		dH += nn.termAT_H
		dS += nn.termAT_S
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
		dS += nn.symmetry_S
		// No symmetry correction for enthalpy.
	}

	// Calculate free energy change:
	e := dH - (nn.T*dS)/1000.0
	// Calculate equlibrium constant:
	k := math.Exp(-e / (1.9872 * nn.T))

	// Store in cache:
	(*nn.kmerCache)[kmer] = k
	return k
}

func (nn NNthermo) GetBindingProfiles(tr Transcripter, rev bool) *BindingProfiles {
	bp := new(BindingProfiles)
	bp.Forward = nn.GetBindingProfile(tr, false)
	if rev {
		bp.Reverse = nn.GetBindingProfile(tr, true)
	}
	return bp
}

func (nn NNthermo) GetBindingProfile(tr Transcripter, reverse bool) BindingProfile {
	var seq string
	if !reverse {
		seq = tr.GetSeq()
	} else {
		seq = tr.GetRevSeq()
	}
	profLen := tr.GetLen() - nn.kmerLength
	p := make(BindingProfile, profLen)

	for i := uint32(0); i < profLen; i++ {
		end := i + nn.kmerLength
		p[i] = nn.KmerAffinity(seq[i:end])
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
