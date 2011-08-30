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
	"rand"
	"math"
	"time"
	"os"
	"fmt"
)

type Rander interface {
	Split() Rander
	Int63() int64
	Int63n(n int64) int64
	Float64() float64
	Float64f(max float64) float64
	NormFloat64() float64
	NormFloat64p(mean float64, sd float64) float64
	Poisson(xm float64) int32
	Binomial(n uint64, p float64) uint64
	TruncNormUint64(mean uint64, sd uint64, low uint64, high uint64) uint64
	SampleIndexUint64(p []uint64) (ind uint64, ok bool)
	SampleIndexFloat64(p []float64) (ind uint64, ok bool)
	ExpFloat64(x float64) float64
	TruncExpInt(mean int, max int) int
}

type RandGen struct {
	rand     *rand.Rand
	source   rand.Source
	Seed     int64
	Suint64  *[]uint64
	Sfloat64 *[]float64
}

func NewRandGen(seed int64) (rg *RandGen) {
	rg = new(RandGen)
	rg.Seed = seed
	rg.source = rand.NewSource(seed)
	rg.rand = rand.New(rg.source)
	tmpUint64 := make([]uint64, 1048576)
	rg.Suint64 = &tmpUint64
	tmpFloat64 := make([]float64, 1048576)
	rg.Sfloat64 = &tmpFloat64
	return
}

func (rg RandGen) Split() Rander {
	seed := rg.Int63()
	return NewRandGen(seed)
}

func (rg RandGen) Int63() int64 {
	return rg.rand.Int63()
}

func (rg RandGen) Int63n(n int64) int64 {
	return rg.rand.Int63n(n)
}

func (rg RandGen) Float64() float64 {
	return rg.rand.Float64()
}

func (rg RandGen) Float64f(max float64) float64 {
	return (rg.rand.Float64() * max)
}

func (rg RandGen) NormFloat64() float64 {
	return rg.rand.NormFloat64()
}

func (rg RandGen) NormFloat64p(mean float64, sd float64) float64 {
	return (rg.rand.NormFloat64()*sd + mean)
}

func Lgamma(x float64) float64 {
	val, _ := math.Lgamma(x)
	return val
}

// Generate Poisson random variates.
// Algorithm from Press, WH; Teukolsky, SA; Vetterling, WT; Flannery, BP:
// Numerical Recipes in C. The Art of Scientific Computing
func (rg RandGen) Poisson(xm float64) int32 {

	var sq, alxm, g float64
	var em, t, y float64

	// Use direct method if mean is small:
	if xm < 12.0 {
		g = math.Exp(-xm)

		em = -1
		t = 1.0

		for {
			em++
			t *= rg.Float64()
			if !(t > g) {
				return int32(em)
			}
		}
	} else {
		// Use rejection method:
		sq = math.Sqrt(2.0 * xm)
		alxm = math.Log(xm)
		g = xm*alxm - Lgamma(xm+1.0)

		for {
		IL:
			for {
				y = math.Tan(math.Pi * rg.Float64())
				em = sq*y + xm
				if !(em < 0.0) {
					break IL
				}
			}
			em = math.Floor(em)
			t = 0.9 * (1.0 + y*y) * math.Exp(em*alxm-Lgamma(em+1.0)-g)
			if !(rg.Float64() > t) {
				return int32(em)
			}
		}
	}
	L.Panic("Your should never reach this point!")
	return 0
}

// Generate Binomila random variates.
// Algorithm after Press, WH; Teukolsky, SA; Vetterling, WT; Flannery, BP:
// Numerical Recipes in C. The Art of Scientific Computing
func (rg RandGen) Binomial(n uint64, pp float64) uint64 {
	var j uint64

	// Transform if pp <= 0.5. Do not forget to transform back!
	p := pp
	if pp <= 0.5 {
		p = 1.0 - pp
	}

	am := float64(n) * p // Mean of the deviate.
	var bnl uint64

	//Use direct method when n is not too large:
	if n < 25 {

		for j = 1; j <= n; j++ {
			if rand.Float64() < p {
				bnl++
			}
		}

	} else {
		// Use rejection method:
		en := float64(n)
		g := Lgamma(en + 1.0)

		pc := 1.0 - p
		plog := math.Log(p)
		pclog := math.Log(pc)

		sq := math.Sqrt(2.0 * am * pc)

		var em, y, t float64
	REJ:
		for {
		LOR:
			for {
				angle := math.Pi * rand.Float64()
				y = math.Tan(angle)
				em = sq*y + am
				// Loop control:
				if em < 0.0 || em >= (en+1.0) {
				} else {
					break LOR
				}
			}

			em = math.Floor(em)

			t = 1.2 * sq * (1.0 + y*y) * math.Exp(g-Lgamma(em+1.0)-Lgamma(en-em+1.0)+em*plog+(en-em)*pclog)

			// Reject:
			if rand.Float64() > t {
			} else {
				break REJ
			}
		}
		bnl = uint64(em)
	} // rejection method

	// Transform bnl if necessary:
	if p != pp {
		bnl = n - bnl
	}
	return bnl
}

func (rg RandGen) TruncNormUint64(mean uint64, sd uint64, low uint64, high uint64) uint64 {
	// FIXME
	h := float64(high)
	l := float64(low)
	if low < 0 {
		L.Fatal("Cannot sample negative lengths!")
	}
	uf := rg.NormFloat64()*float64(sd) + float64(mean)
	for uf < l || uf > h {
		uf = rg.NormFloat64()*float64(sd) + float64(mean)
	}
	return uint64(uf)
}

func (rg RandGen) SampleIndexFloat64(p []float64) (ind uint64, ok bool) {
	// Increase sampler buffer size if necessary:
	size := len(p)
	buff := *rg.Sfloat64
	cacheSize := len(*rg.Sfloat64)
	if size > cacheSize {
		tmp := make([]float64, size-cacheSize)
		buff = append(buff, tmp...)
		*rg.Sfloat64 = buff
	}
	// Reslice sampler buffer:
	csum := buff[:size]
	if len(csum) != len(p) {
		L.Panic("Sampler buffer size mismatch!")
	}

	ind = 0
	ok = false

	csum[0] = p[0]
	for i := 1; i < len(p); i++ {
		csum[i] = csum[i-1] + p[i]
	}

	// Sample index:
	high := csum[len(p)-1]
	if high == 0.0 {
		return
	}

	u := rg.Float64f(high)
	return SearchFloat64(csum, u), true
}

func SearchUints64(a []uint64, x uint64) uint64 {
	var n, i, j uint64
	// Code taken from sort.Search()
	n = uint64(len(a))

	i, j = 0, n
	for i < j {
		h := i + (j-i)/2 // avoid overflow when computing h
		// i ≤ h < j
		if !(a[h] > x) {
			i = h + 1
		} else {
			j = h
		}
	}
	return i
}

func SearchFloat64(a []float64, x float64) uint64 {
	var n, i, j uint64
	// Code taken from sort.Search()
	n = uint64(len(a))

	i, j = 0, n
	for i < j {
		h := i + (j-i)/2 // avoid overflow when computing h
		// i ≤ h < j
		if !(a[h] > x) {
			i = h + 1
		} else {
			j = h
		}
	}
	return i
}

func (rg RandGen) SampleIndexUint64(p []uint64) (ind uint64, ok bool) {
	// Increase sampler buffer size if necessary:
	size := len(p)
	buff := *rg.Suint64
	cacheSize := len(*rg.Suint64)
	if size > cacheSize {
		tmp := make([]uint64, size-cacheSize)
		buff = append(buff, tmp...)
		*rg.Suint64 = buff
	}
	// Reslice sampler buffer:
	csum := buff[:size]
	if len(csum) != len(p) {
		L.Panic("Sampler buffer size mismatch!")
	}

	ind = 0
	ok = false

	csum[0] = p[0]

	for i := 1; i < len(p); i++ {
		csum[i] = csum[i-1] + p[i]
		if csum[i] < csum[i-1] {
			L.Fatal("Integer overflow detected when calculating cumulative sum!")
		}
	}

	// Sample index:
	high := csum[len(p)-1]
	if high == 0 {
		return
	}

	u := uint64(rg.Int63n(int64(high)))
	return SearchUints64(csum, u), true
}

func (rg RandGen) ExpFloat64(x float64) float64 {
	return (rg.rand.ExpFloat64() / x)
}

// Generate random numbers for testing purposes.
func RandTest() {
	time := time.UTC()
	Rg = NewRandGen(time.Seconds())
	const n = 200000

	// Poisson with mean 10:
	f, e := os.Create("poisson_10.tab")
	if e != nil {
		L.Fatalf("Cannot create output file: %s", e.String())
	}
	defer f.Close()
	for i := 0; i < n; i++ {
		u := fmt.Sprintf("%d\n", Rg.Poisson(10))
		f.WriteString(u)
	}

	// Poisson with mean 100:
	f, e = os.Create("poisson_100.tab")
	if e != nil {
		L.Fatalf("Cannot create output file: %s", e.String())
	}
	defer f.Close()
	for i := 0; i < n; i++ {
		u := fmt.Sprintf("%d\n", Rg.Poisson(100))
		f.WriteString(u)
	}

	// Binomial n = 20, p = 0.3:
	f, e = os.Create("binomial_20_0.3.tab")
	if e != nil {
		L.Fatalf("Cannot create output file: %s", e.String())
	}
	defer f.Close()
	for i := 0; i < n; i++ {
		u := fmt.Sprintf("%d\n", Rg.Binomial(20, 0.3))
		f.WriteString(u)
	}

	// Binomial n = 1000, p = 0.51:
	f, e = os.Create("binomial_1000_0.75.tab")
	if e != nil {
		L.Fatalf("Cannot create output file: %s", e.String())
	}
	defer f.Close()
	for i := 0; i < n; i++ {
		u := fmt.Sprintf("%d\n", Rg.Binomial(1000, 0.75))
		f.WriteString(u)
	}

	// Binomial n = 1000, p = 0.75:
	f, e = os.Create("binomial_1000_0.75.tab")
	if e != nil {
		L.Fatalf("Cannot create output file: %s", e.String())
	}
	defer f.Close()
	for i := 0; i < n; i++ {
		u := fmt.Sprintf("%d\n", Rg.Binomial(1000, 0.75))
		f.WriteString(u)
	}

}

func (rg RandGen) TruncExpInt(mean int, max int) int {
	rate := 1.0 / float64(mean)
	var u int
	for {
		u = int(rg.ExpFloat64(rate))
		if u <= max {
			return u
		}
	}
	panic("You should never reach this point!")
	return 0
}
