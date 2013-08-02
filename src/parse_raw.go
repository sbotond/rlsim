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
	"encoding/json"
	"fmt"
	"io/ioutil"
)

// Functions for parsing raw parameter values:

type LenProbStruct struct {
	Length []uint32
	Prob   []float64
}

type RawParams struct {
	ReqFrags    int64
	NrCycles    int64
	TargetProbs LenProbStruct
	GcEffs      []float64
}

func DecodeRawParams(fname string) *RawParams {
	res := new(RawParams)
	data, err := ioutil.ReadFile(fname)
	if err != nil {
		L.Fatal("Error when reading raw parameter file %s: %s", fname, err.Error())
	}

	var u interface{}
	json.Unmarshal(data, &u)

	d := u.(map[string]interface{})

	//number of fragments
	nr_raw, ok := d["nr_frags"]
	if !ok {
		L.Fatalf("The nr_frags key is missing from the raw parameter file!")
	}
	res.ReqFrags = int64(nr_raw.(float64))

	//number of fragments
	cycles_raw, ok := d["nr_cycles"]
	if !ok {
		L.Fatalf("The nr_cycles key is missing from the raw parameter file!")
	}
	res.NrCycles = int64(cycles_raw.(float64))

	//fragment size distribution
	dist_raw, ok := d["frag_dist"]
	if !ok {
		L.Fatalf("The frag_dist key is missing from the raw parameter file!")
	}
	distMap := dist_raw.(map[string]interface{})
	res.TargetProbs = LenProbStruct{}
	res.TargetProbs.Length = make([]uint32, len(distMap))
	res.TargetProbs.Prob = make([]float64, len(distMap))

	// flatten frag_dist map
	var i int64
	var totalFrags float64
	i = 0
	var length uint32
	for l, count := range distMap {
		n, e := fmt.Sscanf(l, "%d", &length)
		if n == 0 || e != nil {
			L.Fatal("Failed to convert frag_dist key %s: %s", l, err.Error())
		}
		res.TargetProbs.Length[i] = length
		res.TargetProbs.Prob[i] = count.(float64)
		totalFrags += res.TargetProbs.Prob[i]
		i++
	}

	// Normalise counts:
	for i, v := range res.TargetProbs.Prob {
		res.TargetProbs.Prob[i] = v / totalFrags
	}

	//gc efficiencies
	gc_raw, ok := d["gc_eff"]
	if !ok {
		L.Fatalf("The gc_eff key is missing from the raw parameter file!")
	}
	gcMap := gc_raw.(map[string]interface{})
	gcSlice := make([]float64, len(gcMap))

	i = -1 // reset i to an invalid value
	for gc, eff := range gcMap {
		n, e := fmt.Sscanf(gc, "%d", &i)
		if n == 0 || e != nil {
			L.Fatalf("Failed to convert gc_eff key %s: %s", gc, err.Error())
		}
		gcSlice[i] = eff.(float64)
	}
	res.GcEffs = gcSlice

	return res
}
