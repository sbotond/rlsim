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
	"fmt"
	"strings"
)

// Mixture component:
type MixComp struct {
	Type     string
	Location float64
	Scale    float64
	Shape    float64
	Low      uint64
	High     uint64
}

func (comp *MixComp) String() (s string) {
	if comp.Type == "sn" {
		s = fmt.Sprintf("%s:(%g, %g, %g, %d, %d)", comp.Type, comp.Location, comp.Scale, comp.Shape, comp.Low, comp.High)
	} else {
		s = fmt.Sprintf("%s:(%g, %g, %d, %d)", comp.Type, comp.Location, comp.Scale, comp.Low, comp.High)
	}
	return
}

func (comp *MixComp) SampleLength(rg Rander) uint32 {
	// Fix the sampled lengths if Low equals High.
	if comp.Low == comp.High {
		return uint32(comp.Low)
	}
	switch comp.Type {
	case "n":
		return rg.TruncNormUint32(comp.Location, comp.Scale, comp.Low, comp.High)
	case "sn":
		return rg.TruncSNormUint32(comp.Location, comp.Scale, comp.Shape, comp.Low, comp.High)
	case "g":
		return rg.TruncGammaUint32(comp.Location, comp.Scale, comp.Low, comp.High)
	default:
		L.Fatalf("Invalid distribution type encountered when sampling length: %s", comp.Type)
	}
	panic("You should never reach this point.")
	return 0
}

type TargetMix struct {
	Components map[*MixComp]float64
}

func (mix *TargetMix) String() (s string) {
	i := 0
	for k, v := range mix.Components {
		s += " " + fmt.Sprintf("\t%.2f", v) + ":" + k.String()
		if i != len(mix.Components)-1 {
			s += "\n"
		}
		i++
	}
	return
}

func parseTargetMixString(s string) *TargetMix {
	// The target mixture string format is:
	// proportion:type:(mean, sd, low, high) + ...

	// Construct mix structure, initialize map:
	mix := &TargetMix{}
	mix.Components = make(map[*MixComp]float64)

	// Split mixture string:
	com_strs := strings.Split(s, "+")

	for _, str := range com_strs {
		str = strings.TrimSpace(str)
		if len(str) == 0 {
			L.Fatal("Empty target mixture string!")
		}
		comp, weight := parseTargetMixComponent(str)
		if comp.High < comp.Low {
			L.Fatal("Maximum fragment size is smaller than minimum fragment size for component \"" + str + "\"!")
		}
		mix.Components[comp] = weight
	}

	return mix
}

func validateMixCompType(s string) {
	valid_types := [...]string{"n", "sn", "g"}
	var found bool
	for _, v := range valid_types {
		if s == v {
			found = true
		}
	}
	if !found {
		L.Fatalf("Invalid mixture type: %s", s)
	}
}

func parseTargetMixComponent(s string) (*MixComp, float64) {
	// The target mixture string format is:
	// proportion:(mean, sd, low, high) + ...
	mixComp := &MixComp{}
	var weight float64
	i := 0

	spl1 := strings.Split(s, ":")

	if len(spl1) != 3 {
		L.Fatal("Something is missing from mixture component: " + s)
	}

	// Parse component weight:
	nr, err := fmt.Sscanf(spl1[i], "%e", &weight)
	if err != nil || nr == 0 {
		L.Fatal("Invalid mixture weight: " + spl1[0])
	}
	i++

	// Parse component type:
	nr, err = fmt.Sscanf(spl1[i], "%s", &mixComp.Type)
	if err != nil || nr == 0 {
		L.Fatal("Invalid mixture type: " + spl1[1])
	}
	validateMixCompType(mixComp.Type)
	i++

	// Check component string:
	tmp := spl1[i]
	if tmp[0] != '(' || tmp[len(tmp)-1] != ')' {
		L.Fatal("Missing paranthesis in mixture component string: " + tmp)
	}

	// Trim parantheses:
	tmp = strings.Trim(tmp, "()")
	spl2 := strings.Split(tmp, ",")
	i = 0

	if mixComp.Type == "sn" {
		if len(spl2) != 5 {
			L.Fatalf("Not enough parameters in skew normal mixture string: \"%s\"", string(tmp))
		}
	} else if len(spl2) != 4 {
		L.Fatalf("Not enough parameters in mixture string: \"%s\"", string(tmp))
	}

	// Parse location:
	nr, err = fmt.Sscanf(spl2[i], "%f", &mixComp.Location)
	if err != nil || nr == 0 {
		L.Fatalf("Invalid location \"%s\" in mixture: \"%s\"", spl2[0], tmp)
	}
	i++

	// Parse scale:
	nr, err = fmt.Sscanf(spl2[i], "%f", &mixComp.Scale)
	if err != nil || nr == 0 {
		L.Fatalf("Invalid scale \"%s\" in mixture: \"%s\"", spl2[1], tmp)
	}
	i++

	// Parse shape:
	if mixComp.Type == "sn" {
		nr, err = fmt.Sscanf(spl2[i], "%f", &mixComp.Shape)
		if err != nil || nr == 0 {
			L.Fatalf("Invalid shape \"%s\" in mixture: \"%s\"", spl2[i], tmp)
		}
		i++
	}

	// Parse lower bound:
	nr, err = fmt.Sscanf(spl2[i], "%d", &mixComp.Low)
	if err != nil || nr == 0 {
		L.Fatalf("Invalid lower bound \"%s\" in mixture: \"%s\"", spl2[i], tmp)
	}
	i++

	// Parse upper bound:
	nr, err = fmt.Sscanf(spl2[i], "%d", &mixComp.High)
	if err != nil || nr == 0 {
		L.Fatalf("Invalid upper bound \"%s\" in mixture: \"%s\"", spl2[i], tmp)
	}

	validateMixComp(mixComp, s)
	return mixComp, weight
}

func validateMixComp(m *MixComp, s string) {
	if m.Location < 0 {
		L.Fatalf("Location parameter is negative in mixture string: %s", s)
	}
	if m.Scale < 0 {
		L.Fatalf("Scale parameter is negative in mixture string: %s", s)
	}

	if m.Low < 0 {
		L.Fatalf("Lower boundary parameter is negative in mixture string: %s", s)
	}
	if m.High < 0 {
		L.Fatalf("Upper boundary parameter is negative in mixture string: %s", s)
	}
}

func (mix *TargetMix) SampleMixComp(rand Rander) *MixComp {
	p := make([]float64, len(mix.Components))
	comps := make([]*MixComp, len(mix.Components))

	i := 0
	for k, v := range mix.Components {
		p[i] = v
		comps[i] = k
		i++
	}

	index, _ := rand.SampleIndexFloat64(p)
	return comps[index]
}

func (mix *TargetMix) SampleMixLen(rand Rander) uint32 {
	comp := mix.SampleMixComp(rand)
	return comp.SampleLength(rand)
}
