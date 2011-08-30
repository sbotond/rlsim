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

import "fmt"

type Fragment interface {
	GetName() string
	GetStart() uint32
	GetEnd() uint32
	GetSeq() string
	GetStrand() string
	SetStrand(strand string) Fragment
	GetId() uint64
	SetId(id uint64) Fragment
	GetTranscript() Transcripter
	String() string
}

type Frag struct {
	start  uint32
	end    uint32
	strand string
	id     uint64
	tr     Transcripter
}

func NewFrag(tr Transcripter, start uint32, end uint32) (f *Frag) {
	f = new(Frag)
	f.tr = tr
	f.start = start
	f.end = end
	f.strand = "+"
	return
}

func (f Frag) GetStart() uint32 {
	return f.start
}

func (f Frag) GetEnd() uint32 {
	return f.end
}

func (f Frag) GetSeq() string {
	if f.tr == nil {
		return "<nil>"
	}

	start := f.GetStart()
	end := f.GetEnd()
	// Maximal coordinate:
	final := f.tr.GetLen() - 1
	var seq string
	if f.strand == "+" {
		seq = f.tr.GetSeq()
	} else {
		seq = f.tr.GetRevSeq()
		// Transform coordinate system:
		start, end = final-end, final-start
	}
	return seq[start:end]
}

func (f Frag) GetStrand() string {
	return f.strand
}

func (f Frag) SetStrand(strand string) Fragment {
	if strand != "+" && strand != "-" {
		L.Fatal("Invalid strand!")
	}
	f.strand = strand
	return f
}

func (f Frag) GetId() uint64 {
	return f.id
}

func (f Frag) SetId(id uint64) Fragment {
	f.id = id
	return f
}

func (f Frag) GetTranscript() Transcripter {
	return f.tr
}

func (f Frag) GetName() string {
	if f.tr == nil {
		return "<nil>"
	}
	return f.tr.GetName()
}

func (f Frag) String() string {
	s := fmt.Sprintf(">Frag_%d %s (Strand %s Offset %d -- %d)\n%s", f.GetId(), f.GetName(), f.GetStrand(), f.GetStart(), f.GetEnd(), f.GetSeq())
	return s
}
