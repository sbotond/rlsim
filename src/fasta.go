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
	"bytes"
	"fmt"
	"io"
	"os"
	"strings"
	"testing"
)

type Fasta struct {
	Record string
}

func NewFasta(str string) *Fasta {
	return &Fasta{str}
}

type FastaReader interface {
	NextFasta() *Fasta
}

type FastaWriter interface {
	WriteFasta(f *Fasta)
}

type SeqReader interface {
	NextSeq() *Seq
}

type SeqWriter interface {
	WriteSeq(*Seq)
}

type Validator interface {
	NextSeq() *Seq
}

type Seq struct {
	Name string
	Seq  string
}

func NewSeq(name, seq string) *Seq {
	s := &Seq{name, seq}
	return s
}

func (s Seq) String() string {
	return fmt.Sprintf("Name: %s\nSeq: %s\n\n", s.Name, s.Seq)
}

type fastaBuffer struct {
	c     [1]byte
	pos   int
	bytes []byte
}

type ReadFasta struct {
	reader io.Reader
	buff   *fastaBuffer
}

func (r ReadFasta) String() string {
	return "[fasta reader]"
}

func OpenFasta(file string) FastaReader {
	f, err := os.Open(file)
	if err != nil {
		L.Fatalf("Could not open file: %s", err.Error())
	}
	r := NewFastaReader(f)
	return r
}

func NewFastaReader(r io.Reader) FastaReader {
	buff := new(fastaBuffer)
	buff.bytes = make([]byte, 2000)
	fr := ReadFasta{r, buff}
	return fr
}

func (r ReadFasta) NextFasta() (f *Fasta) {

	if r.buff.pos == -1 {
		return nil
	}

	for {
		// Get next character:
		_, e := r.reader.Read(r.buff.c[:])
		switch e {
		case io.EOF:
			if r.buff.pos == 0 {
				return nil
			}
			// Return the last record:
			str := string(r.buff.bytes[:r.buff.pos])
			r.buff.pos = -1
			return NewFasta(str)
		case nil:
			// do nothing
		default:
			L.Fatalf("Error when reading bytes: %s", e.Error())
		}

		// Next record:
		if r.buff.c[0] == byte('>') && r.buff.pos > 1 {
			str := string(r.buff.bytes[:r.buff.pos])
			// Place '>' in buffer:
			r.buff.bytes[0] = '>'
			r.buff.pos = 1
			// Return previous record:
			return NewFasta(str)
		}

		// Grow buffer if necessary:
		if r.buff.pos >= len(r.buff.bytes) {
			tmp := make([]byte, len(r.buff.bytes))
			r.buff.bytes = append(r.buff.bytes, tmp...)
		}

		// Place character in buffer:
		r.buff.bytes[r.buff.pos] = r.buff.c[0]
		r.buff.pos++

	}

	return nil
}

func isSpace(c byte) bool {
	switch c {
	case ' ', '\t', '\r', '\n':
		return true
	}
	return false
}

type FastaToSeq struct {
	fastaReader FastaReader
}

func NewFastaToSeq(fr FastaReader) SeqReader {
	fs := FastaToSeq{fr}
	return fs
}

func (sr FastaToSeq) NextSeq() (s *Seq) {
	frec := sr.fastaReader.NextFasta()
	if frec == nil {
		return nil
	}
	// Find out the end of the first line:
	pos := strings.Index(frec.Record, "\n")
	// Parse name:
	if frec.Record[0] != '>' {
		L.Fatalf("Illegal fasta record:\n%s", frec.Record)
	}
	name := frec.Record[1:pos]
	// Remove newlines from sequence:
	seq := despace([]byte(frec.Record[(pos + 1):]))
	frec.Record = ""
	// Convert to uppercase:
	seq = bytes.ToUpper(seq)
	s = NewSeq(name, string(seq))
	return s
}

func despace(b []byte) []byte {
	tmp := make([]byte, len(b))
	i := 0
	for _, c := range b {
		cc := byte(c)
		if isSpace(cc) {
			continue
		}
		tmp[i] = cc
		i++
	}
	return tmp[:i]
}

type WriteFasta struct {
	writer io.Writer
}

func NewFastaWriter(w io.Writer) *WriteFasta {
	fw := &WriteFasta{w}
	return fw
}

func CreateFasta(file string) *WriteFasta {
	f, err := os.Create(file)
	if err != nil {
		L.Fatalf("Cannot create fasta file: %s", err.Error())
	}
	return NewFastaWriter(f)
}

var nl, start []byte

func (fw WriteFasta) WriteSeq(s *Seq) {
	// Write '>'
	_, err := fw.writer.Write([]byte(start))
	if err != nil {
		L.Fatalf("Error when writing \">\" for sequence:\n%s", s.String())
	}
	// Write name:
	_, err = fw.writer.Write([]byte(s.Name))
	if err != nil {
		L.Fatalf("Error when writing name for sequence:\n%s", s.String())
	}
	// Write first newline:
	_, err = fw.writer.Write(nl)
	if err != nil {
		L.Fatalf("Error when writing first newline for sequence:\n%s", s.String())
	}
	// Write sequence newline:
	_, err = fw.writer.Write([]byte(s.Seq))
	if err != nil {
		L.Fatalf("Error when writing sequence for:\n%s", s.String())
	}
	// Write second newline:
	_, err = fw.writer.Write(nl)
	if err != nil {
		L.Fatalf("Error when writing second newline for sequence:\n%s", s.String())
	}
}

func init() {
	start = []byte{'>'}
	nl = []byte{'\n'}
}

type ValidateSeq struct {
	sr   SeqReader
	vtab vTable
	vset string
}

type vTable []int
type vIndex map[string]vTable

var vindex vIndex

func NewValidator(sr SeqReader, set string) Validator {
	v := new(ValidateSeq)
	v.sr = sr
	vt, ok := vindex[set]
	if !ok {
		L.Fatalf("Invalid validation set %s!", set)
	}
	v.vtab = vt
	v.vset = set
	return v
}

func (v ValidateSeq) NextSeq() (s *Seq) {
NEXT_SEQ:
	s = v.sr.NextSeq()
	if s == nil {
		return nil
	}

	for _, c := range s.Seq {
		trouble := true
	BASES:
		for _, val := range v.vtab {
			if c == rune(val) {
				trouble = false
				break BASES
			}
		}
		if trouble {
			L.Printf("[%s] Invalid base \"%s\" in sequence \"%s\"\n", v.vset, string(c), s.Name)
			goto NEXT_SEQ
		}
	}
	return s
}

func init() {
	vindex = make(vIndex, 1)
	vindex["DNA_strict"] = []int{'A', 'T', 'G', 'C'}
	vindex["DNA_ambig"] = []int{'A', 'T', 'G', 'C', 'R', 'Y', 'K', 'M', 'S', 'W', 'B', 'D', 'H', 'V', 'N'}
}

func TestOpen(t *testing.T) {
	t.Log("Testing fasta open")
	fas := OpenFasta("test.fas")
	_ = fas
}

func TestFastaReader(t *testing.T) {
	t.Log("Testing Fasta reader")
	fr := OpenFasta("test.fas")
	for f := fr.NextFasta(); f != nil; f = fr.NextFasta() {
		//fmt.Println(f.Record)
	}
}

func TestSeqReader(t *testing.T) {
	t.Log("Testing Seq reader")
	fr := OpenFasta("test.fas")
	sr := NewFastaToSeq(fr)
	for s := sr.NextSeq(); s != nil; s = sr.NextSeq() {
		//fmt.Print(s)
	}
}

func TestFastaWriter(t *testing.T) {
	t.Log("Testing Fasta writer")
	fr := OpenFasta("test.fas")
	fw := CreateFasta("test_out.fas")
	sr := NewFastaToSeq(fr)
	for s := sr.NextSeq(); s != nil; s = sr.NextSeq() {
		fw.WriteSeq(s)
	}
}

func TestValidator(t *testing.T) {
	t.Log("Testing Fasta writer")
	fr := OpenFasta("test.fas")
	fw := CreateFasta("test_out.fas")
	sr := NewValidator(NewFastaToSeq(fr), "DNA_strict")
	for s := sr.NextSeq(); s != nil; s = sr.NextSeq() {
		fw.WriteSeq(s)
	}
}
