package fasta

import (
	"io"
	"os"
	"strings"
	"bytes"
)

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
		L.Fatalf("Could not open file: %s", err.String())
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
		case os.EOF:
			// Return the last record:
			str := string(r.buff.bytes[:r.buff.pos])
			r.buff.pos = -1
			return NewFasta(str)
		case nil:
			// do nothing
		default:
			L.Fatalf("Error when reading bytes: %s", e.String())
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
