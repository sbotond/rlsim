package fasta

import (
	"io"
	"os"
)

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
		L.Fatalf("Cannot create fasta file: %s", err.String())
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
