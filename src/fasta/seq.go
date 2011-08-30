package fasta

import (
	"fmt"
)

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
