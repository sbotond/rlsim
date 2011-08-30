package fasta

import ()

type Fasta struct {
	Record string
}

func NewFasta(str string) *Fasta {
	return &Fasta{str}
}
