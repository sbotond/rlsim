package fasta

import "testing"

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
