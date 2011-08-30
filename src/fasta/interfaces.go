package fasta

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
