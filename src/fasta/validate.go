package fasta

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
			if c == val {
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
