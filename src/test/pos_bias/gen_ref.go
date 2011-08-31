package main

import (
	"flag"
	"fmt"
	"math/rand"
	"strconv"
	"strings"
)

type base_comp struct {
	Nucs   []string
	Probs  []float64
	Cprobs []float64
}

func main() {
	length, level, nr, bc := parse_args()
	_ = bc
	for n := int64(0); n < nr; n++ {
		gen_tr(length, level, n, bc)
	}
}

func gen_tr(length, level, nr int64, bc *base_comp) {
	seq := gen_seq(length, bc)
	fmt.Printf(">t%d$%d\n%s\n", nr, level, seq)
}

func parse_args() (int64, int64, int64, *base_comp) {
	var length int64
	var level int64
	var nr int64
	var p string

	// Process simple command line parameters:
	flag.Int64Var(&length, "l", 10000, "Length of the generated random sequence.")
	flag.Int64Var(&level, "e", 10, "Expression level.")
	flag.Int64Var(&nr, "n", 1000, "Number of transcripts.")
	flag.StringVar(&p, "p", "A:1.0, T:1.0, G:1.0, C:1.0", "Base composition.")
	flag.Parse()
	bc := parse_base_comp(p)
	return length, level, nr, bc
}

func gen_seq(length int64, bc *base_comp) string {
	seq := make([]string, length)
	for i := int64(0); i < length; i++ {
		seq[i] = sample_base(bc)
	}
	return strings.Join(seq, "")
}

func parse_base_comp(s string) *base_comp {
	res := &base_comp{make([]string, 4), make([]float64, 4), make([]float64, 4)}
	units := strings.Split(s, ",")
	for i := 0; i < len(units); i++ {
		u := units[i]
		tmp := strings.Split(u, ":")
		bs, ps := tmp[0], tmp[1]
		bs = strings.TrimSpace(bs)
		ps = strings.TrimSpace(ps)
		res.Nucs[i] = bs
		res.Probs[i], _ = strconv.ParseFloat(ps, 64)
	}
	res.Cprobs = cumsum(res.Probs)
	return res
}

func sample_base(bc *base_comp) string {
	u := rand.Float64()
	var i int
	for i = 0; i < len(bc.Cprobs); i++ {
		if bc.Cprobs[i] > u {
			break
		}
	}
	return bc.Nucs[i]
}

func cumsum(p []float64) []float64 {
	l := len(p)
	c := make([]float64, l)
	c[0] = p[0]
	for i := 1; i < l; i++ {
		c[i] = c[i-1] + p[i]
	}
	head := c[l-1]
	for i := 0; i < l; i++ {
		c[i] = c[i] / head
	}
	return c
}
