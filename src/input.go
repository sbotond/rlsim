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

import (
	"os"
	"io"
	"bufio"
	"fasta"
	"strings"
	"fmt"
)

type Seq *fasta.Seq

type Inputer interface {
	NextSeq() Seq
	GetTranscriptChan(tmpDir string, polyAmax int) (c chan *Transcript)
}

type Input struct {
	InputFiles []string
	Funnel     io.Reader
	SeqReader  fasta.SeqReader
}

func NewInput(InputFiles []string) Input {
	buffSize := 1024 * 4 * 10
	res := Input{}

	if len(InputFiles) == 0 {
		res.Funnel = os.Stdin
	} else {
		tmp := make([]io.Reader, 0)
		for _, fileName := range InputFiles {
			file, err := os.Open(fileName)
			if err != nil {
				L.Fatal("Could not open input file " + fileName + ": " + err.String())
			}
			tmp = append(tmp, file)
		}
		res.Funnel = io.MultiReader(tmp...)
	}

	buffReader := bufio.NewReader(res.Funnel)
	var err os.Error
	buffReader, err = bufio.NewReaderSize(buffReader, buffSize)
	if err != nil {
		L.Fatal("Cannot increase read buffer size!")
	}
	fr := fasta.NewFastaReader(buffReader)
	sr := fasta.NewFastaToSeq(fr)
	res.SeqReader = sr

	return res
}

func (input Input) NextSeq() Seq {
	sr := input.SeqReader
	return sr.NextSeq()
}

func (input Input) ParseSeqName(s string) (name string, level uint64, ok bool) {
	err_fmt := "Skipping transcript with malformed name: %s"
	spl := strings.Split(s, "$")
	ok = false
	if len(spl) != 2 || len(spl[0]) == 0 {
		L.Printf(err_fmt, s)
		return
	}
	name = spl[0]

	spl[1] = strings.TrimSpace(spl[1])
	_, err := fmt.Sscanf(spl[1], "%d", &level)
	if err != nil {
		L.Printf(err_fmt, s)
		return
	}

	ok = true
	return
}

func (input Input) GetTranscriptChan(tmpDir string, polyAmax int) (c chan *Transcript) {
	c = make(chan *Transcript, 1000)

	go func() {
		for seq := input.NextSeq(); seq != nil; seq = input.NextSeq() {
			name, level, ok := input.ParseSeqName(seq.Name)
			if !ok {
				continue
			}
			tr := NewTranscript(name, seq.Seq, level, polyAmax, tmpDir)
			c <- tr
		}
		close(c)
	}()
	return
}
