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
	"bufio"
	"fmt"
	"io"
	"os"
	"strings"
)

type Inputer interface {
	NextSeq() *Seq
	GetTranscriptChan(tmpDir string, polyAmax int, st FragStater, exprMul float64) (c chan *Transcript)
}

type Input struct {
	InputFiles []string
	Funnel     io.Reader
	SeqReader  SeqReader
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
				L.Fatal("Could not open input file " + fileName + ": " + err.Error())
			}
			tmp = append(tmp, file)
		}
		res.Funnel = io.MultiReader(tmp...)
	}

	buffReader := bufio.NewReader(res.Funnel)
	buffReader = bufio.NewReaderSize(buffReader, buffSize)
	fr := NewFastaReader(buffReader)
	sr := NewFastaToSeq(fr)
	res.SeqReader = sr

	return res
}

func (input Input) NextSeq() *Seq {
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

func (input Input) GetTranscriptChan(tmpDir string, polyAmax int, st FragStater, exprMul float64) (c chan *Transcript) {
	c = make(chan *Transcript, 1000)

	go func() {
		for seq := input.NextSeq(); seq != nil; seq = input.NextSeq() {
			name, level, ok := input.ParseSeqName(seq.Name)
			if !ok {
				continue
			}
			// Apply the global expression level multiplier:
			level = uint64(float64(level) * exprMul)
			// Update transcript length distribution:
			st.UpdateTrLengths(uint32(len(seq.Seq)), level)
			// Update expression level distribution:
			st.UpdateExprLevels(uint32(level))
			tr := NewTranscript(name, seq.Seq, level, polyAmax, tmpDir)
			c <- tr
		}
		close(c)
	}()
	return
}
