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
	"encoding/json"
	"fmt"
	"os"
)

type Reporter interface {
	ReportMapInt32t64(m map[uint32]uint64, xl string, yl string, title string, vis string)
	ReportSliceInt32t64(x []uint32, y []uint64, xl string, yl string, title string, vis string)
	ReportSliceInt32f64(x []uint32, y []float64, xl string, yl string, title string, vis string)
	ReportSliceFloat64f64(x []float64, y []float64, xl string, yl string, title string, vis string)
	WriteJSON()
}

type JSONEntry struct {
	Xl   string
	Yl   string
	Vis  string
	Pos  int
	Data interface{}
}

type Report struct {
	ReportFile  string
	cursor      *int
	fileHandler *os.File
	Writer      *bufio.Writer
	JSONPool    map[string]JSONEntry
}

func NewReport(File string) (rep Report) {
	f, err := os.Create(File)
	if err != nil {
		err_str := fmt.Sprintf("Could not create report file %s: %s", File, err.Error())
		L.Fatal(err_str)
	}

	rep = Report{}
	rep.ReportFile = File
	rep.fileHandler = f
	rep.Writer = bufio.NewWriter(f)
	rep.JSONPool = make(map[string]JSONEntry)
	rep.cursor = new(int)
	return
}

func (r Report) IncPos() {
	(*r.cursor)++
}

func (r Report) ReportMapInt32t64(m map[uint32]uint64, xl string, yl string, title string, vis string) {
	data := make(map[string]uint64)
	cont := JSONEntry{Xl: xl, Yl: yl, Vis: vis, Pos: *r.cursor}

	for x, y := range m {
		s := fmt.Sprintf("%d", x)
		data[s] = y
	}
	cont.Data = data
	r.JSONPool[title] = cont
	r.IncPos()
}

func (r Report) ReportSliceInt32t64(x []uint32, y []uint64, xl string, yl string, title string, vis string) {
	if len(x) != len(y) {
		L.Fatalf("Cannot report %s: x/y length mismatch!\n", title)
	}
	data := make(map[string]uint64)

	for i := 0; i < len(x); i++ {
		s := fmt.Sprintf("%d", x[i])
		data[s] = y[i]
	}
	cont := JSONEntry{Xl: xl, Yl: yl, Data: data, Vis: vis, Pos: *r.cursor}
	r.JSONPool[title] = cont
	r.IncPos()
}

func (r Report) ReportSliceInt32f64(x []uint32, y []float64, xl string, yl string, title string, vis string) {
	if len(x) != len(y) {
		L.Fatalf("Cannot report %s: x/y length mismatch!\n", title)
	}
	data := make(map[string]float64)

	for i := 0; i < len(x); i++ {
		s := fmt.Sprintf("%d", x[i])
		data[s] = y[i]
	}
	cont := JSONEntry{Xl: xl, Yl: yl, Data: data, Vis: vis, Pos: *r.cursor}
	r.JSONPool[title] = cont
	r.IncPos()
}

func (r Report) ReportSliceFloat64f64(x []float64, y []float64, xl string, yl string, title string, vis string) {
	if len(x) != len(y) {
		L.Fatalf("Cannot report %s: x/y length mismatch!\n", title)
	}
	data := make(map[string]float64)

	for i := 0; i < len(x); i++ {
		s := fmt.Sprintf("%g", x[i])
		data[s] = y[i]
	}
	cont := JSONEntry{Xl: xl, Yl: yl, Data: data, Vis: vis, Pos: *r.cursor}
	r.JSONPool[title] = cont
	r.IncPos()
}

func (r Report) WriteJSON() {
	bytes, err := json.Marshal(r.JSONPool)
	if err != nil {
		L.Fatal("Failed to marshal report!")
	}
	r.Writer.Write(bytes)
	r.Writer.Flush()
	r.fileHandler.Close()
}
