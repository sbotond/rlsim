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

import "os"
import "bufio"
import "fmt"

type Reporter interface {
	ReportMapInt32t64(m map[uint32]uint64, xl string, yl string, title string, vis string)
	ReportSliceInt32t64(x []uint32, y []uint64, xl string, yl string, title string, vis string)
	ReportSliceInt32f64(x []uint32, y []float64, xl string, yl string, title string, vis string)
	ReportSliceFloat64f64(x []float64, y []float64, xl string, yl string, title string, vis string)
}

type Report struct {
	ReportFile string
	Writer     *bufio.Writer
}

func NewReport(File string) (rep Report) {
	f, err := os.Create(File)
	if err != nil {
		err_str := fmt.Sprintf("Could not create report file %s: %s", File, err.String())
		_ = err_str
		L.Fatal(err_str)
	}

	rep = Report{}
	rep.ReportFile = File
	rep.Writer = bufio.NewWriter(f)
	return
}

func (r Report) ReportMapInt32t64(m map[uint32]uint64, xl string, yl string, title string, vis string) {
	w := r.Writer
	// Write title:
	s := fmt.Sprintf("# Title: %s\n", title)
	w.WriteString(s)
	// Write header:
	s = fmt.Sprintf("#\t%s\t%s\t%s\n", xl, yl, vis)
	w.WriteString(s)
	w.WriteString("\n")

	for x, y := range m {
		s = fmt.Sprintf("\t%d\t%d\n", x, y)
		w.WriteString(s)
	}
	w.WriteString("\n")
	w.Flush()
}

func (r Report) ReportSliceInt32t64(x []uint32, y []uint64, xl string, yl string, title string, vis string) {
	if len(x) != len(y) {
		L.Fatalf("Cannot report %s: x/y length mismatch!\n", title)
	}
	w := r.Writer
	// Write title:
	s := fmt.Sprintf("# Title: %s\n", title)
	w.WriteString(s)
	// Write header:
	s = fmt.Sprintf("#\t%s\t%s\t%s\n", xl, yl, vis)
	w.WriteString(s)
	w.WriteString("\n")

	for i := 0; i < len(x); i++ {
		s = fmt.Sprintf("\t%d\t%d\n", x[i], y[i])
		w.WriteString(s)
	}
	w.WriteString("\n")
	w.Flush()
}

func (r Report) ReportSliceInt32f64(x []uint32, y []float64, xl string, yl string, title string, vis string) {
	if len(x) != len(y) {
		L.Fatalf("Cannot report %s: x/y length mismatch!\n", title)
	}
	w := r.Writer
	// Write title:
	s := fmt.Sprintf("# Title: %s\n", title)
	w.WriteString(s)
	// Write header:
	s = fmt.Sprintf("#\t%s\t%s\t%s\n", xl, yl, vis)
	w.WriteString(s)
	w.WriteString("\n")

	for i := 0; i < len(x); i++ {
		s = fmt.Sprintf("\t%d\t%f\n", x[i], y[i])
		w.WriteString(s)
	}
	w.WriteString("\n")
	w.Flush()
}

func (r Report) ReportSliceFloat64f64(x []float64, y []float64, xl string, yl string, title string, vis string) {
	if len(x) != len(y) {
		L.Fatalf("Cannot report %s: x/y length mismatch!\n", title)
	}
	w := r.Writer
	// Write title:
	s := fmt.Sprintf("# Title: %s\n", title)
	w.WriteString(s)
	// Write header:
	s = fmt.Sprintf("#\t%s\t%s\t%s\n", xl, yl, vis)
	w.WriteString(s)
	w.WriteString("\n")

	for i := 0; i < len(x); i++ {
		s = fmt.Sprintf("\t%f\t%f\n", x[i], y[i])
		w.WriteString(s)
	}
	w.WriteString("\n")
	w.Flush()
}
