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

func RevCompDNA(s string) string {
	size := len(s)
	tmp := make([]byte, size)
	var inBase byte
	var outBase byte
	for i := 0; i < size; i++ {
		inBase = s[i]
		switch inBase {
		case 'A':
			outBase = 'T'
		case 'T':
			outBase = 'A'
		case 'G':
			outBase = 'C'
		case 'C':
			outBase = 'G'
		default:
			L.Fatalf("Invalid base: %s", inBase)
		}
		tmp[size-1-i] = outBase
	}
	return string(tmp)
}
