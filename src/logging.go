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

import "io"
import "log"
import "os"

type Logger interface {
	SetLogLevel(level int)
	SetWriter(w io.Writer)
	SetPrefix(s string)
	Fatal(v ...interface{})
	Fatalf(format string, v ...interface{})
	Panic(v ...interface{})
	Panicf(format string, v ...interface{})
	Printf(format string, v ...interface{})
	Println(v ...interface{})
	PrintfV(format string, v ...interface{})
	PrintlnV(v ...interface{})
}

type Log struct {
	log    *log.Logger
	w      io.Writer
	level  int
	flag   int
	prefix string
}

func NewLog(level int) (l *Log) {
	l = new(Log)
	l.w = os.Stderr
	l.level = level
	l.flag = 4
	l.prefix = ""
	l.log = log.New(l.w, l.prefix, l.flag)
	return
}

func (l Log) Fatal(v ...interface{}) {
	l.log.Fatal(v...)
}

func (l Log) Fatalf(format string, v ...interface{}) {
	l.log.Fatalf(format, v...)
}

func (l Log) Panic(v ...interface{}) {
	l.log.Panic(v...)
}

func (l Log) Panicf(format string, v ...interface{}) {
	l.log.Panicf(format, v...)
}

func (l Log) Printf(format string, v ...interface{}) {
	l.log.Printf(format, v...)
}

func (l Log) Println(v ...interface{}) {
	l.log.Println(v...)
}

func (l Log) PrintfV(format string, v ...interface{}) {
	if l.level != 1 {
		return
	}
	l.log.Printf(format, v...)
}

func (l Log) PrintlnV(v ...interface{}) {
	if l.level != 1 {
		return
	}
	l.log.Println(v...)
}

func (l Log) SetLogLevel(level int) {
	switch level {
	case -1, 0, 1:
		l.level = level
	default:
		l.Fatalf("Invalid log level: %d\n", level)
	}
}

func (l Log) SetWriter(w io.Writer) {
	l.w = w
	l.log = log.New(l.w, l.prefix, l.flag)
}

func (l Log) SetPrefix(prefix string) {
	l.prefix = prefix
	l.log = log.New(l.w, l.prefix, l.flag)
}
