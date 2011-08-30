package fasta

import (
	"log"
	"os"
)

var L *log.Logger

func init() {
	// Initilaize logging:
	L = log.New(os.Stderr, "[fasta] ", 0)
}
