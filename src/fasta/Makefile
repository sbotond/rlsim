include $(GOROOT)/src/Make.inc

TARG=fasta
GOFILES=\
    common.go \
    interfaces.go \
    seq.go \
    fasta.go \
    read.go \
    write.go \
    validate.go

include $(GOROOT)/src/Make.pkg

t:
	gotest
fmt:
	for f in *.go; do gofmt $$f > tmp.fmt; mv tmp.fmt $$f; done
com: fmt
	git commit -a

