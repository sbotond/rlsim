.PHONY: ct  com push fetch gt fmt

ct:
	git log --graph
com: fmt
	git commit -a
push:
	git push --all
fetch:
	git fetch --all
gt:
	gitk --all
t: rlsim
	./rlsim -n 2000 ../test/basic/test_transcripts.fas > /dev/null
st: rlsim
	./rlsim -n 200 ../test/basic/test_transcripts.fas 
bt: rlsim
	./rlsim -n 30000000 -t 4 -prof prof.out ../test/big_test/test_transcripts_full.fas > /dev/null
lt:
	./sel ../test/basic/raw_transcripts.fas
prof:
	gopprof ./rlsim prof.out


include $(GOROOT)/src/Make.inc

GC:=$(GC) -I ./fasta/_obj
LD:=$(LD) -L ./fasta/_obj
 
TARG=rlsim
GOFILES=\
	main.go\
	args.go\
	frag.go\
	fragmentor.go\
	nnthermo.go\
	pool.go\
	report.go\
	sampler.go\
	target.go\
	transcript.go\
	utils.go\
	input.go\
	thermocycler.go\
	random.go\
	logging.go\
	fragstats.go\
    version.go\

include $(GOROOT)/src/Make.cmd 

