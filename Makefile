
include $(GOROOT)/src/Make.inc

GC:=$(GC) -I ./src/fasta/_obj
LD:=$(LD) -L ./src/fasta/_obj
 
TARG=rlsim
GOFILES=\
	src/main.go\
	src/args.go\
	src/frag.go\
	src/fragmentor.go\
	src/nnthermo.go\
	src/pool.go\
	src/report.go\
	src/sampler.go\
	src/target.go\
	src/transcript.go\
	src/utils.go\
	src/input.go\
	src/thermocycler.go\
	src/random.go\
	src/logging.go\
	src/fragstats.go\
    src/version.go\

include $(GOROOT)/src/Make.cmd 

