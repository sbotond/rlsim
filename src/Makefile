all: rlsim tools

.PHONY: ct  com push fetch gt fmt fix

TOOLS=../tools

# Utility targets:
ct:
	git log --graph
fmt: *.go
	go fmt *.go
com: fmt
	git commit -a
push:
	git push --all
fetch:
	git fetch --all
gt:
	gitk --all
fix:
	go fix .
# Run standard test:
t: rlsim
	./rlsim -v -j test/basic/raw_params.json -n 2000 ./test/basic/test_transcripts.fas > /dev/null
	$(TOOLS)/plot_rlsim_report
# Run small test:
st: rlsim
	./rlsim -v -n 10 ./test/basic/test_transcripts.fas 
	$(TOOLS)/plot_rlsim_report
# Run big test with profiling:
rprof: rlsim
	./rlsim -si 7567 -v -n 3000000 -t 4 -f "after_prim_double" -prof prof.out ./test/basic/test_transcripts_full.fas > /dev/null
	$(TOOLS)/plot_rlsim_report
# View profiling results:
vprof:
	go tool pprof ./rlsim prof.out

# Targets for simple comparison of CPU usage:
NR_TIME_FRAGS	=1000000
TIME_EXPR_MUL	=1.0
TIME_CORES		=2
gtime: rlsim
	(/usr/bin/time -f "%U" ./rlsim -si 7567 -n $(NR_TIME_FRAGS) -m $(TIME_EXPR_MUL) -t $(TIME_CORES) -f "after_prim_double" ./test/basic/test_transcripts_full.fas > /dev/null) 2>>gtimes.tab 
gcctime: rlsim.gcc
	(/usr/bin/time -f "%U" ./rlsim.gcc -si 7567 -n $(NR_TIME_FRAGS) -m $(TIME_EXPR_MUL) -t $(TIME_CORES) -f "after_prim_double" ./test/basic/test_transcripts_full.fas > /dev/null) 2>>gcctimes.tab
times: gtime gcctime

# Main target:
TARG=rlsim
GOFILES=\
	main.go\
	args.go\
	frag.go\
	fragmentor.go\
	fragfilter.go\
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
	target_mix.go\
    version.go\
	fasta.go\
	parse_raw.go\
	raw_target.go\

rlsim: $(GOFILES)
	go build -o $(TARG) $(GOFILES)

# Build rlsim.gcc using gccgo:
GCC_FLAGS	=  -static-libgo -static-libgcc -O3 -march=native -funroll-loops -fomit-frame-pointer
gccbuild:
	go build -compiler gccgo -gccgoflags '$(GCC_FLAGS)' -o $(TARG).gcc $(GOFILES)

tools:
	@cd ../tools;make; cd ..

clean:
	rm -f $(TARG) $(TARG).gcc
	@cd ../tools;make clean; cd ..
