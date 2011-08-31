
.PHONY: t st rprof vprof

TOOLS=../../../tools
RLSIM=../../rlsim

# Run standard test:
t: 
	$(RLSIM) -v -d "1:n:(300,20,100,5000)" -eg "(8.0,0.5,0.95)" -el "(0.4,0.5,1.0)" -n 200000 test_transcripts.fas > /dev/null
	$(TOOLS)/plot_rlsim_report

# Run small test:
st: 
	$(RLSIM) -v -n 10 test_transcripts.fas 
	$(TOOLS)/plot_rlsim_report
# Run big test with profiling:
rprof: 
	$(RLSIM) -v -n 3000000 -t 4 -f "after_prim_double" -prof prof.out test_transcripts_full.fas > /dev/null
	$(TOOLS)/plot_rlsim_report
# View profiling results:
vprof:
	go tool pprof $(RLSIM) prof.out

clean:
	rm -fr rlsim_report.json prof.out

