.PHONY: t bt

all: 
	@echo "making effest"
	@cd ./effest_src; make; cp effest ../

t_sel:
	./sel ../src/test/basic/test_transcripts.fas

t_pb_plot:
	./pb_plot -f test/ref.fas -i test/iso.list test/aln2.sam

t_cov_cmp:
	./cov_cmp -c 700 -t 10 -l 6000 -p cov_cmp.pk -f test/ref.fas -i test/iso.list test/aln1.sam test/aln2.sam

t: t_sel t_pb_plot t_cov_cmp
	@echo
clean:
	rm -fr effest
	@cd effest_src; make clean; cd ..
