# Generate and plot data under different priming methods (and PCR amplification).
# Additional dependencies: simNGS, bwa, samtools.
# 

.PHONY: all 

NR_FRAGS		= 5000000
READ_LENGTH		= 75

RLSIM			= ../../rlsim
PLOT_COV		= ../../../tools/plot_cov
PLOT_COV		= echo 
PB_PLOT			= ../../../tools/pb_plot -w 15
SIMNGS			= simNGS
RUNFILE			= ../cov/s_4_0066.runfile
REF				= ref.fas

gen_plot = \
	@echo Generating test data $(1); \
	$(RLSIM) $(1) $(REF) > frags.fas; \
	rm rlsim_report.json;\
	$(SIMNGS) -n $(READ_LENGTH) -p paired -o fastq -O reads $(RUNFILE) < frags.fas;rm frags.fas 2>/dev/null; \
	bwa aln ref.fas reads_end1.fq > end1.sai ;\
	bwa aln ref.fas reads_end2.fq > end2.sai ;\
	bwa sampe ref.fas end1.sai end2.sai reads_end1.fq reads_end2.fq > aln.sam ;\
	rm reads_end*.fq end*.sai;\
	samtools view -b -S aln.sam > aln.bam ;rm aln.sam;\
	samtools sort aln.bam alns;\
	samtools sort -n aln.bam aln_ns ;rm aln.bam;\
	samtools view -h aln_ns.bam > aln_ns.sam;\
	$(PB_PLOT) -f ref.fas -r pb_$(3).pdf aln_ns.sam;rm aln_ns.sam aln_ns.bam;\
	mv alns.bam $(2);samtools index $(2);\
	$(PLOT_COV) -r $(REF) -b $(2) -o cov_$(3).pdf;\
	rm $(2) $(2).bai;\

NO_POLYA			= -a "1.0:g:(1,1,0,1)"
COM_PARAMS			= -p 40.0 -t 4 -v -n $(NR_FRAGS) -d  "1.0:n:(200, 20, 150, 100000)" $(NO_POLYA)

AFTER_PRIM_P		= -c 0 -e 1.0 $(COM_PARAMS) -f "after_prim" $(NO_POLYA)
AFTER_NOPRIM_P		= -c 0 -e 1.0 $(COM_PARAMS) -f "after_noprim" $(NO_POLYA)
AFTER_PRIM_DOUBLE_P	= -e 1.0 $(COM_PARAMS) -f "after_prim_double" $(NO_POLYA)
AFTER_NOPRIM_DOUBLE_P	= -e 1.0 $(COM_PARAMS) -f "after_noprim_double" $(NO_POLYA)
PRE_PRIM_P			= -e 1.0 $(COM_PARAMS) -f "pre_prim:500000" $(NO_POLYA)
PRIM_JUMP_P			= -e 1.0 $(COM_PARAMS) -f "prim_jump" $(NO_POLYA)

PCR_P				= -eg "(0.2, 0.5, 0.95)" $(COM_PARAMS) -f "after_noprim" $(NO_POLYA)
FULL_P				= -eg "(0.2, 0.5, 0.95)" $(COM_PARAMS) -f "after_prim" $(NO_POLYA)
FULL_A_P			= -eg "(0.2, 0.5, 0.95)" $(COM_PARAMS) -f "after_prim" -a "1.0:g:(200,1,1,250)"

all: t_after_prim t_after_noprim t_after_prim_double t_after_noprim_double t_pre_prim t_prim_jump t_pcr t_full t_full_polya

gen_ref: gen_ref.go
	go build -o gen_ref gen_ref.go

ref.fas: gen_ref 
	@echo Generating a *huge* random reference transcriptome.
	@./gen_ref -l 1000000 -e 100 -n 500 > ref.fas
	@echo Building index for reference transcriptome.
	@bwa index ref.fas

t_after_prim: ref.fas
	$(call gen_plot,$(AFTER_PRIM_P),alns_after_prim.bam,after_prim)

t_after_noprim: ref.fas
	$(call gen_plot,$(AFTER_NOPRIM_P),alns_after_noprim.bam,after_noprim)

t_after_prim_double: ref.fas
	$(call gen_plot,$(AFTER_PRIM_DOUBLE_P),alns_after_prim_double.bam,after_prim_double)

t_after_noprim_double: ref.fas
	$(call gen_plot,$(AFTER_NOPRIM_DOUBLE_P),alns_after_noprim_double.bam,after_noprim_double)

t_pre_prim: ref.fas
	$(call gen_plot,$(PRE_PRIM_P),alns_pre_prim.bam,pre_prim)

t_prim_jump: ref.fas
	$(call gen_plot,$(PRIM_JUMP_P),alns_prim_jump.bam,prim_jump)

t_pcr: ref.fas
	$(call gen_plot,$(PCR_P),alns_pcr.bam,pcr)

t_full: ref.fas
	$(call gen_plot,$(FULL_P),alns_full.bam,full)

t_full_polya: ref.fas
	$(call gen_plot,$(FULL_A_P),alns_full.bam,full_polya)

clean:
	rm ref.fas*
