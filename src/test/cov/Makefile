
# Generate and plot data under different priming methods (and PCR amplification).
# Additional dependencies: simNGS, bwa, samtools.
# 

.PHONY: all 

NR_FRAGS		= 40000
READ_LENGTH		= 50
PLOT_COV        = ../../../tools/plot_cov
RLSIM			= ../../rlsim
SIMNGS			= simNGS
RUNFILE			= s_4_0066.runfile
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
	rm aln.bam;\
	mv alns.bam $(2);samtools index $(2);\
	$(PLOT_COV) -r $(REF) -b $(2) -o cov_$(3).pdf;\
	rm $(2) $(2).bai;\

COM_PARAMS			= -p 1.0 -b 0.0 -g -t 4 -v -n $(NR_FRAGS) -d "1.0:n:(200,20,100,20000)" -a "1.0:g:(1,1,0,1)"

AFTER_PRIM_P		= -c 0  $(COM_PARAMS) -f "after_prim"
AFTER_NOPRIM_P		= -c 0 $(COM_PARAMS) -f "after_noprim"
AFTER_PRIM_DOUBLE_P	= -c 0 $(COM_PARAMS) -f "after_prim_double"
AFTER_NOPRIM_DOUBLE_P	= -c 0 $(COM_PARAMS) -f "after_noprim_double"
PRE_PRIM_P			= -c 0 $(COM_PARAMS) -f "pre_prim:1500"
PRIM_JUMP_P			= -c 0 $(COM_PARAMS) -f "prim_jump"

PCR_P				= -c 30 -eg "(8.0,0.2,0.95)" $(COM_PARAMS) -f "after_noprim"
FULL_P				= -c 30 -eg "(8.0,0.2,0.95)" $(COM_PARAMS) -f "after_prim"

all: t_after_prim t_after_noprim t_after_prim_double t_after_noprim_double t_pre_prim t_prim_jump t_pcr t_full

ref.fas: gen_ref
	@./gen_ref > ref.fas; bwa index ref.fas

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

clean:
	rm -f ref.fas*	
