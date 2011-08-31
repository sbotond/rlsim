# 
# Generate test data used by the additional tools.
# Additional dependencies: simNGS, bwa, samtools.
# 

.PHONY: test_data

RLSIM			= ../../src/rlsim
SIMNGS			= simNGS
RLSIM_PARAMS1	= -t 4 -n 50000 -f "after_prim_double" -d "1.0:n:(450,100,310,6000)" -eg "(8.0,0.5,1.0)"
RLSIM_PARAMS2	= -t 4 -n 50000 -f "after_prim_double" -d "1.0:n:(450,100,310,6000)" -e 1.0
RUNFILE			= ../../src/test/cov/s_4_0066.runfile
REF				= ref.fas

gen_data = \
	@echo Generating test data $(1); \
	$(RLSIM) $(1) $(REF) > frags.fas; \
	$(SIMNGS) -p paired -o fastq -O reads $(RUNFILE) < frags.fas;rm frags.fas 2>/dev/null; \
	bwa index ref.fas ;\
	bwa aln ref.fas reads_end1.fq > end1.sai ;\
	bwa aln ref.fas reads_end2.fq > end2.sai ;\
	bwa sampe ref.fas end1.sai end2.sai reads_end1.fq reads_end2.fq > aln.sam ;\
	rm reads_end*.fq ;\
	samtools view -b -S aln.sam > aln.bam ;\
	samtools sort -n aln.bam alns ;\
	rm aln.sam aln.bam *.sai ref.fas.*;\
	samtools view -h alns.bam > alns.sam ;\
	mv alns.sam $(2) ; rm alns.bam;\

test_data:
	$(call gen_data,$(RLSIM_PARAMS1),aln1.sam)
	$(call gen_data,$(RLSIM_PARAMS2),aln2.sam)

