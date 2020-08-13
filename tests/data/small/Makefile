# Helper Makefile to build the small example.

# FASTA file to extract region from.
REF ?= /fast/projects/cubit/current/static_data/reference/GRCh37/hs37d5/hs37d5.fa
# The region to extract.
REGION ?= "1:58,929,877-59,028,960"
# Variants per sample to pull.
NUM_VARS ?= 20

SHELL := bash
RESULTS := vars_all.vcf.gz \
	vars_all.vcf.gz.tbi \
	vars_ontarget.vcf.gz \
	vars_ontarget.vcf.gz.tbi \
	ref.fasta.fai \
	oma1.bed \
	vars_singleton1.vcf.gz \
	vars_singleton1.vcf.gz.tbi
TMP := vars_all.vcf \
	vars_singleton1.vcf \
	reads_singleton1_1.fq \
	reads_singleton1_2.fq

.PHONY: all
all: $(RESULTS)

.PHONY: clean
clean:
	rm -f $(RESULTS) $(TMP)

rand_%.bin:
	head -c 1000 /dev/urandom >$@

oma1.bed: oma1.bed.in
	awk \
		-v offset=58929876 \
		-F $$'\t' \
		'BEGIN { OFS = FS } ($$6 == "OMA1") { print "contig", $$2 - offset, $$3 - offset }' \
		$< \
	| sort -k 2,2n \
	> $@

ref.fasta:
	samtools faidx "$(REF)" "$(REGION)" \
	| sed -e 's/^>.*/>contig/' \
	> $@

ref.fasta.fai: ref.fasta
	samtools faidx $<

# Simulate all variants along the genome.
vars_all.vcf.gz vars_all.vcf.gz.tbi: ref.fasta.fai
	mason_variator -ir ref.fasta -ov vars_all.vcf --snp-rate 0.05
	sed -i -e 's/sim_snp/snp/' -e 's/1$$/0\/1/' vars_all.vcf
	bgzip -c vars_all.vcf > vars_all.vcf.gz
	rm -f vars_all.vcf
	tabix -f vars_all.vcf.gz

# Extract on-target variants.
vars_ontarget.vcf.gz vars_ontarget.vcf.gz.tbi: vars_all.vcf.gz oma1.bed
	bedtools intersect -header -a vars_all.vcf.gz -b oma1.bed \
	| bgzip -c \
	> vars_ontarget.vcf.gz
	tabix -f vars_ontarget.vcf.gz

# Sample on-target variants for one singleton.
vars_singleton1.vcf.gz vars_singleton1.vcf.gz.tbi: vars_ontarget.vcf.gz rand_1.bin
	zgrep '^#' $< \
	| sed -e 's/simulated/singleton1/' \
	> vars_singleton1.vcf

	zgrep -v '^#' $< \
	| shuf --random-source rand_1.bin -n $(NUM_VARS) \
	| sort -k 1,1g -k2,2n \
	| awk -v seed=1 -F $$'\t' 'BEGIN { OFS = FS; srand(seed); rate = 0.5; } { \
			x = rand(); \
			if (x < rate * rate) { \
				$$NF = "1|1"; \
			} else if (x < 0.5 + rate * rate / 2) { \
				$$NF = "0|1"; \
			} else { \
				$$NF = "1|0"; \
			} \
			print; \
		}' \
	>> vars_singleton1.vcf
	bgzip -c vars_singleton1.vcf > vars_singleton1.vcf.gz
	rm -f vars_singleton1.vcf
	tabix -f vars_singleton1.vcf.gz

# Generate reads and BAM file for singleton.
reads_singleton1_1.fq.gz reads_singleton1_2.fq.gz reads_singleton1.bam reads_singleton1.bam.bai: vars_singleton1.vcf.gz ref.fasta.fai
	len=$$(grep ^contig ref.fasta.fai | cut -f 2); \
	reads=$$(($$len * 40 / 100)); \
	mason_simulator \
		-ir ref.fasta \
		-n $$reads \
		-iv vars_singleton1.vcf.gz \
		-o reads_singleton1_1.fq \
		-or reads_singleton1_2.fq \
		-oa reads_singleton1.unsorted.bam

	gzip -c reads_singleton1_1.fq >reads_singleton1_1.fq.gz
	gzip -c reads_singleton1_2.fq >reads_singleton1_2.fq.gz
	rm -f reads_singleton1_1.fq reads_singleton1_2.fq

	samtools sort -O BAM -o reads_singleton1.bam reads_singleton1.unsorted.bam
	samtools index reads_singleton1.bam
	rm -f reads_singleton1.unsorted.bam

