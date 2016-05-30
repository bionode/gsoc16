SPECIES="Salmonella enterica"
READS=2492428

THREADS=2
MEMORYGB=4
TMPDIR=/Users/jmazz/Desktop/kmc-temp
TMPFILE=./tmp

# === MAIN PIPELINE ===

all: reads.vcf

# === GET REFERENCE ===

reference-urls.txt:
	@# This can take several minutes
	@# Can be made faster using an accession, but the point is the user need only
	@# know the species name
	@echo "=== Running $@ ==="
	@export start=`date +%s`; \
	bionode-ncbi urls assembly ${SPECIES} > $@; \
	export end=`date +%s`; \
	export runtime=$$(($$end-$$start)); \
	echo "Target $@ took $$runtime seconds"

reference.genomic.fna.gz: reference-urls.txt
	@# Decently sketchy. bionode-ncbi giving a lot of rna_from_genomic instead of genomic under fna key
	@# Take the first result that is not rna_from_genomic and has genomic.fna.gz
	@echo "=== Running $@ ==="
	@export start=`date +%s`; \
	export REFERENCE_URL=`cat $< | grep -v rna_from_genomic | head -n 1 | json | grep genomic.fna.gz | cut -d \" -f4`; \
	curl -s $$REFERENCE_URL -o $@; \
	export end=`date +%s`; \
	export runtime=$$(($$end-$$start)); \
	echo "Target $@ took $$runtime seconds"

# === GET READS ===

reads.sra:
	@echo "=== Running $@ ==="
	@export start=`date +%s`; \
	bionode-ncbi download sra ${READS} > ${TMPFILE}; \
	cp ${READS}/*.sra $@; \
	export end=`date +%s`; \
	export runtime=$$(($$end-$$start)); \
	echo "Target $@ took $$runtime seconds"

fastq-dump.log: reads.sra
	@echo "=== Running $@ ==="
	@export start=`date +%s`; \
	fastq-dump --split-files --skip-technical --gzip $< > $@; \
	export end=`date +%s`; \
	export runtime=$$(($$end-$$start)); \
	cat $@; \
	echo "Target $@ took $$runtime seconds"; \

# === FILTER ===
# TODO

# === MAP READS TO REFERENCE ===

bwa-index.log: reference.genomic.fna.gz
	@echo "=== Running $@ ==="
	@export start=`date +%s`; \
	bwa index $< > $@; \
	export end=`date +%s`; \
	export runtime=$$(($$end-$$start)); \
	echo "Target $@ took $$runtime seconds"

reads.sam: bwa-index.log fastq-dump.log
	@echo "=== Running $@ ==="
	@export start=`date +%s`; \
	bwa mem -t ${THREADS} reference.genomic.fna.gz *.fastq.gz > $@; \
	export end=`date +%s`; \
	export runtime=$$(($$end-$$start)); \
	echo "Target $@ took $$runtime seconds"

reads.unsorted.bam: reads.sam
	@echo "=== Running $@ ==="
	@export start=`date +%s`; \
	samtools view -@ ${THREADS} -bh $< > $@; \
	export end=`date +%s`; \
	export runtime=$$(($$end-$$start)); \
	echo "Target $@ took $$runtime seconds"

reads.bam: reads.unsorted.bam
	@echo "=== Running $@ ==="
	@export start=`date +%s`; \
	samtools sort -@ ${THREADS} $< -o $@; \
	export end=`date +%s`; \
	export runtime=$$(($$end-$$start)); \
	echo "Target $@ took $$runtime seconds"

reads.bam.bai: reads.bam
	@echo "=== Running $@ ==="
	@export start=`date +%s`; \
	samtools index $< $@; \
	export end=`date +%s`; \
	export runtime=$$(($$end-$$start)); \
	echo "Target $@ took $$runtime seconds"

reference.genomic.fna: reference.genomic.fna.gz
	@echo "=== Running $@ ==="
	@export start=`date +%s`; \
	bgzip -d $< --stdout > $@; \
	export end=`date +%s`; \
	export runtime=$$(($$end-$$start)); \
	echo "Target $@ took $$runtime seconds"

reads.vcf: reads.bam.bai reference.genomic.fna
	@echo "=== Running $@ ==="
	@export start=`date +%s`; \
	samtools mpileup -uf reference.genomic.fna reads.bam | bcftools call -c - > $@; \
	export end=`date +%s`; \
	export runtime=$$(($$end-$$start)); \
	echo "Target $@ took $$runtime seconds"

# === RULE TEMPLATE ===

time:
	@echo "=== Running $@ ==="
	@export start=`date +%s`; \
	sleep 2; \
	export end=`date +%s`; \
	export runtime=$$(($$end-$$start)); \
	echo "Target $@ took $$runtime seconds"

# === TEST and CLEAN ===

test:
	# Check that we have all the required programs
	@bionode-ncbi -v
	@fastq-dump --help
	# No --help or -v for this..
	#@bwa

clean:
	@# rm reference-urls.txt
	-rm -rf ${READS}
	-rm *.sra
	-rm *.fna
	-rm *.fna.*
	-rm *.fastq.gz
	-rm *.log
	-rm *.sam
	-rm *.bam.bai
	-rm *.bam
	-rm *.vcf
	-rm ${TMPFILE}
