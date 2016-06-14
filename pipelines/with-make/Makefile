SPECIES="Salmonella enterica"
READS=2492428
REFERENCE_URL=http://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000988525.2_ASM98852v2/GCA_000988525.2_ASM98852v2_genomic.fna.gz

# Toggle these variable declarations to switch between

# 1. no trimming/filtering
# PRECONDITION_TO_USE=fastq-dump.log
# READS_1=reads_1.fastq.gz
# READS_2=reads_2.fastq.gz

# 2. kmc
# PRECONDITION_TO_USE=reads.trim.pe.kmc.fastq.gz
# READS_1=reads.trim.pe.kmc.fastq.gz

# 3. khmer
PRECONDITION_TO_USE=reads.trim.pe.khmer.fastq.gz
READS_1=reads.trim.pe.khmer.fastq.gz

KMERSIZE=20
PLOTXMAX=60
PLOTYMAX=1200000
# hardcoded, or don't use the make all rule.
MINCOVERAGE=5

THREADS=2
MEMORYGB=4
# need to mkdir this
TMPDIR=${HOME}/kmc-temp
TMPFILE=./tmp

# File paths
TRIMMOMATIC=../bin/trimmomatic-0.36.jar
ADAPTERS=../adapters
PLOT_KMC=../bin/plotKMC.R
PLOT_KHMER=../bin/plotKhmer.R

# === MAIN PIPELINE ===

all: reads.vcf

# === GET REFERENCE ===

reference.genomic.fna.gz:	
	@echo "=== Running $@ ==="
	@export start=`date +%s`; \
	curl -s ${REFERENCE_URL} -o $@; \
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

# === TRIMMING ===

trim.happened: fastq-dump.log
	@echo "=== Running $@ ==="
	@export start=`date +%s`; \
	java -jar ${TRIMMOMATIC} PE -phred33 \
	  reads_1.fastq.gz reads_2.fastq.gz \
	  reads_1.trim.pe.fastq.gz reads_1.trim.se.fastq.gz \
	  reads_2.trim.pe.fastq.gz reads_2.trim.se.fastq.gz \
	  ILLUMINACLIP:${ADAPTERS}/TruSeq3-PE.fa:2:30:10 \
	  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 > $@; \
	export end=`date +%s`; \
	export runtime=$$(($$end-$$start)); \
	cat $@; \
	echo "Target $@ took $$runtime seconds"

reads.trim.pe.fastq.gz: trim.happened
	@echo "=== Running $@ ==="
	@export start=`date +%s`; \
	seqtk mergepe reads_1.trim.pe.fastq.gz reads_2.trim.pe.fastq.gz | gzip > $@; \
	export end=`date +%s`; \
	export runtime=$$(($$end-$$start)); \
	echo "Target $@ took $$runtime seconds"

# === FILTERING ===

# With refresh-bio/KMC
reads.trim.pe.kmc: reads.trim.pe.fastq.gz
	@echo "=== Running $@ ==="
	@export start=`date +%s`; \
	kmc -k${KMERSIZE} -m${MEMORYGB} -t${THREADS} $< $@ ${TMPDIR}; \
	export end=`date +%s`; \
	export runtime=$$(($$end-$$start)); \
	echo "Target $@ took $$runtime seconds"

reads.trim.pe.kmc.hist.txt: reads.trim.pe.kmc
	@echo "=== Running $@ ==="
	@export start=`date +%s`; \
	kmc_tools histogram $< $@; \
	export end=`date +%s`; \
	export runtime=$$(($$end-$$start)); \
	echo "Target $@ took $$runtime seconds"

reads.trim.pe.kmc.hist.txt.png: reads.trim.pe.kmc.hist.txt
	@echo "=== Running $@ ==="
	@export start=`date +%s`; \
	${PLOT_KMC} $< $@ ${PLOTXMAX} ${PLOTYMAX}; \
	export end=`date +%s`; \
	export runtime=$$(($$end-$$start)); \
	echo "Target $@ took $$runtime seconds"

reads.trim.pe.kmc.fastq.gz: reads.trim.pe.kmc.hist.txt.png reads.trim.pe.fastq.gz
	@echo "=== Running $@ ==="
	@export start=`date +%s`; \
	kmc_tools filter reads.trim.pe.kmc -cx${MINCOVERAGE} reads.trim.pe.fastq.gz -ci0 -cx0 $@; \
	export end=`date +%s`; \
	export runtime=$$(($$end-$$start)); \
	echo "Target $@ took $$runtime seconds"

# With khmer
reads.trim.pe.fastq.gz.kh: reads.trim.pe.fastq.gz
	@echo "=== Running $@ ==="
	@export start=`date +%s`; \
	load-into-counting.py -N 4 -k ${KMERSIZE} -M ${MEMORYGB}e9 -T ${THREADS} $@ $<; \
	export end=`date +%s`; \
	export runtime=$$(($$end-$$start)); \
	echo "Target $@ took $$runtime seconds"

reads.trim.pe.fastq.gz.kh.hist: reads.trim.pe.fastq.gz.kh
	@echo "=== Running $@ ==="
	@export start=`date +%s`; \
	abundance-dist.py $< reads.trim.pe.fastq.gz $@; \
	export end=`date +%s`; \
	export runtime=$$(($$end-$$start)); \
	echo "Target $@ took $$runtime seconds"

reads.trim.pe.fastq.gz.kh.hist.png: reads.trim.pe.fastq.gz.kh.hist
	@echo "=== Running $@ ==="
	@export start=`date +%s`; \
	${PLOT_KHMER} $< $@ ${PLOTXMAX} ${PLOTYMAX}; \
	export end=`date +%s`; \
	export runtime=$$(($$end-$$start)); \
	echo "Target $@ took $$runtime seconds"

reads.trim.pe.khmer.fastq.gz: reads.trim.pe.fastq.gz.kh reads.trim.pe.fastq.gz
	@echo "=== Running $@ ==="
	@export start=`date +%s`; \
	filter-abund.py -T ${THREADS} -C ${MINCOVERAGE} $< -o $@ reads.trim.pe.fastq.gz; \
	export end=`date +%s`; \
	export runtime=$$(($$end-$$start)); \
	echo "Target $@ took $$runtime seconds"

# === MAP READS TO REFERENCE ===

bwa-index.log: reference.genomic.fna.gz
	@echo "=== Running $@ ==="
	@export start=`date +%s`; \
	bwa index $< > $@; \
	export end=`date +%s`; \
	export runtime=$$(($$end-$$start)); \
	echo "Target $@ took $$runtime seconds"

reads.sam: bwa-index.log ${PRECONDITION_TO_USE}
	@echo "=== Running $@ ==="
	@export start=`date +%s`; \
	bwa mem -t ${THREADS} reference.genomic.fna.gz ${READS_1} ${READS_2} > $@; \
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
	samtools sort -@ ${THREADS} $< -o ANYTHING > $@; \
	export end=`date +%s`; \
	export runtime=$$(($$end-$$start)); \
	echo "Target $@ took $$runtime seconds"

reads.bam.bai: reads.bam
	@echo "=== Running $@ ==="
	@export start=`date +%s`; \
	samtools index $<; \
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
	-rm *.fastq.gz.*
	-rm *.log
	-rm trim.happened
	-rm *.kmc
	-rm *.kmc.*
	-rm *.sam
	-rm *.bam.bai
	-rm *.bam
	-rm reads.vcf
	-rm ${TMPFILE}
