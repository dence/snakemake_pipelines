rule all:
	input:
		"trinity_guided_assembly/Trinity-GG.fasta"

rule sra_prefetch:
	output:
		temp("sra/sra/{sample}.sra")	
	params:
		SRA="{sample}"
	log:
		"logs/sra_prefetch/{sample}.log"
	shell:
		#This command-line should work, but results in a job failure message for some reason. DENCE 03/04/2019
		#"module load sra/2.8.2.1; echo $PWD >> {log} ; vdb-config --set /repository/user/main/public/root=$PWD/sra &>> {log}; prefetch {params.SRA} &>> {log}"
		"module load sra/2.8.2.1; echo $PWD >> {log} ; prefetch {params.SRA} &>> {log}"

rule dump_accessions:
	input:
		#"/home/d.ence/projects/Pnigra_Sempervirens/supernova_genome_annotation/evidence/RNAseq/SRA_downloads/sra/{sample}.sra"
		"sra/sra/{sample}.sra"	
	output:
		temp("tmp_fastqs/{sample}.fastq.bz2"),
		"logs/fastq_dump/{sample}.log"
	log:
		"logs/fastq_dump/{sample}.log"
	shell:
		"module load sra; fastq-dump --bzip2 --outdir tmp_fastqs {input} &> {log}"
rule pretrim_fastqc:
	input:
		"tmp_fastqs/{sample}.fastq.bz2"
	output:
		"fastqc_out/pre_process/{sample}/",
		"logs/fastqc_pre_out/{sample}.log"
	log:
		"logs/fastqc_pre_out/{sample}.log"
	benchmark:
		"benchmarks/{sample}.fastqc_pre.benchmark.txt"
	shell:
		"module load fastqc ; fastqc -f fastq -o {output} {input} &> {log}"

rule trim_reads:
	input:
		fq1="tmp_fastqs/{sample}.fastq.bz2",
		prev_log="logs/fastqc_pre_out/{sample}.log"
	output:
		fq1=temp("trimmed_reads/{sample}/{sample}.trimmed.R1.fastq.gz")
	log:
		"logs/cutadapt/{sample}.log"
	benchmark:
		"benchmarks/{sample}.fastqc.benchmark.txt"
	shell:
		"module load gcc/5.2.0; module load cutadapt; cutadapt -u -10 -q30,30 --minimum-length=25 -o {output.fq1} {input.fq1} &> {log}"
rule posttrim_fastqc:
	input:
		"trimmed_reads/{sample}/{sample}.trimmed.R1.fastq.gz"
	output:
		"fastqc_out/post_trim/{sample}/",
		"logs/fastqc_post_out/{sample}.log"
	log:
		"logs/fastqc_post_out/{sample}.log"
	benchmark:
		"benchmarks/{sample}.fastqc_post.benchmark.txt"
	shell:
		"module load fastqc ; fastqc -f fastq -o {output} {input} &> {log}"


rule star_align:
	input:
		fq1="trimmed_reads/{sample}/{sample}.trimmed.R1.fastq.gz",
		old_log="logs/fastqc_post_out/{sample}.log",
		ref="/home/d.ence/projects/Pnigra_Sempervirens/EVRGREEN/supernova_ref_assembly/supernova_PA/"
	output:
		filename=temp("mapped_reads/{sample}.Aligned.out.bam")
	params:
		prefix="mapped_reads/{sample}."
	log:
		"logs/rnastar/{sample}.log"
	benchmark:
		"benchmarks/{sample}.rnastar.benchmark.txt"
	threads:
		10
	shell:
		"module load samtools; module load gcc/5.2.0; module load star/2.6.0a; STAR --runThreadN {threads} --genomeDir {input.ref} --readFilesIn {input.fq1} --readFilesCommand zcat -c --outFileNamePrefix {params.prefix} --outSAMtype BAM Unsorted > {log}"

rule samtools_sort:
	input:
		#"mapped_reads/{sample}.hisat2.bam"	
		"mapped_reads/{sample}.Aligned.out.bam"
	output:
		#"sorted_reads/{sample}.hisat2.sorted.bam"
		temp("sorted_reads/{sample}.star.sorted.bam")
	log:
		"logs/samtools_sort/{sample}.log"
	benchmark:
		"benchmarks/{sample}.sort.benchmark.txt"
	threads: 4		
	shell:
		"module load samtools; samtools sort --threads {threads} -T sorted_reads/{wildcards.sample} -O bam {input} > {output}"

rule samtools_rmdup:
	input:
		#"sorted_reads/{sample}.hisat2.sorted.bam"
		"sorted_reads/{sample}.star.sorted.bam"
	output:
		#"rmduped_reads/{sample}.hisat2.sorted.rmdup.bam"
		protected("rmduped_reads/{sample}.star.sorted.rmdup.bam")
	log:
		"logs/samtools_rmdup/{sample}.log"
	benchmark:
		"benchmarks/{sample}.rmdup.benchmark.txt"
	shell:
		"module load samtools; samtools rmdup {input} {output}"

rule samtools_index_rmduped:
	input:
		#"rmduped_reads/{sample}.hisat2.sorted.rmdup.bam"
		"rmduped_reads/{sample}.star.sorted.rmdup.bam"
	output:
		#"rmduped_reads/{sample}.hisat2.sorted.rmdup.bam.bai"
		"rmduped_reads/{sample}.star.sorted.rmdup.bam.bai"
	log:
		"logs/samtools_index_sorted.{sample}.log"
	benchmark:
		"benchmarks/{sample}.index_sorted.benchmark.txt"
	shell:
		"module load samtools; samtools index {input}"

rule samtools_index_sorted:
	input:
		"sorted_reads/{sample}.star.sorted.rmdup.bam"
		#"sorted_reads/{sample}.hisat2.sorted.rmdup.bam"
	output:
		"sorted_reads/{sample}.star.sorted.rmdup.bam.bai"
		#"sorted_reads/{sample}.hisat2.sorted.rmdup.bam.bai"
	log:
		"logs/samtools_index_sorted/{sample}.log"
	benchmark:
		"benchmarks/{sample}.index_sorted.benchmark.txt"	
	shell:
		"module load samtools; samtool index {input}"
rule merge_bams:
	input:
		#expand("rmduped_reads/{sample}.hisat2.sorted.rmdup.bam",sample=config["samples"])
		expand("rmduped_reads/{sample}.star.sorted.rmdup.bam",sample=config["samples"])
	output:
		#"merged_bams/" + config["accession"][0] + ".hisat2.sorted.rmdup.merged.bam"
		"merged_bams/" + config["accession"][0] + ".star.sorted.rmdup.merged.bam"
	log:
		"logs/picard_mergesamfiles.log"
	params:
		"SORT_ORDER=coordinate ASSUME_SORTED=true"
	#wrapper:
	#	"0.23.1/bio/picard/mergesamfiles"
	run:
		#hacking the wrapper to get the module loaded 
		from snakemake.shell import shell
		inputs = " ".join("INPUT={}".format(in_) for in_ in input)
		#log = snakemake.log_gmt_shell(stdout=False,stderr=True)
		shell(
			"module load picard;"
			" picard"
			" MergeSamFiles"
			" {params}"
			" {inputs}"
			" OUTPUT={output[0]}"
			" &> {log}")
		
#rule picard_merge_sams:
#	input:
#		bam=expand("rmduped_reads/{sample}.hisat2.sorted.rmdup.bam",sample=config["samples"]),
#		bai=expand("rmduped_reads/{sample}.hisat2.sorted.rmdup.bam",sample=config["samples"])
#	output:
#		"merged_bams/" + config["accession"][0] + ".hisat2.sorted.rmdup.merged.bam"
#	log:
#		"logs/picard_merge_bams/" + config["accession"][0] + ".merge_bam.log"
#	shell:
#		"module load samtools; module load picard; picard MergeSamFiles INPUT={input.bam} OUTPUT={output} SORT_ORDER=coordinate ASSUME_SORTED=true &> {log} ;samtools index {output}"
		
rule trinity_genome_guided:
	input:
		#merged_bam="merged_bams/" + config["accession"] + ".hisat2.sorted.rmdup.merged.bam"
		merged_bam="merged_bams/" + config["accession"][0] + ".star.sorted.rmdup.merged.bam"
#		#fq1=expand("trimmed_reads/{sample}/{sample}.trimmed.R1.fastq.gz",sample=config["samples"])	
#		fq1=expand("trimmed_reads/{sample}/{sample}.trimmed.R1.fastq.gz",sample=config["samples"]),	
#		fq2=expand("trimmed_reads/{sample}/{sample}.trimmed.R1.fastq.gz",sample=config["samples"])	
	params:
		"--seqType fq --SS_lib_type " + config["layout"][0] + " --max_memory 50G --genome_guided_max_intron 200000"
	output:
		output_dir="trinity_guided_assembly",
		output_fasta="trinity_guided_assembly/Trinity-GG.fasta",	
		logs="logs/trinity_genome_guided/" + config["accession"][0] +  ".log"
	log:
		"logs/trinity_genome_guided/" + config["accession"][0] + ".log"
	benchmark:
		"benchmarks/" + config["accession"][0] + ".trinity.benchmarks.txt"
	threads: 
		10
	shell:
		"module load gcc/5.2.0; module load trinity/r20180213-2.6.5; Trinity --CPU {threads} --genome_guided_bam {input.merged_bam} --output {output.output_dir} {params} &> {log}"

