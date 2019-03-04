rule all:
	input:
		"freebayes_calls.report.html"
rule pretrim_fastqc:
	input:
		"/ufrc/kirst/d.ence/pinus_taeda_L/wegrzyn_elite_transcriptome_raw_data/test.{sample}.R1.fastq",
		"/ufrc/kirst/d.ence/pinus_taeda_L/wegrzyn_elite_transcriptome_raw_data/test.{sample}.R2.fastq"
		#"/ufrc/kirst/d.ence/pinus_taeda_L/wegrzyn_elite_transcriptome_raw_data/{sample}-Run02_S5_L001_R1_001.fastq.gz",	
		#"/ufrc/kirst/d.ence/pinus_taeda_L/wegrzyn_elite_transcriptome_raw_data/{sample}-Run02_S5_L001_R2_001.fastq.gz"
	output:
		"fastqc_out/pre_process/{sample}/",
		"logs/fastqc_pre_out/{sample}.log"
	log:
		"logs/fastqc_pre_out/{sample}.log"
	benchmark:
		"benchmarks/{sample}.fastqc.benchmark.txt"
	shell:
		"module load fastqc ; fastqc -f fastq -o {output} {input} &> {log}"			
rule trim_reads:
	input:
		fq1="/ufrc/kirst/d.ence/pinus_taeda_L/wegrzyn_elite_transcriptome_raw_data/test.{sample}.R1.fastq",
		fq2="/ufrc/kirst/d.ence/pinus_taeda_L/wegrzyn_elite_transcriptome_raw_data/test.{sample}.R2.fastq",
		#fq1="/ufrc/kirst/d.ence/pinus_taeda_L/wegrzyn_elite_transcriptome_raw_data/{sample}-Run02_S5_L001_R1_001.fastq.gz",	
		#fq2="/ufrc/kirst/d.ence/pinus_taeda_L/wegrzyn_elite_transcriptome_raw_data/{sample}-Run02_S5_L001_R2_001.fastq.gz",
		prev_log="logs/fastqc_pre_out/{sample}.log"
	output:
		fq1="trimmed_reads/{sample}/{sample}.trimmed.R1.fastq.gz",
		fq2="trimmed_reads/{sample}/{sample}.trimmed.R2.fastq.gz"		
	log:
		"logs/cutadapt/{sample}.log"
	benchmark:
		"benchmarks/{sample}.fastqc.benchmark.txt"
	shell:
		"module load cutadapt; cutadapt -u 10 -U 10 -q30,30 --minimum-length=50 -o {output.fq1} -p {output.fq2} {input.fq1} {input.fq2} &> {log}"
rule posttrim_fastqc:
	input:
		"trimmed_reads/{sample}/{sample}.trimmed.R1.fastq.gz",
		"trimmed_reads/{sample}/{sample}.trimmed.R2.fastq.gz"
	output:
		"fastqc_out/post_trim/{sample}/",
		"logs/fastqc_post_out/{sample}.log"
	log:
		"logs/fastqc_post_out/{sample}.log"
	benchmark:
		"benchmarks/{sample}.fastqc.benchmark.txt"
	shell:
		"module load fastqc ; fastqc -f fastq -o {output} {input} &> {log}"
rule hisat2_align:
	input:
		"trimmed_reads/{sample}/{sample}.trimmed.R1.fastq.gz",
		"trimmed_reads/{sample}/{sample}.trimmed.R2.fastq.gz",
		"logs/fastqc_post_out/{sample}.log"
	output:
		"mapped_reads/{sample}.bam"
	log:
		"logs/hisat2/{sample}.log"
	params:
		rg="-p 20 --rg-id={sample} --rg-sample={sample} --rg-library={sample} --rg-platform=illumina",
		ref=config["reference"]["V1_01"]["full"]	
	benchmark:
		"benchmakrs/{sample}.hisat2.benchmark.txt"
	shell:
		"module load hisat2; module load samtools; hisat2 --time --summary-file --met-file {log} -x {params.ref} {params.rg} {input} | samtools view - > {output}"	
rule tophat_process_align:
	input:
		#"/ufrc/kirst/d.ence/pinus_taeda_L/wegrzyn_elite_transcriptome_raw_data/{sample}-Run02_S5_L001_R1_001.fastq.gz",	
		#"/ufrc/kirst/d.ence/pinus_taeda_L/wegrzyn_elite_transcriptome_raw_data/{sample}-Run02_S5_L001_R2_001.fastq.gz",
		"trimmed_reads/{sample}/{sample}.trimmed.R1.fastq.gz",
		"trimmed_reads/{sample}/{sample}.trimmed.R2.fastq.gz",
		"logs/fastqc_post_out/{sample}.log"
	output:
		"tmp_tophat_{sample}/accepted_hits.bam"
	params:
		rg='-p 20 --rg-id={sample} --rg-sample={sample} --rg-library={sample} --rg-platform=illumina',
		ref_dir=config["reference"]["V1_01"]["dir"],
		ref_full=config["reference"]["V1_01"]["full"],
		output_dir='tmp_tophat_{sample}'
	log:
		"logs/tophat/{sample}.log"
	benchmark:
		"benchmarks/{sample}.tophat.benchmark.txt"
	shell:
		"module load tophat;  export BOWTIE2_INDEXES=""; tophat {params.rg} -o {params.output_dir} {params.ref_full} {input}"

rule samtools_sort:
	input:
		"tmp_tophat_{sample}/accepted_hits.bam"
	output:
		"sorted_reads/{sample}.sorted.bam"
	log:
		"logs/samtools_sort/{sample}.log"
	benchmark:
		"benchmarks/{sample}.sort.benchmark.txt"
	threads: 4		
	shell:
		"module load samtools; samtools sort --threads {threads} -T sorted_reads/{wildcards.sample} -O bam {input} > {output}"

rule samtools_rmdup:
	input:
		"sorted_reads/{sample}.sorted.bam"
	output:
		"rmduped_reads/{sample}.sorted.rmdup.bam"
	log:
		"logs/samtools_rmdup/{sample}.log"
	benchmark:
		"benchmarks/{sample}.rmdup.benchmark.txt"
	shell:
		"module load samtools; samtools rmdup {input} {output}"
rule stringtie:
	input:
		"realigned_bams/{sample}.sorted.rmdup.realigned.bam"
	output:
		"stringtie_gtfs/{sample}.stringtie.out.gtf"
	log:
		"logs/stringtie/{sample}.log"
	benchmark:
		"benchmarks/{sample}.stringtie.benchmark.txt"
	shell:
		"module load stringtie ; stringtie -o {output} {input}"


rule samtools_index_rmduped:
	input:
		"rmduped_reads/{sample}.sorted.rmdup.bam"
	output:
		"rmduped_reads/{sample}.sorted.rmdup.bam.bai"
	log:
		"logs/samtools_index_sorted{sample}.log"
	benchmark:
		"benchmarks/{sample}.index_sorted.benchmark.txt"
	shell:
		"module load samtools; samtools index {input}"
		

rule samtools_index_sorted:
	input:
		"sorted_reads/{sample}.sorted.rmdup.bam"
	output:
		"sorted_reads/{sample}.sorted.rmdup.bam.bai"
	log:
		"logs/samtools_index_sorted/{sample}.log"
	benchmark:
		"benchmarks/{sample}.index_sorted.benchmark.txt"	
	shell:
		"module load samtools; samtool index {input}"
rule gatk_indel_creator:
	input:
		"rmduped_reads/{sample}.sorted.rmdup.bam"
	output:
		"realigner_intervals/{sample}.intervals"
	params:
		ref=config["reference"]["V1_01"]["full"]
	log:
		"logs/realigner_intervals/{sample}.intervals.log"
	benchmark:
		"benchmarks/{sample}.intervals.benchmark.txt"
	shell:
		"module load gatk;  java -jar -Xmx4g /apps/gatk/3.7.0/GenomeAnalysisTK.jar -T RealignerTargetCreator -I {input} -o {output} -R {params.ref}"
rule gatk_indel_realign:
	input:
		bam="rmduped_reads/{sample}.sorted.rmdup.bam",
		interval="realigner_intervals/{sample}.intervals"
	output:
		"realigned_bams/{sample}.sorted.rmdup.realigned.bam"
	params:
		#ref=config["reference"]["V1_01"]["full"]
		ref="/home/d.ence/projects/ref_genomes/loblolly_pine/V1_1/greater_than_3k_ref/ptaeda.v1.01.fa.masked.longer_than3k.fasta"
	log:
		"logs/realigner/{sample}.realigner.log"
	benchmark:
		"benchmarks/{sample}.intervals.benchmark.txt"
	shell:
		"module load gatk; java -jar -Xmx4g /apps/gatk/3.7.0/GenomeAnalysisTK.jar -T IndelRealigner -R {params.ref} -I {input.bam} -targetIntervals {input.interval} -o {output}"
rule samtools_index_realigned:
	input:
		"realigned_bams/{sample}.sorted.rmdup.realigned.bam"
	output:
		"realigned_bams/{sample}.sorted.rmdup.realigned.bam.bai"
	log:
		"logs/samtools_index_realigned/{sample}.log"
	benchmark:
		"benchmarks/{sample}.index_realigned.benchmark.txt"
	shell:
		"module load samtools; samtool index {input}"

rule freebayes:
        input:
                ref=config["reference"]["V1_01"]["full"],
                bam=expand("rmduped_reads/{sample}.sorted.rmdup.realigned.bam",sample=config["samples"]),
		bai=expand("rmduped_reads/{sample}.sorted.rmdup.bam.bai",sample=config["samples"])                
                #bam=expand("indel_realigned/{sample}.sorted.rmdup.bam",sample=config["samples"])
                #pop_file=config["pop_file"]
        output:
                "calls/all_samples.ptaeda.wegrzyn_transcriptomes.freebayes_populations.vcf"
        threads: 4
        log:
                "logs/freebayes/all_samples.log"
        benchmark:
                "benchmarks/all_samples.freebayes.benchmark.txt"
        shell:
                #"module load freebayes; freebayes --populations {input.pop_file} -f {input.ref} {input.bam} > {output}"
                "module load freebayes; freebayes -f {input.ref} {input.bam} > {output}"

rule report:
	input:
		#T1="calls/all_samples.ptaeda.wegrzyn_transcriptomes.freebayes_populations.vcf",
		expand("stringtie_gtfs/{sample}.stringtie.out.gtf",sample=config["samples"])
	output:
		"freebayes_calls.report.html"
	run:
		from snakemake.utils import report
		with open(input) as gtf:
			n_calls = sum(1 for l in gtf if not l.startswith("#"))

		report("""
		A first test of a reference-based transcript assembly workflow
		=====================================================
		
		Reads were mapped to the Ptaeda v1.1 reference assembly 
		reference genome and variants were called jointly with
		freebayes.

		This resulted in {n_calls} variants (see Table T2_).
		""", output[0], T2=input[0])
