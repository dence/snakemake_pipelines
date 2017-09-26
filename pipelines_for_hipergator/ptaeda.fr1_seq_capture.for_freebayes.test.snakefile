rule all:
	input:
		"freebayes_calls.report.html"
rule pretrim_fastqc:
	input:
		fq1="/home/d.ence/projects/pinus_taeda_L/Fr1_project/aligning_Fr1_samples/aligning_to_V1_1/test.{sample}.R1.fastq",
		fq2="/home/d.ence/projects/pinus_taeda_L/Fr1_project/aligning_Fr1_samples/aligning_to_V1_1/test.{sample}.R2.fastq"
                #"/ufrc/kirst/d.ence/pinus_taeda_L/wegrzyn_elite_transcriptome_raw_data/test.{sample}.R1.fastq",
                #"/ufrc/kirst/d.ence/pinus_taeda_L/wegrzyn_elite_transcriptome_raw_data/test.{sample}.R2.fastq"
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
		fq1="/home/d.ence/projects/pinus_taeda_L/Fr1_project/aligning_Fr1_samples/aligning_to_V1_1/test.{sample}.R1.fastq",
		fq2="/home/d.ence/projects/pinus_taeda_L/Fr1_project/aligning_Fr1_samples/aligning_to_V1_1/test.{sample}.R2.fastq",
                #fq1="/ufrc/kirst/d.ence/pinus_taeda_L/wegrzyn_elite_transcriptome_raw_data/test.{sample}.R1.fastq",
                #fq2="/ufrc/kirst/d.ence/pinus_taeda_L/wegrzyn_elite_transcriptome_raw_data/test.{sample}.R2.fastq",
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
rule bowtie2_map:
	input:
		ref=config["reference"]["V1_01"]["full"],
		fq1="trimmed_reads/{sample}/{sample}.trimmed.R1.fastq.gz",
		fq2="trimmed_reads/{sample}/{sample}.trimmed.R2.fastq.gz",
		old_log="logs/fastqc_post_out/{sample}.log"
		#fq1="/home/d.ence/projects/pinus_taeda_L/Fr1_project/aligning_Fr1_samples/aligning_to_V1_1/test.{sample}.R1.fastq",
		#fq2="/home/d.ence/projects/pinus_taeda_L/Fr1_project/aligning_Fr1_samples/aligning_to_V1_1/test.{sample}.R2.fastq"
		#fq1="/home/d.ence/projects/slash_pine/downsampled_slash_pines_for_development/prefix_{sample}.R1.fastq",
		#fq2="/home/d.ence/projects/slash_pine/downsampled_slash_pines_for_development/prefix_{sample}.R2.fastq"
	output:
		temp(tmp_sam="mapped_reads/{sample}.bowtie2.sam"),
		bam="mapped_reads/{sample}.bowtie2.bam"
	params:
		rg="--rg-id {sample} --rg 'SM:{sample}\\tPL:illumina'",
	log:
		"logs/bowtie2/{sample}.log"
	benchmark:
		"benchmarks/{sample}.bowtie2.benchmark.txt"
	threads: 4
	shell:
		"module load bowtie2; bowtie2 -q -t --mm -p {threads} {params.rg} -x {input.ref} -1 {input.fq1} -2 {input.fq2} -S {output.tmp_sam}  ; samtools view -Sb {output.tmp_sam} > {output.bam}"

rule samtools_sort:
	input:
		"mapped_reads/{sample}.bowtie2.bam"
	output:
		"sorted_reads/{sample}.bowtie2.sorted.bam"
	log:
		"logs/samtools_sort/{sample}.log"
	benchmark:
		"benchmarks/{sample}.sort.benchmark.txt"
	threads: 4		
	shell:
		"module load samtools; samtools sort --threads {threads} -T sorted_reads/{wildcards.sample} -O bam {input} > {output}"

rule samtools_rmdup:
	input:
		"sorted_reads/{sample}.bowtie2.sorted.bam"
	output:
		"rmduped_reads/{sample}.bowtie2.sorted.rmdup.bam"
	log:
		"logs/samtools_rmdup/{sample}.log"
	benchmark:
		"benchmarks/{sample}.rmdup.benchmark.txt"
	shell:
		"module load samtools; samtools rmdup {input} {output}"

rule samtools_index_rmduped:
	input:
		"rmduped_reads/{sample}.bowtie2.sorted.rmdup.bam"
	output:
		"rmduped_reads/{sample}.bowtie2.sorted.rmdup.bam.bai"
	log:
		"logs/smatools_index_sorted{sample}.log"
	benchmark:
		"benchmarks/{sample}.index_sorted.benchmark.txt"
	shell:
		"module load samtools; samtools index {input}"

rule samtools_index_sorted:
	input:
		"sorted_reads/{sample}.bowtie2.sorted.rmdup.bam"
	output:
		"sorted_reads/{sample}.bowtie2.sorted.rmdup.bam.bai"
	log:
		"logs/samtools_index_sorted/{sample}.log"
	benchmark:
		"benchmarks/{sample}.index_sorted.benchmark.txt"	
	shell:
		"module load samtools; samtool index {input}"
rule gatk_indel_creator:
        input:
                "rmduped_reads/{sample}.bowtie2.sorted.rmdup.bam"
        output:
                "realigner_intervals/{sample}.bowtie2.intervals"
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
                bam="rmduped_reads/{sample}.bowtie2.sorted.rmdup.bam",
                interval="realigner_intervals/{sample}.bowtie2.intervals"
        output:
                "realigned_bams/{sample}.bowtie2.sorted.rmdup.realigned.bam"
        params:
                ref=config["reference"]["V1_01"]["full"]
        log:
                "logs/realigner/{sample}.realigner.log"
        benchmark:
                "benchmarks/{sample}.intervals.benchmark.txt"
        shell:
                "module load gatk; java -jar -Xmx4g /apps/gatk/3.7.0/GenomeAnalysisTK.jar -T IndelRealigner -R {params.ref} -I {input.bam} -targetIntervals {input.interval} -o {output}"
rule samtools_index_realigned:
        input:
                "realigned_bams/{sample}.bowtie2.sorted.rmdup.realigned.bam"
        output:
                "realigned_bams/{sample}.bowtie2.sorted.rmdup.realigned.bam.bai"
        log:
                "logs/samtools_index_realigned/{sample}.log"
        benchmark:
                "benchmarks/{sample}.index_realigned.benchmark.txt"
        shell:
                "module load samtools; samtool index {input}"

rule freebayes:
        input:
                ref=config["reference"]["V1_01"]["full"],
                #bam=expand("realigned_bams/{sample}.bowtie2.sorted.rmdup.realigned.bam",sample=config["samples"])
                #bai=expand("realigned_bams/{sample}.bowtie2.sorted.rmdup.realigned.bai",sample=config["samples"])
                bam=expand("rmduped_reads/{sample}.bowtie2.sorted.rmdup.bam",sample=config["samples"]),
                bai=expand("rmduped_reads/{sample}.bowtie2.sorted.rmdup.bam.bai",sample=config["samples"])
        output:
                "calls/all_samples.P.nigra.bowtie2.freebayes_populations.vcf"
        threads: 20
        log:
                "logs/freebayes/all_samples.log"
        benchmark:
                "benchmarks/all_samples.freebayes.benchmark.txt"
        shell:
                "module load freebayes; freebayes -f {input.ref} {input.bam} > {output} 2> freebayes.error"

rule report:
	input:
		T1="calls/all_samples.P.nigra.bowtie2.freebayes_populations.vcf",
		T2=expand("benchmarks/{sample}.bowtie2.benchmark.txt",sample=config["samples"])
	output:
		"freebayes_calls.report.html"
	run:
		from snakemake.utils import report
		with open(input.T1) as vcf:
			n_calls = sum(1 for l in vcf if not l.startswith("#"))

		report("""
		A first test of a structural variant calling workflow
		=====================================================
		
		Reads were mapped to the P.nigra
		reference genome and variants were called jointly with
		gatk unified genotyper.

		This resulted in {n_calls} variants (see Table T1_).
		Benchmark results for BWA can be found in the tables T2_.
		""", output[0], T1=input[0])
