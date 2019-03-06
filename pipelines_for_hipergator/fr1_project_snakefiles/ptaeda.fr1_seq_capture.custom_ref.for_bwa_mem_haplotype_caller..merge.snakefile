rule all:
	input:
		"freebayes_calls.report.html"
rule pretrim_fastqc:
	input:
		#fq1="/ufrc/kirst/d.ence/pinus_taeda_L/fr1_sequencing//RAPiD-Genomics_HJYNGBBXX_UFL_104301_{sample}_R1_001.fastq.gz",
		fq1="/orange/kirst/d.ence/pinus_taeda_L/fr1_sequencing/RAPiD-Genomics_HJYNGBBXX_UFL_104301_{sample}_R1_001.fastq.gz",
		#fq2="/ufrc/kirst/d.ence/pinus_taeda_L/fr1_sequencing//RAPiD-Genomics_HJYNGBBXX_UFL_104301_{sample}_R2_001.fastq.gz"
		fq2="/orange/kirst/d.ence/pinus_taeda_L/fr1_sequencing/RAPiD-Genomics_HJYNGBBXX_UFL_104301_{sample}_R2_001.fastq.gz"
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
		fq1="/orange/kirst/d.ence/pinus_taeda_L/fr1_sequencing/RAPiD-Genomics_HJYNGBBXX_UFL_104301_{sample}_R1_001.fastq.gz",
		fq2="/orange/kirst/d.ence/pinus_taeda_L/fr1_sequencing/RAPiD-Genomics_HJYNGBBXX_UFL_104301_{sample}_R2_001.fastq.gz"
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
rule bwamem_align:
	input:
		fq1="trimmed_reads/{sample}/{sample}.trimmed.R1.fastq.gz",
		fq2="trimmed_reads/{sample}/{sample}.trimmed.R2.fastq.gz"
	params:
		RG="\"\@RG\tID:{sample}\tSM:{sample}\"",
		ref=config["reference"]["V1_01"]["custom"]
	output:
		"mapped_reads/{sample}.bwa_mem.bam"
	log:
		"logs/bwa_mem/{sample}.bwa_mem.log"
	shell:
		"module load bwa/0.7.17 ; bwa mem {config ref} -R {params.RG} {input.fq1} {input.fq2} > {output} 2> {log}"	

rule samtools_sort:
	input:
		"mapped_reads/{sample}.bwa_mem.bam"
	output:
		"sorted_reads/{sample}.bwa_mem.sorted.bam"
	log:
		"logs/samtools_sort/{sample}.log"
	benchmark:
		"benchmarks/{sample}.sort.benchmark.txt"
	threads: 4		
	shell:
		"module load samtools; samtools sort --threads {threads} -T sorted_reads/{wildcards.sample} -O bam {input} > {output}"

rule samtools_rmdup:
	input:
		"sorted_reads/{sample}.bwa_mem.sorted.bam"
	output:
		"rmduped_reads/{sample}.bwa_mem.sorted.rmdup.bam"
	log:
		"logs/samtools_rmdup/{sample}.log"
	benchmark:
		"benchmarks/{sample}.rmdup.benchmark.txt"
	shell:
		"module load samtools; samtools rmdup {input} {output}"

rule samtools_index_rmduped:
	input:
		"rmduped_reads/{sample}.bwa_mem.sorted.rmdup.bam"
	output:
		"rmduped_reads/{sample}.bwa_mem.sorted.rmdup.bam.bai"
	log:
		"logs/samtools_index_sorted.{sample}.log"
	benchmark:
		"benchmarks/{sample}.index_sorted.benchmark.txt"
	shell:
		"module load samtools; samtools index {input}"

rule samtools_index_sorted:
	input:
		"sorted_reads/{sample}.bwa_mem.sorted.rmdup.bam"
	output:
		"sorted_reads/{sample}.bwa_mem.sorted.rmdup.bam.bai"
	log:
		"logs/samtools_index_sorted/{sample}.log"
	benchmark:
		"benchmarks/{sample}.index_sorted.benchmark.txt"	
	shell:
		"module load samtools; samtool index {input}"
rule gatk_indel_creator:
	input:
		"rmduped_reads/{sample}.bwa_mem.sorted.rmdup.bam"
	output:
		"realigner_intervals/{sample}.bwa_mem.intervals"
	params:
		ref=config["reference"]["V1_01"]["custom"]
	log:
		"logs/realigner_intervals/{sample}.intervals.log"
	benchmark:
		"benchmarks/{sample}.intervals.benchmark.txt"
	shell:
		"module load gatk;  java -jar -Xmx4g /apps/gatk/3.7.0/GenomeAnalysisTK.jar -T RealignerTargetCreator -I {input} -o {output} -R {params.ref}"
rule gatk_indel_realign:
	input:
		bam="rmduped_reads/{sample}.bwa_mem.sorted.rmdup.bam",
		interval="realigner_intervals/{sample}.bwa_mem.intervals"
	output:
		"realigned_bams/{sample}.bwa_mem.sorted.rmdup.realigned.bam"
	params:
		ref=config["reference"]["V1_01"]["custom"]
	log:
		"logs/realigner/{sample}.realigner.log"
	benchmark:
		"benchmarks/{sample}.intervals.benchmark.txt"
	shell:
		"module load gatk; java -jar -Xmx4g /apps/gatk/3.7.0/GenomeAnalysisTK.jar -T IndelRealigner -R {params.ref} -I {input.bam} -targetIntervals {input.interval} -o {output}"
rule samtools_index_realigned:
	input:
		"realigned_bams/{sample}.bwa_mem.sorted.rmdup.realigned.bam"
	output:
		"realigned_bams/{sample}.bwa_mem.sorted.rmdup.realigned.bai"
	log:
		"logs/samtools_index_realigned/{sample}.log"
	benchmark:
		"benchmarks/{sample}.index_realigned.benchmark.txt"
	shell:
		"module load samtools; samtool index {input}"

rule freebayes:
        input:
                ref=config["reference"]["V1_01"]["custom"],
                bam=expand("realigned_bams/{sample}.bwa_mem.sorted.rmdup.realigned.bam",sample=config["samples"]),
                bai=expand("realigned_bams/{sample}.bwa_mem.sorted.rmdup.realigned.bai",sample=config["samples"])
                #bam=expand("rmduped_reads/{sample}.hisat2.sorted.rmdup.bam",sample=config["samples"]),
                #bai=expand("rmduped_reads/{sample}.hisat2.sorted.rmdup.bam.bai",sample=config["samples"])
        output:
                "calls/all_samples.bwa_mem.freebayes_populations.vcf"
        threads: 20
        log:
                "logs/freebayes/all_samples.log"
        benchmark:
                "benchmarks/all_samples.freebayes.benchmark.txt"
        shell:
                "module load freebayes; freebayes -f {input.ref} {input.bam} > {output} 2> freebayes.error"

rule merge_lane_bams:
	input:
		#making a dumb assumption about the names of the bams to merged. specific to the Fr1 project. DE
		L4_bam="rmduped_reads/{sample}_L004.bwa_mem.sorted.rmdup.bam",
		L5_bam="rmduped_reads/{sample}_L005.bwa_mem.sorted.rmdup.bam"
	output:
		"merged_lane_bams/{sample}.bwa_mem.sorted.rmdup.merged.bam"
	log:
		"logs/picard_merge_sam_files/{sample}.log"
	shell:
		"module load picard; java -jar $HPC_PICARD_DIR/picard.jar MergeSamFiles I={input.L4_bam} I={input.L5_bam} O={output} &> {log}"

rule index_merged_bams:
	input:
		"merged_lane_bams/{sample}.bwa_mem.sorted.rmdup.merged.bam"
	output:
		"merged_lane_bams/{sample}.bwa_mem.sorted.rmdup.merged.bam.bai"
	shell:
		"module load samtools; samtools index {input}"

rule samtools_index_replaced:
	input:
		"RG_replaced_bams/{sample}.bwa_mem.sorted.rmdup.merged.bam"
	output:
		"RG_replaced_bams/{sample}.bwa_mem.sorted.rmdup.merged.bam.bai"
	log:
		"logs/samtools_index_realigned/{sample}.log"
	benchmark:
		"benchmarks/{sample}.index_realigned.benchmark.txt"
	shell:
		"module load samtools; samtools index {input}"

rule ReplaceRG_merged:
	input:
		"merged_lane_bams/{sample}.bwa_mem.sorted.rmdup.merged.bam"
	output:
		"RG_replaced_bams/{sample}.bwa_mem.sorted.rmdup.merged.bam"
	params:
		RG_fields="RGID={sample} RGLB={sample} RGPL=illumina RGPU={sample} RGSM={sample}"
	log:
		"logs/picard_replaceRG/{sample}.log"
	shell:
		"module load picard; java -jar $HPC_PICARD_DIR/picard.jar AddOrReplaceReadGroups I={input} O={output} {params.RG_fields} &> {log}"

rule haplotype_caller_unique_samples:
	input:
		bam="RG_replaced_bams/{sample}.bwa_mem.sorted.rmdup.merged.bam",
		bai="RG_replaced_bams/{sample}.bwa_mem.sorted.rmdup.merged.bam.bai",
		#bam="merged_lane_bams/{sample}.hisat2.sorted.rmdup.merged.bam",
		#bai="merged_lane_bams/{sample}.hisat2.sorted.rmdup.merged.bam.bai",
		ref=config["reference"]["V1_01"]["custom"]
	output:
		"calls/indv_HaplCalls/{sample}.HC_individual.gvcf"
	shell:
		"module load gatk/3.7.0; java -jar $HPC_GATK_DIR/GenomeAnalysisTK.jar -R {input.ref} -T HaplotypeCaller -I {input.bam}  -rf BadCigar -o {output}"
rule combine_GVCFs:
	input:
		gvcf=expand("calls/indv_HaplCalls/{sample}.HC_individual.gvcf",sample=config["samples"]),
		ref=config["reference"]["V1_01"]["custom"]
	output:
		"calls/all_samples.combined.haplotypecaller.gvcf"
	shell:
		"module load gatk/3.7.0; java -jar $HPC_GATK_DIR/GenomeAnalysisTK.jar -R {input.ref} --variant {input} -o {output}"
rule report:
	input:
		T1="calls/all_samples.combined.haplotypecaller.gvcf"
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
		Benchmark results for BWA can be found in the tables.
		""", output[0], T1=input[0])



