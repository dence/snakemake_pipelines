rule all:
	input:
		"freebayes_calls.report.html"

rule pretrim_fastqc:
	input:
		fq1="/orange/kirst/d.ence/pinus_taeda_L/fr1_sequencing/RAPiD-Genomics_HJYNGBBXX_UFL_104301_{sample}_R1_001.fastq.gz",
		fq2="/orange/kirst/d.ence/pinus_taeda_L/fr1_sequencing/RAPiD-Genomics_HJYNGBBXX_UFL_104301_{sample}_R2_001.fastq.gz"
	output:
		out="fastqc_out/pre_process/{sample}/",
		log_file="logs/fastqc_pre_out/{sample}.log"
	log:
		"logs/fastqc_pre_out/{sample}.log"
	benchmark:
		"benchmarks/{sample}.fastqc.benchmark.txt"
	shell:
		"module load fastqc ; fastqc -f fastq -o {output.out} {input} &> {log}"
rule trim_reads:
	input:
		fq1="/orange/kirst/d.ence/pinus_taeda_L/fr1_sequencing/RAPiD-Genomics_HJYNGBBXX_UFL_104301_{sample}_R1_001.fastq.gz",
		fq2="/orange/kirst/d.ence/pinus_taeda_L/fr1_sequencing/RAPiD-Genomics_HJYNGBBXX_UFL_104301_{sample}_R2_001.fastq.gz",
		prev_log="logs/fastqc_pre_out/{sample}.log"
	output:
		fq1="trimmed_reads/{sample}/{sample}.trimmed.R1.fastq.gz",
		fq2="trimmed_reads/{sample}/{sample}.trimmed.R2.fastq.gz"
	log:
		"logs/cutadapt/{sample}.log"
	benchmark:
		"benchmarks/{sample}.fastqc.benchmark.txt"
	shell:
		"module load python/2.7.14; module load gcc/5.2.0; module load cutadapt/1.18; cutadapt -u 10 -U 10 -q30,30 --minimum-length=50 -o {output.fq1} -p {output.fq2} {input.fq1} {input.fq2} &> {log}"
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
		sample="{sample}",
		ref=config["reference"]["V1_01"]["custom"]
	output:
		"mapped_reads/{sample}.bwa_mem.bam"
	log:
		"logs/bwa_mem/{sample}.bwa_mem.log"
	shell:
		"module load bwa/0.7.17 ; bwa mem {params.ref} " 
		+ "-R '@RG\\tID:{params.sample}\\tSM:{params.sample}\\tPL:{params.sample}'" 
		+ "{input.fq1} {input.fq2} > {output} 2> {log}"	

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
		bam="rmduped_reads/{sample}.bwa_mem.sorted.rmdup.bam",
		bai="rmduped_reads/{sample}.bwa_mem.sorted.rmdup.bam.bai"
	output:
		"realigner_intervals/{sample}.bwa_mem.intervals"
	params:
		ref=config["reference"]["V1_01"]["custom"]
	log:
		"logs/realigner_intervals/{sample}.intervals.log"
	benchmark:
		"benchmarks/{sample}.intervals.benchmark.txt"
	shell:
		"module load gatk;  java -jar -Xmx10g /apps/gatk/3.7.0/GenomeAnalysisTK.jar -T RealignerTargetCreator -I {input.bam} -o {output} -R {params.ref}"
rule gatk_indel_realign:
	input:
		bam="rmduped_reads/{sample}.bwa_mem.sorted.rmdup.bam",
		bai="rmduped_reads/{sample}.bwa_mem.sorted.rmdup.bam.bai",
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
		L4_bam="realigned_bams/{sample}_L004.bwa_mem.sorted.rmdup.realigned.bam",
		L5_bam="realigned_bams/{sample}_L005.bwa_mem.sorted.rmdup.realigned.bam"
	output:
		"merged_lane_bams/{sample}.bwa_mem.sorted.rmdup.realigned.merged.bam"
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
		"merged_lane_bams/{sample}.bwa_mem.sorted.rmdup.realigned.merged.bam"
	output:
		"RG_replaced_bams/{sample}.bwa_mem.sorted.rmdup.realigned.merged.bam"
	params:
		RG_fields="RGID={sample} RGLB={sample} RGPL=illumina RGPU={sample} RGSM={sample}"
	log:
		"logs/picard_replaceRG/{sample}.log"
	shell:
		"module load picard; java -jar $HPC_PICARD_DIR/picard.jar AddOrReplaceReadGroups I={input} O={output} {params.RG_fields} &> {log}"

rule samtools_index_merged:
        input:
                "RG_replaced_bams/{sample}.bwa_mem.sorted.rmdup.realigned.merged.bam"
        output:
                "RG_replaced_bams/{sample}.bwa_mem.sorted.rmdup.realigned.merged.bam.bai"
        log:
                "logs/samtools_index_merged/{sample}.log"
        benchmark:
                "benchmarks/{sample}.index_merged.benchmark.txt"
        shell:
                "module load samtools; samtools index {input}"

rule freebayes_haploid:
	input:
		bam=expand("RG_replaced_bams/{sample}.bwa_mem.sorted.rmdup.realigned.merged.bam",sample=config["Fr1_meg_samples"]),
		bai=expand("RG_replaced_bams/{sample}.bwa_mem.sorted.rmdup.realigned.merged.bam.bai",sample=config["Fr1_meg_samples"]),
		bam_list="/home/d.ence/projects/pinus_taeda_L/Fr1_project/aligning_Fr1_samples/aligning_to_custom_R_gene_V1_1_ref/calls/bam_list.megs.txt",
		ref=config["reference"]["V1_01"]["custom"]
	log:
		"logs/freebayes/Fr1_megs.freebayes.log"
	output:
		"calls/freebayes/Fr1_megs.bwa_mem.freebayes.vcf"
	shell:
		"module load freebayes; freebayes-v1.3.1 --skip-coverage 200 --ploidy 2  --min-coverage 50 --bam-list {input.bam_list} -f {input.ref} --vcf {output} --targets /home/d.ence/projects/pinus_taeda_L/Fr1_project/probe_files/RG_0610_merged.fixed.bed &> {log}"

rule freebayes_10_5_prog:
        input:
                bam=expand("RG_replaced_bams/{sample}.bwa_mem.sorted.rmdup.realigned.merged.bam",sample=config["Fr1_prog_samples"]),
                bai=expand("RG_replaced_bams/{sample}.bwa_mem.sorted.rmdup.realigned.merged.bam.bai",sample=config["Fr1_prog_samples"]),
                bam_list="/home/d.ence/projects/pinus_taeda_L/Fr1_project/aligning_Fr1_samples/aligning_to_custom_R_gene_V1_1_ref/calls/bam_list.Fr1_prog.txt",
                ref=config["reference"]["V1_01"]["custom"]
        log:
                "logs/freebayes/Fr1_prog.freebayes.log"
        output:
                "calls/freebayes/Fr1_prog.bwa_mem.freebayes.vcf"
        shell:
                "module load freebayes; freebayes-v1.3.1 --ploidy 2 --skip-coverage 200 --min-coverage 50 --bam-list {input.bam_list} -f {input.ref} --vcf {output} --targets /home/d.ence/projects/pinus_taeda_L/Fr1_project/probe_files/RG_0610_merged.fixed.bed &> {log}"

rule freebayes_non_10_5_indv:
        input:
                bam=expand("RG_replaced_bams/{sample}.bwa_mem.sorted.rmdup.realigned.merged.bam",sample=config["non_10_5_indv"]),
                bai=expand("RG_replaced_bams/{sample}.bwa_mem.sorted.rmdup.realigned.merged.bam.bai",sample=config["non_10_5_indv"]),
                bam_list="/home/d.ence/projects/pinus_taeda_L/Fr1_project/aligning_Fr1_samples/aligning_to_custom_R_gene_V1_1_ref/calls/bam_list.non_10_5_indv.txt",
                ref=config["reference"]["V1_01"]["custom"]
        log:
                "logs/freebayes/non_10_5_indv.freebayes.log"
        output:
                "calls/freebayes/non_10_5_indv.bwa_mem.freebayes.vcf" 
        shell:
                "module load freebayes; freebayes-v1.3.1 --ploidy 2 --skip-coverage 200 --min-coverage 50 --bam-list {input.bam_list} -f {input.ref} --vcf {output} --targets /home/d.ence/projects/pinus_taeda_L/Fr1_project/probe_files/RG_0610_merged.fixed.bed &> {log}"

rule freebayes_non_10_5_pooled_single_sample:
	input:
                bam="RG_replaced_bams/{sample}.bwa_mem.sorted.rmdup.realigned.merged.bam",
                bai="RG_replaced_bams/{sample}.bwa_mem.sorted.rmdup.realigned.merged.bam.bai",
                ref=config["reference"]["V1_01"]["custom"]
	output:
		vcf="calls/freebayes/{sample}.non_10_5_pooled.bwa_mem.freebayes.vcf"
	log:
                "logs/freebayes/{sample}.non_10_5_pooled.freebayes.log"
	shell:
                "module load freebayes; freebayes-v1.3.1 --ploidy 20 --skip-coverage 200 --use-best-n-alleles 2 --pooled-discrete  --min-coverage 50 --bam {input.bam} -f {input.ref} --vcf {output.vcf} --targets /home/d.ence/projects/pinus_taeda_L/Fr1_project/probe_files/RG_0610_merged.fixed.bed &> {log}"

rule bgzip_and_tabix_vcf:
	input:
		vcf="calls/freebayes/{sample}.non_10_5_pooled.bwa_mem.freebayes.vcf"
	output:
		bgzipped="calls/freebayes/{sample}.non_10_5_pooled.bwa_mem.freebayes.vcf.gz",
		tabix_indexed="calls/freebayes/{sample}.non_10_5_pooled.bwa_mem.freebayes.gz.tbi"
	shell:
		"module load vcflib; cat {input.vcf} | bgzip -c > {output.bgzipped} ; tabix -p vcf {output.bgzipped} ; touch {output.tabix_indexed}"

rule merge_freebayes_non_10_5_pooled_singles:
	input:
		bgzipped_vcf=expand("calls/freebayes/{sample}.non_10_5_pooled.bwa_mem.freebayes.vcf.gz",sample=config["non_10_5_pooled"])
	output:
		vcf="calls/freebayes/non_10_5_pooled.bwa_mem.freebayes.vcf"
	log:
		"logs/freebayes/non_10_5_pooled.freebayes.log"
	shell:
		"module load vcftools; vcf-merge -d --ref-for-missing 0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0 {input.bgzipped_vcf} > {output.vcf} 2> {log}"

rule mosdepth:
	input: bam="RG_replaced_bams/{sample}.bwa_mem.sorted.rmdup.realigned.merged.bam",
		bai="RG_replaced_bams/{sample}.bwa_mem.sorted.rmdup.realigned.merged.bam.bai",
		ref=config["reference"]["V1_01"]["custom"]
	params:
		" --fast-mode -b /home/d.ence/projects/pinus_taeda_L/Fr1_project/probe_files/RG_0610_merged.fixed.bed --no-per-base {sample}"
	output:
		summary="mosdepth/{sample}.mosdepth.summary.txt"
	log:
		"logs/mosdepth/{sample}.mosdepth.log"
	shell:
		"time mosdepth {params} {input.bam} "	

#rule freebayes_non_10_5_pooled:
#        input:
#                bam=expand("RG_replaced_bams/{sample}.bwa_mem.sorted.rmdup.realigned.merged.bam",sample=config["non_10_5_pooled"]),
#                bai=expand("RG_replaced_bams/{sample}.bwa_mem.sorted.rmdup.realigned.merged.bam.bai",sample=config["non_10_5_pooled"]),
#                bam_list="/home/d.ence/projects/pinus_taeda_L/Fr1_project/aligning_Fr1_samples/aligning_to_custom_R_gene_V1_1_ref/calls/bam_list.non_10_5_pooled.txt",
#                ref=config["reference"]["V1_01"]["custom"]
#        log:
#                "logs/freebayes/non_10_5_pooled.freebayes.log"
#        output:
#                "calls/freebayes/non_10_5_pooled.bwa_mem.freebayes.vcf"
#        shell:
#                "module load freebayes; freebayes-v1.3.1 --ploidy 20  --use-best-n-alleles 20 --pooled-discrete  --min-coverage 50 --bam-list {input.bam_list} -f {input.ref} --vcf {output} --targets /home/d.ence/projects/pinus_taeda_L/Fr1_project/probe_files/RG_0610_merged.fixed.bed &> {log}"

rule report:
	input:
		megs_vcf="calls/freebayes/Fr1_megs.bwa_mem.freebayes.vcf",
		prog_vcf="calls/freebayes/Fr1_prog.bwa_mem.freebayes.vcf",
		non_10_5_indv_vcf="calls/freebayes/non_10_5_indv.bwa_mem.freebayes.vcf",
		non_10_5_pooled_vcf="calls/freebayes/non_10_5_pooled.bwa_mem.freebayes.vcf",
		mosdepth_files="mosdepth/{sample}.summary.txt"
	output:
		"freebayes_calls.report.html"
	run:
		from snakemake.utils import report
		with open(input.megs_vcf) as vcf:
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



