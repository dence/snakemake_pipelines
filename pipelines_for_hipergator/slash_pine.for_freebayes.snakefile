rule all:
	input:
		"freebayes_calls.report.html"

rule bwa_map:
	input:
		config["reference"]["V1_01"],	
		#"/ufrc/kirst/share/Genomes/Ptaeda_v101/bwa/ptaeda.v1.0.scaffolds.fasta",	
		"/ufrc/kirst/share/slash_pine_CFRGP/RAPiD-Genomics_HG5Y5BBXX_MAT_068502_P01_{sample}_i5-501_i7-27_S2_L001_R1_001.fastq.gz",
		"/ufrc/kirst/share/slash_pine_CFRGP/RAPiD-Genomics_HG5Y5BBXX_MAT_068502_P01_{sample}_i5-501_i7-27_S2_L001_R2_001.fastq.gz"
	output:
		"mapped_reads/{sample}.bam"
	params:
		rg="@RG\tID:{sample}\tSM:{sample}"
	log:
		"logs/bwa_mem/{sample}.log"
	benchmark:
		"benchmarks/{sample}.bwa.benchmark.txt"
	threads: 4
	shell:
		"module load bwa; bwa mem -R \"{params.rg}\" -t {threads} {input} | samtools view -Sb > {output} 2> {log}"

rule samtools_sort:
	input:
		"mapped_reads/{sample}.bam"
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
		"module load samtools; samtools rmdup {input} {output} "

rule samtools_index_rmduped:
	input:
		"rmduped_reads/{sample}.sorted.rmdup.bam"
	output:
		"rmduped_reads/{sample}.sorted.rmdup.bam.bai"
	log:
		"logs/smatools_index_sorted{sample}.log"
	benchmark:
		"benchmarks/{sample}.index_sorted.benchmark.txt"
	shell:
		"module load samtools; samtools index {input} "
		

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
		"module load samtools; samtool index {input} ;"

rule freebayes:
	input:	
		ref=config["reference"]["V1_01"],
		bam=expand("rmduped_reads/{sample}.sorted.rmdup.bam",sample=config["samples"])
	output:
		"calls/all_samples.pinus_taeda_L.v1_1.freebayes_populations.vcf"
	threads: 4
	log:
		"logs/freebayes/all_samples.log"
	benchmark:
		"benchmarks/all_samples.freebayes.benchmark.txt"
	shell:
		"module load freebayes; freebayes -f {input.ref} {input.bam} > {output}"

rule report:
	input:
		T1="calls/all_samples.pinus_taeda_L.v1_1.freebayes_populations.vcf",
		T2=expand("benchmarks/{sample}.bwa.benchmark.txt",sample=config["samples"])
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
