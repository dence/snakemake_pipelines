#Daniel Ence
#08/23/2020


rule pretrim_fastqc:
	input:
		fq1="/orange/kirst/d.ence/pinus_taeda_L/fr1_sequencing/RAPiD-Genomics_HJYNGBBXX_UFL_104301_{long_sample}_{unit}_R1_001.fastq.gz",
		fq2="/orange/kirst/d.ence/pinus_taeda_L/fr1_sequencing/RAPiD-Genomics_HJYNGBBXX_UFL_104301_{long_sample}_{unit}_R2_001.fastq.gz"
	output:
		out="results/fastqc_out/pre_process/{sample}-{unit}/",
		log_file="results/logs/fastqc_pre_out/{sample}-{unit}.log"
	log:
		"results/logs/fastqc_pre_out/{sample}-{unit}.log"
	benchmark:
		"benchmarks/{sample}-{unit}.fastqc.benchmark.txt"
	shell:
		"module load fastqc ; fastqc -f fastq -o {output.out} {input} &> {log}"

rule posttrim_fastqc:
	input:
		"results/trimmed_reads/{sample}-{unit}/{sample}-{unit}.trimmed.R1.fastq.gz",
		"results/trimmed_reads/{sample}-{unit}/{sample}-{unit}.trimmed.R2.fastq.gz"
	output:
		"results/fastqc_out/post_trim/{sample}-{unit}/",
		"results/logs/fastqc_post_out/{sample}-{unit}.log"
	log:
		"results/logs/fastqc_post_out/{sample}-{unit}.log"
	benchmark:
		"benchmarks/{sample}-{unit}.fastqc.benchmark.txt"
	shell:
		"module load fastqc ; fastqc -f fastq -o {output} {input} &> {log}"
