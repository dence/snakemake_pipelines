#Daniel Ence
#08/23/2020


rule pretrim_fastqc:
	input:
		get_fastqs
	output:
		out=directory("results/fastqc_out/pre_process/{sample}-{unit}/"),
		log_file="results/logs/fastqc_pre_out/{sample}-{unit}.log"
	log:
		"results/logs/fastqc_pre_out/{sample}-{unit}.log"
	benchmark:
		"benchmarks/{sample}-{unit}.fastqc.benchmark.txt"
	shell:
		"module load fastqc ; fastqc -f fastq -o {output.out} {input} &> {log}"

rule posttrim_fastqc:
	input:
		get_trimmed
	output:
		directory("results/fastqc_out/post_trim/{sample}-{unit}/"),
		"results/logs/fastqc_post_out/{sample}-{unit}.log"
	log:
		"results/logs/fastqc_post_out/{sample}-{unit}.log"
	benchmark:
		"benchmarks/{sample}-{unit}.fastqc.benchmark.txt"
	shell:
		"module load fastqc ; fastqc -f fastq -o {output[0]} {input} &> {log}"
