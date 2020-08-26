
#Daniel Ence, 08/20/2020
rule trim_reads:
	input:
		get_fastqs
	output:
		fq1=temp("results/trimmed/{sample}-{unit}.trimmed.R1.fastq.gz"),
		fq2=temp("results/trimmed/{sample}-{unit}.trimmed.R2.fastq.gz")
	log:
		"results/logs/cutadapt/{sample}-{unit}.log"
	benchmark:
		"results/benchmarks/{sample}-{unit}.fastqc.benchmark.txt"
	shell:
		"module load python/2.7.14; module load gcc/5.2.0; module load cutadapt; cutadapt -u 10 -U 10 -q30,30 --minimum-length=50 -o {output.fq1} -p {output.fq2} {input[0]} {input[1]} &> {log}"
