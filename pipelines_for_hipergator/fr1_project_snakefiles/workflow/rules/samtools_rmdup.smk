rule samtools_rmdup:
	input:
		"results/sorted/{sample}-{unit}.sorted.bam"
	output:
		"results/rmduped/{sample}-{unit}.sorted.rmdup.bam"
	log:
		"logs/samtools_rmdup/{sample}-{unit}.log"
	benchmark:
		"benchmarks/{sample}-{unit}.rmdup.benchmark.txt"
	shell:
		"module load samtools; samtools rmdup {input} {output}"
