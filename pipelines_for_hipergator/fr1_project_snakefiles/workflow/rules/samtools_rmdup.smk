rule samtools_rmdup:
	input:
		"results/sorted/{sample}.sorted.bam"
	output:
		"results/rmduped/{sample}.sorted.rmdup.bam"
	log:
		"logs/samtools_rmdup/{sample}.log"
	benchmark:
		"benchmarks/{sample}.rmdup.benchmark.txt"
	shell:
		"module load samtools; samtools rmdup {input} {output}"
