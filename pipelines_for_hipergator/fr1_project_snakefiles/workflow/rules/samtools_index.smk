rule samtools_index_rmduped:
	input:
		"results/rmduped/{sample}.sorted.rmdup.bam"
	output:
		"results/rmduped/{sample}.sorted.rmdup.bam.bai"
	log:
		"logs/samtools_index_sorted.{sample}.log"
	benchmark:
		"benchmarks/{sample}.index_sorted.benchmark.txt"
	shell:
		"module load samtools; samtools index {input}"

rule samtools_index_realigned:
	input:
		"results/realigned/{sample}.sorted.rmdup.realigned.bam"		
	output:
		"results/realigned/{sample}.sorted.rmdup.realigned.bam.bai"
	log:
		"logs/samtools_index_realigned.{sample}.log"
	benchmark:
		"benchmarks/{sample}.index_realigned.benchmark.txt"
	shell:
		"module load samtools; samtools index {input}"
