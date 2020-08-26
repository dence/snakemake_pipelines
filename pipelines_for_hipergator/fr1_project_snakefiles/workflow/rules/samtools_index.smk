rule samtools_index_rmduped:
	input:
		"rmduped{sample}-{unit}.sorted.rmdup.bam"
	output:
		"rmduped/{sample}-{unit}.sorted.rmdup.bam.bai"
	log:
		"logs/samtools_index_sorted.{sample}-{unit}.log"
	benchmark:
		"benchmarks/{sample}-{unit}.index_sorted.benchmark.txt"
	shell:
		"module load samtools; samtools index {input}"

rule samtools_index_realigned:
	input:
		"realigned/{sample}-{unit}.sorted.rmdup.realigned.bam"
	output:
		"realigned/{sample}-{unit}.sorted.rmdup.realigned.bam.bai"
	log:
		"logs/samtools_index_realigned.{sample}-{unit}.log"
	benchmark:
		"benchmarks/{sample}-{unit}.index_realigned.benchmark.txt"
	shell:
		"module load samtools; samtools index {input}"
