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
		"unset TMPDIR; module load samtools; samtools index {input}"

rule samtools_index_RG_replaced:
	input:
		"results/RG_replaced_bams/{sample}.merged.bam"
	output:
		"results/RG_replaced_bams/{sample}.merged.bam.bai"
	log:
		"logs/samtools_index_RG_replaced.{sample}.log"
	benchmark:
		"benchmarks/{sample}.index_RG_replaced.benchmark.txt"
	shell:
		"unset TMPDIR; module load samtools; samtools index {input}"

rule samtools_index_realigned:
	input:
		"results/realigned/{sample}.realigned.bam"
	output:
		"results/realigned/{sample}.realigned.bam.bai"
	log:
		"logs/samtools_index_realigned.{sample}.log"
	benchmark:
		"benchmarks/{sample}.index_realigned.benchmark.txt"
	shell:
		"unset TMPDIR; module load samtools; samtools index {input}"
