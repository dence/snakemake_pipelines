rule samtools_sort:
	input:
		"results/RG_replaced_bams/{sample}.merged.bam"
	output:
		temp("results/sorted/{sample}.sorted.bam")
	log:
		"logs/samtools_sort/{sample}.log"
	benchmark:
		"benchmarks/{sample}.sort.benchmark.txt"
	params:
		"-T sorted_reads/{sample}"
	threads: 4
	shell:
		"unset TMPDIR; module load samtools; samtools sort --threads {threads} {params} -O bam {input} > {output}"


#def get_bams(wildcards):
#	return expand("results/mapped_bwamem/{sample}-{unit}.bam",sample=units["sample"],unit=units["unit"])
