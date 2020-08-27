rule samtools_sort_bwa_mem:
	input:
        "results/RG_replaced_bams/{sample}.bwa_mem.sorted.rmdup.merged.bam"
	output:
		temp("results/sorted/{sample}.sorted.bam")
	log:
		"logs/samtools_sort/{sample}.log"
	benchmark:
		"benchmarks/{sample}.sort.benchmark.txt"
	threads: 4
	shell:
		"module load samtools; samtools sort --threads {threads} -T sorted_reads/{sample}-{unit} -O bam {input} > {output}"


#def get_bams(wildcards):
#	return expand("results/mapped_bwamem/{sample}-{unit}.bam",sample=units["sample"],unit=units["unit"])
