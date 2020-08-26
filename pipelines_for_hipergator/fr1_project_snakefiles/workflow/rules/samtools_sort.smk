rule samtools_sort_bwa_mem:
	input:
        "results/mapped_bwamem/{sample}-{unit}.bam"
	output:
		temp("results/sorted/{sample}-{unit}.sorted.bam")
	log:
		"logs/samtools_sort/{sample}-{unit}.log"
	benchmark:
		"benchmarks/{sample}-{unit}.sort.benchmark.txt"
	threads: 4
	shell:
		"module load samtools; samtools sort --threads {threads} -T sorted_reads/{sample}-{unit} -O bam {input} > {output}"
