
rule hisat2_align:
	input:
		get_trimmed
	output:
		temp("mapped_reads/{sample}-{unit}.bam")
	log:
		"logs/hisat2/{sample}-{unit}.log"
	params:
		rg="--rg-id={sample}-{unit} --rg \"SM:{sample}-{unit}\" --rg \"LB:{sample}-{unit}\" --rg \"PL:ILLUMINA\"",
        ref=get_reference
        #ref=config["resources"]["reference"]["V1_01"]["custom"]
		#ref=config["reference"]["V2_01"]["custom"]
	benchmark:
		"benchmarks/{sample}.hisat2.benchmark.txt"
	shell:
		"module load hisat2; module load samtools; hisat2 -p 4 --mm --time -x {params.ref} {params.rg} -1 {input.fq1} -2 {input.fq2} | samtools view -bS - > {output} 2> {log}"
