
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
		"benchmarks/{sample}-{unit}.hisat2.benchmark.txt"
	shell:
		"module load hisat2; module load samtools; hisat2 -p 4 --mm --time -x {params.ref} {params.rg} -1 {input.fq1} -2 {input.fq2} | samtools view -bS - > {output} 2> {log}"

rule mosaik_align:
    input:
        get_trimmed
    output:
        temp("mapped_reads/{sample}-{unit}.bam")
    log:
        "logs/mosaid_align/{sample}-{unit}.log"
    params:
        ref=get_reference
    benchmark:
        "benchmarks/{sample}-{unit}.mosaik_align.benchmark.txt"
    shell:
        "module load mosaik; module load samtools; MosaikAligner -q {input[0]} -q2 {input[1]} -out {output} -ia {params.ref}"
