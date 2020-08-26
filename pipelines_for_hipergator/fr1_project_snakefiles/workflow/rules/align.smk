rule bwa_mem_align:
	input:
		get_trimmed
	params:
		sample="{sample}-{unit}",
		ref=get_reference
	output:
		temp("results/mapped_bwamem/{sample}-{unit}.bam")
	log:
		"logs/bwa_mem/{sample}-{unit}.bwa_mem.log"
	shell:
		"module load bwa/0.7.17 ; bwa mem {params.ref} "
		+ "-R '@RG\\tID:{params.sample}\\tSM:{params.sample}\\tPL:ILLUMINA'"
		+ "{input[0]} {input[1]} > {output} 2> {log}"


rule hisat2_align:
	input:
		get_trimmed
	output:
		temp("results/mapped_hisat2/{sample}-{unit}.bam")
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
		temp("results/mapped_mosaik/{sample}-{unit}.bam"),
		temp("results/mosaik_temp/{sample}-{unit}.mosaik")
	log:
		"logs/mosaik_align/{sample}-{unit}.log"
	params:
		ref=get_mosaik_reference,
		seq=config["settings"]["mosaik_seq_param"]
	benchmark:
		"benchmarks/{sample}-{unit}.mosaik_align.benchmark.txt"
	shell:
		"module load mosaik; MosaikBuild -q {input[0]} -q2 {input[1]} -out {output[1]} -st {params.seq}; MosaikAligner -in {output[1]} -out {output[0]} -ia {params.ref}"
