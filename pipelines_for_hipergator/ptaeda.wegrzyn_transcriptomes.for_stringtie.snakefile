rule all:
	input:
		"freebayes_calls.report.html"

rule gmap:
	input:
		fasta="/home/d.ence/projects/pinus_taeda_L/wegrzyn_elite_transcriptome/trinityassemblies/{sample}_transcriptassembly.fasta",
	output:
		"mapped_reads/{sample}.bam"
	params:
		rg='--read-group={sample} --read-group-name={sample}={sample} --read-group-library={sample} --read-group-platform=trinity',
		ref_param="-D /home/d.ence/projects/pinus_taeda_L/aligning_unigenes_w_gmap/CCLONES_all_seq/CCLONES_all_seq/ -d CCLONES_all_seq" 
	log:
		"logs/gmap/{sample}.log"
	benchmark:
		"benchmarks/{sample}.gmap.benchmark.txt"
	shell:
		"module load gmap; gmap {params.ref_param} {params.rg} --format=samse --cross-species -K 350000 {input.fasta} | samtools view -Sb > {output}"

rule samtools_sort:
	input:
		"mapped_reads/{sample}.bam"
	output:
		"sorted_reads/{sample}.sorted.bam"
	log:
		"logs/samtools_sort/{sample}.log"
	benchmark:
		"benchmarks/{sample}.sort.benchmark.txt"
	threads: 4		
	shell:
		"module load samtools; samtools sort --threads {threads} -T sorted_reads/{wildcards.sample} -O bam {input} > {output}"

rule samtools_rmdup:
	input:
		"sorted_reads/{sample}.sorted.bam"
	output:
		"rmduped_reads/{sample}.sorted.rmdup.bam"
	log:
		"logs/samtools_rmdup/{sample}.log"
	benchmark:
		"benchmarks/{sample}.rmdup.benchmark.txt"
	shell:
		"module load samtools; samtools rmdup {input} {output}"

rule samtools_index_rmduped:
	input:
		"rmduped_reads/{sample}.sorted.rmdup.bam"
	output:
		"rmduped_reads/{sample}.sorted.rmdup.bam.bai"
	log:
		"logs/smatools_index_sorted{sample}.log"
	benchmark:
		"benchmarks/{sample}.index_sorted.benchmark.txt"
	shell:
		"module load samtools; samtools index {input}"
		

rule samtools_index_sorted:
	input:
		"sorted_reads/{sample}.sorted.rmdup.bam"
	output:
		"sorted_reads/{sample}.sorted.rmdup.bam.bai"
	log:
		"logs/samtools_index_sorted/{sample}.log"
	benchmark:
		"benchmarks/{sample}.index_sorted.benchmark.txt"	
	shell:
		"module load samtools; samtool index {input}"

rule freebayes:
        input:
                ref="/home/d.ence/projects/pinus_taeda_L/aligning_unigenes_w_gmap/CCLONES_all_seq.fasta",
                bam=expand("rmduped_reads/{sample}.sorted.rmdup.bam",sample=config["samples"])
                #bam=expand("indel_realigned/{sample}.sorted.rmdup.bam",sample=config["samples"])
                #pop_file=config["pop_file"]
        output:
                "calls/all_samples.ptaeda.wegrzyn_transcriptomes.freebayes_populations.vcf"
        threads: 4
        log:
                "logs/freebayes/all_samples.log"
        benchmark:
                "benchmarks/all_samples.freebayes.benchmark.txt"
        shell:
                #"module load freebayes; freebayes --populations {input.pop_file} -f {input.ref} {input.bam} > {output}"
                "module load freebayes; freebayes -f {input.ref} {input.bam} > {output}"

rule report:
	input:
		T1="calls/all_samples.ptaeda.wegrzyn_transcriptomes.freebayes_populations.vcf",
		T2=expand("benchmarks/{sample}.gmap.benchmark.txt",sample=config["samples"])
	output:
		"freebayes_calls.report.html"
	run:
		from snakemake.utils import report
		with open(input.T1) as vcf:
			n_calls = sum(1 for l in vcf if not l.startswith("#"))

		report("""
		A first test of a structural variant calling workflow
		=====================================================
		
		Reads were mapped to the CCLONES EST contigs 
		reference genome and variants were called jointly with
		freebayes.

		This resulted in {n_calls} variants (see Table T1_).
		Benchmark results for gmap can be found in the tables T2_.
		""", output[0], T1=input[0])
