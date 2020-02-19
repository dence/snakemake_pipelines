rule all:
	input:
		"freebayes_calls.report.html"
rule prepare_superreads:
	input:
		#fastqs to assemble
                fq1="/ufrc/kirst/d.ence/pinus_taeda_L/wegrzyn_elite_transcriptome_raw_data/test.{sample}.R1.fastq",
                fq2="/ufrc/kirst/d.ence/pinus_taeda_L/wegrzyn_elite_transcriptome_raw_data/test.{sample}.R1.fastq",
#		#template config file
		template=config["conf_template"],
		script_base=config["script_base"]
	output:
#		#prepared config file (prepared from template)
		"superreads/sr_config.{sample}.txt"
	shell:
#		#make a directory for the superreads run
#		#copy/customize the config file
		#"cat {input.template} | perl -ane 's/FQ1/{input.fq1}/; s/FQ2/{input.fq2}/; print $_' > {output} "
		"module load python/2.7.10 ; python {input.script_base}/make_sr_config_file.py {input.fq1} {input.fq2} {input.template}  > {output}"

rule assemble_superreads:
	input:
                fq1="/ufrc/kirst/d.ence/pinus_taeda_L/wegrzyn_elite_transcriptome_raw_data/test.{sample}.R1.fastq",
                fq2="/ufrc/kirst/d.ence/pinus_taeda_L/wegrzyn_elite_transcriptome_raw_data/test.{sample}.R2.fastq",
		sr_conf="superreads/sr_config.{sample}.txt" #prepeared config file
	output:
		long_reads="superreads/{sample}.LongReads.fq",
		tmp_dir="tmp_superreads_{sample}"
	log:
		"logs/superreads/{sample}.log"
	benchmark:
		"benchmarks/{sample}.superreads.benchmark.txt"
	shell:
		"module load masurca; mkdir {output.tmp_dir} ; cd {output.tmp_dir} ; superreads.pl {input.fq1} {input.fq2} -c ../{input.sr_conf} -l ../{output.long_reads}"
#
rule tophat_supereads_align:
	input:
		fq1="/ufrc/kirst/d.ence/pinus_taeda_L/wegrzyn_elite_transcriptome_raw_data/test.{sample}.R1.fastq",
		long_reads="superreads/{sample}.LongReads.fq",
                fq2="/ufrc/kirst/d.ence/pinus_taeda_L/wegrzyn_elite_transcriptome_raw_data/test.{sample}.R2.fastq"
	output:
		"tmp_tophat_{sample}/accepted_hits.bam"
	params:
		rg='-p 20 --rg-id={sample} --rg-sample={sample} --rg-library={sample} --rg-platform=illumina',
		ref_dir=config["reference"]["V1_01"]["dir"],
		ref_full=config["reference"]["V1_01"]["full"],
		output_dir='tmp_tophat_{sample}'
	log:
		"logs/tophat/{sample}.log"
	benchmark:
		"benchmarks/{sample}.tophat.benchmark.txt"
	shell:
		"module load tophat ; export BOWTIE2_INDEXES=""; tophat {params.rg} -o {params.output_dir} {params.ref_full} {input.fq1},{input.long_reads} {input.fq2}"

#rule tophat_raw_align:
#	input:
#		#"superreads/{sample}.LongReads.fq",	
#		"/ufrc/kirst/d.ence/pinus_taeda_L/wegrzyn_elite_transcriptome_raw_data/test.{sample}.R1.fastq",
#		"/ufrc/kirst/d.ence/pinus_taeda_L/wegrzyn_elite_transcriptome_raw_data/test.{sample}.R2.fastq",
#		
#	output#:
#		"mapped_reads/{sample}.bam"
#	params:
#		rg='-p 20 --rg-id={sample} --rg-sample={sample} --rg-library={sample} --rg-platform=illumina',
#		ref_dir=config["reference"]["V1_01"]["dir"],
#		ref_full=config["reference"]["V1_01"]["full"],
#		output_dir='tophat_{sample}'
#	log:
#		"logs/tophat/{sample}.log"
#	benchmark:
#		"benchmarks/{sample}.tophat.benchmark.txt"
#	shell:
#		"module load tophat;  export BOWTIE2_INDEXES=""; tophat {params.rg} -o {params.output_dir} {params.ref_full} {input} ; samtools view -b {params.output_dir}/accepted_hits.bam > {output} "

rule samtools_sort:
	input:
		"tmp_tophat_{sample}/accepted_hits.bam"
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
rule stringtie:
	input:
		"rmduped_reads/{sample}.sorted.rmdup.bam"
	output:
		"stringtie_gtfs/{sample}.stringtie.out.gtf"
	log:
		"logs/stringtie/{sample}.log"
	benchmark:
		"benchmarks/{sample}.stringtie.benchmark.txt"
	shell:
		"module load stringtie ; stringtie -o {output} {input}"


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
                ref=config["reference"]["V1_01"]["full"],
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
		#T1="calls/all_samples.ptaeda.wegrzyn_transcriptomes.freebayes_populations.vcf",
		T2=expand("stringtie_gtfs/{sample}.stringtie.out.gtf",sample=config["samples"])
	output:
		"freebayes_calls.report.html"
	run:
		from snakemake.utils import report
		with open(input.T2) as gtf:
			n_calls = sum(1 for l in gtf if not l.startswith("#"))

		report("""
		A first test of a reference-based transcript assembly workflow
		=====================================================
		
		Reads were mapped to the Ptaeda v1.1 reference assembly 
		reference genome and variants were called jointly with
		freebayes.

		This resulted in {n_calls} variants (see Table T2_).
		""", output[0], T2=input[0])
