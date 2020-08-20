
#Daniel Ence, 08/20/2020


#Need to get a

rule trim_reads:
	input:
        get_fastqs
		fq1="/orange/kirst/d.ence/pinus_taeda_L/fr1_sequencing/RAPiD-Genomics_HJYNGBBXX_UFL_104301_{sample}_R1_001.fastq.gz",
		fq2="/orange/kirst/d.ence/pinus_taeda_L/fr1_sequencing/RAPiD-Genomics_HJYNGBBXX_UFL_104301_{sample}_R2_001.fastq.gz",
		prev_log="logs/fastqc_pre_out/{sample}.log"
	output:
		fq1=temp("trimmed_reads/{sample}/{sample}.trimmed.R1.fastq.gz"),
		fq2=temp("trimmed_reads/{sample}/{sample}.trimmed.R2.fastq.gz")
	log:
		"logs/cutadapt/{sample}.log"
	benchmark:
		"benchmarks/{sample}.fastqc.benchmark.txt"
	shell:
		"module load python/2.7.14; module load gcc/5.2.0; module load cutadapt/1.18; cutadapt -u 10 -U 10 -q30,30 --minimum-length=50 -o {output.fq1} -p {output.fq2} {input.fq1} {input.fq2} &> {log}"
