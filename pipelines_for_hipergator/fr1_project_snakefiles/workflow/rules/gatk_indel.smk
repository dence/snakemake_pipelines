rule gatk_indel_creator:
	input:
		bam="results/rmduped/{sample}-{unit}.sorted.rmdup.bam",
		bai="results/rmduped/{sample}-{unit}.sorted.rmdup.bam.bai"
	output:
		temp("results/realigner_intervals/{sample}-{unit}.intervals")
	params:
		ref=get_reference
	log:
		"logs/realigner_intervals/{sample}-{unit}.intervals.log"
	benchmark:
		"benchmarks/{sample}-{unit}.intervals.benchmark.txt"
	shell:
		"module load gatk;  java -jar -Xmx9g /apps/gatk/3.7.0/GenomeAnalysisTK.jar -T RealignerTargetCreator --filter_reads_with_N_cigar -I {input.bam} -o {output} -R {params.ref}"

rule gatk_indel_realign:
	input:
		bam="results/rmduped/{sample}-{unit}.sorted.rmdup.bam",
		bai="results/rmduped/{sample}-{unit}.sorted.rmdup.bam.bai",
		interval="results/realigner_intervals/{sample}-{unit}.intervals"
	output:
		"results/realigned/{sample}-{unit}.sorted.rmdup.realigned.bam"
	params:
		ref=get_reference
	log:
		"logs/realigner/{sample}-{unit}.realigner.log"
	benchmark:
		"benchmarks/{sample}-{unit}.intervals.benchmark.txt"
	shell:
		"module load gatk; java -jar -Xmx4g /apps/gatk/3.7.0/GenomeAnalysisTK.jar -T IndelRealigner --filter_reads_with_N_cigar -R {params.ref} -I {input.bam} -targetIntervals {input.interval} -o {output} &> {log}"
