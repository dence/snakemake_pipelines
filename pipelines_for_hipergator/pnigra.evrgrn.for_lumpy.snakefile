#samples=["CH1","CH2","CH3","CH4","PA1","PA2","PA3","PA4"]


rule all:
	input:
		"lumpy.report.html"

rule bwa_map:
	input:
		#"/ufrc/kirst/d.ence/Pnigra_Sempervirens/ref_genomes/P.nigra/Consensus_71077-308_version0915-MER.fa",
		"/ufrc/kirst/d.ence/ref_genomes/populus/ref_genomes/Potr/v3/Ptrichocarpa_v3.0_210.fa",
		"/ufrc/kirst/share/Pnigra_Sempervirens/raw/{sample}/{sample}_R2_reverse_paired.trim.fq.gz",
		"/ufrc/kirst/share/Pnigra_Sempervirens/raw/{sample}/{sample}_R1_forward_paired.trim.fq.gz",
	output:
		protected("mapped_reads/{sample}.bam")
	params:
		rg="@RG\\tID:{sample}\\tSM:{sample}"
	log:
		"logs/bwa_mem/{sample}.log"
	benchmark:
		"benchmarks/{sample}.bwa.benchmark.txt"
	threads: 10
	shell:
		"module load bwa; module load samtools; (bwa mem -R '{params.rg}' -t {threads} {input} | samtools view -Sb - > {output}) 2> {log}"
rule extract_discordant:
	input:
		"mapped_reads/{sample}.bam"
	output:
		protected("discordants/{sample}.discordants.unsorted.bam")
	shell:
		"module load samtools; samtools view -b -F 1294 {input} > {output}"

rule extract_split_reads:
	input:
		"discordants/{sample}.discordants.unsorted.bam"
	output:
		protected("split_reads/{sample}.splitters.unsorted.bam")
	params:
		sample="{sample}"
	shell:
		"module load intel/2016.0.109; module load python/2.7.11; module load lumpy; module load samtools; "
		" samtools view -h {input} > discordants/{params.sample}.discordants.unsorted.sam ; $HPC_LUMPY_DIR/scripts/extractSplitReads_BwaMem --inFile discordants/{params.sample}.discordants.unsorted.sam | "
		" samtools view -Sb - > {output}"

rule samtools_sort_bams:
	input:
		"mapped_reads/{sample}.bam"
	output:
		"sorted_bams/{sample}.sorted.bam"
	log:
		"logs/sorted/{sample}.log"
	benchmark:
		"benchmarks/{sample}.sorted.benchmark.txt"
	shell:
		"module load samtools; samtools sort {input} > {output}"

rule samtools_rmdup:
	input:
		"sorted_bams/{sample}.sorted.bam"
	output:
		"rmduped_reads/{sample}.sorted.rmdup.bam"
	log:
		"logs/samtools_rmdup/{sample}.log"
	benchmark:
		"benchmarks/{sample}.rmdup.benchmark.txt"
	shell:
		"module load samtools; samtools rmdup {input} {output}"

rule samtools_index_rmdup:
	input:
		"rmduped_reads/{sample}.sorted.rmdup.bam"
	output:
		"rmduped_reads/{sample}.sorted.rmdup.bam.bai"
	shell:
		"module load samtools; samtools index {input}"

rule gatk_indel_creator:
	input:
		bam="rmduped_reads/{sample}.sorted.rmdup.bam",
		bai="rmduped_reads/{sample}.sorted.rmdup.bam.bai"
	output:
		"realigner_intervals/{sample}.intervals"
	params:
		 ref=config["ref"]["V3"]["full"]
	log:
		"logs/realigner_intervals/{sample}.intervals.log"
	benchmark:
		"benchmarks/{sample}.intervals.benchmark.txt"
	shell:
		"module load gatk;  java -jar -Xmx4g /apps/gatk/3.7.0/GenomeAnalysisTK.jar -T RealignerTargetCreator -I {input.bam} -o {output} -R {params.ref} &> {log}"

rule gatk_indel_realign:
	input:
		bam="rmduped_reads/{sample}.sorted.rmdup.bam",
		interval="realigner_intervals/{sample}.intervals"
	output:
		"realigned_bams/{sample}.sorted.rmdup.realigned.bam"
	params:
		ref=config["ref"]["V3"]["full"]
	log:
		"logs/realigner/{sample}.realigner.log"
	benchmark:
		"benchmarks/{sample}.intervals.benchmark.txt"
	shell:
		"module load gatk; java -jar -Xmx4g /apps/gatk/3.7.0/GenomeAnalysisTK.jar -T IndelRealigner -R {params.ref} -I {input.bam} -targetIntervals {input.interval} -o {output} &> {log}"

rule samtools_index_realigned:
	input:
		"realigned_bams/{sample}.sorted.rmdup.realigned.bam"
	output:
		"realigned_bams/{sample}.sorted.rmdup.realigned.bam.bai"
	log:
		"logs/samtools_index_realigned/{sample}.log"
	benchmark:
		"benchmarks/{sample}.index_realigned.benchmark.txt"
	shell:
		"module load samtools; samtools index {input}"
rule freebayes:
	input:
		ref=config["ref"]["V3"]["full"],
		bam=expand("realigned_bams/{sample}.sorted.rmdup.realigned.bam",sample=config["samples"]),
		bai=expand("realigned_bams/{sample}.sorted.rmdup.realigned.bam.bai",sample=config["samples"])
	output:
		"calls/all_samples.bwa.freebayes.vcf"
	log:
		"logs/freebayes/all_samples.log"
	benchmark:
		"benchmark/all_samples.freebayes.benchmark.txt"
	shell:
		"module load freebayes; freebayes -f {input.ref} {input.bam} > {output} 2> {log}"

rule samtools_sort_discordants:
	input:
		"discordants/{sample}.discordants.unsorted.bam"
	output:
		protected("discordants/{sample}.discordants.bam")
	shell:
		"module load samtools; samtools sort {input} > {output}"

rule samtools_sort_splitters:
	input:
		"split_reads/{sample}.splitters.unsorted.bam"
	output:
		protected("split_reads/{sample}.splitters.bam")
	shell:
		"module load samtools; samtools sort {input} > {output}"

rule insert_size_stats:
	input:
		"mapped_reads/{sample}.bam"
	params:
		rg="{sample}"
	output:
		histo=protected("histos/{sample}.histo"),
		output=protected("histos/{sample}.script_out")
	shell:
		"module load samtools; module load lumpy; samtools view -r {params.rg} {input} | "
		"tail -n+100000 | "
		" /apps/perl/perls/perl-5.16.0/bin/perl -w pairend_distro.pl -r 101 -X 4 -N 10000 -o {output.histo} > {output.output}"

#rule run_lumpy:
#	input:
#		tmp_sample="{sample}",
#		histo="histos/{sample}.histo",
#		sorted="sorted_reads/{sample}.sorted.bam",
#		discordants="discordants/{sample}.discordants.bam",
#		splitters="split_reads/{sample}.splitters.bam"
#	output:
#		"calls/multi_sample.lumpy.vcf"
#	run:
#		command_line = "lumpy -mw 4 -tt 0 "	
#		for tmp in ["CH1","CH2","CH3","CH4","PA1","PA2","PA3","PA4"]:
#			tmp_filename = "histos/" + tmp + ".script_out"
#			tmp_file = open(tmp_filename,'r')
#			mean=0
#			stdev=0
#			for line in tmp_file:
#				parts=line.split("\s+")
#				mean = parts[0].split(":")[1]
#				stdev = parts[1].split(":")[1]
#			close tmp_file
#
#			line = "-pe id:" + tmp + "bam_file:mapped_reads/" + tmp + ".bam,histo_file:histos/" + tmp + ".histo.histo,mean:" + mean + ",stdev:" + stdev + ",read_length:101,min_non_overlap:101,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20 "
#			line = line + " -sr id:" + tmp + ",bam_file:mapped_reads/" + tmp + ".bam,back_distance:10,weight:1,min_mapping_threshold:20 "
#			command_line = command_line + line
#		command_line=command_line + " > {output}"
#		print(command_line)

rule report:
	input:
		T1="calls/all_samples.bwa.freebayes.vcf",
		T2=expand("benchmarks/{sample}.bwa.benchmark.txt",sample=config["samples"])
	output:
		"lumpy.report.html"
	run:
		from snakemake.utils import report
		with open(input.T1) as vcf:
			n_calls = sum(1 for l in vcf if not l.startswith("#"))
		
		report("""
		A first test of a structural variant calling workflow
		=====================================================
		
		Reads were mapped to the P.nigra
		reference genome and variants were called jointly with
		gatk unified genotyper.

		This resulted in {n_calls} variants (see Table T1_).
		Benchmark results for BWA can be found in the tables T2_.
		""", output[0], T1=input[0])
