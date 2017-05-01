rule all:
	input:
		"lumpy.report.html"

rule bwa_map:
	input:
		"/ufrc/kirst/d.ence/Pnigra_Sempervirens/ref_genomes/P.nigra/Consensus_71077-308_version0915-MER.fa",
		"/ufrc/kirst/share/Pnigra_Sempervirens/raw/{sample}/{sample}_R2_reverse_paired.trim.fq.gz",
		"/ufrc/kirst/share/Pnigra_Sempervirens/raw/{sample}/{sample}_R1_forward_paired.trim.fq.gz",
	output:
		"mapped_reads/{sample}.bam"
	params:
		rg="@RG\tID:{sample}\tSM:{sample}"
	log:
		"logs/bwa_mem/{sample}.log"
	benchmark:
		"benchmarks/{sample}.bwa.benchmark.txt"
	threads: 4
	shell:
		"module load bwa; bwa mem -R '{params.rg}' -t {threads} {input} | samtools view -Sb > {output}"
rule extract_discordant:
	input:
		"mapped_reads/{sample}.bam"
	output:
		"discordants/{sample}.discordants.unsorted.bam"
	shell:
		"module load samtools; samtools view -b -F 1294 {input} > {outpu}"

rule extract_split_reads:
	input:
		"discordants/{sample}.discordants.unsorted.bam"
	output:
		"split_reads/{sample}.splitters.unsorted.bam"
	shell:
		"module load intel/2016.0.109; module load python/2.7.11; module load lumpy; module load samtools; | "
		" $HPC_LUMPY_DIR/scripts/extractSplitReads_BwaMem --inFile {input} | "
		" samtools view -Sb - > {output}"

rule samtools_sort_splitters:
	input:
		"split_reads/{sample}.splitters.unsorted.bam"
	output:
		"split_reads/{sample}.splitters.bam"
	shell:
		"module load samtools samtools sort {input} > {output}"

rule insert_size_stats:
	input:
		"mapped_reads/{sample}.bam"
	output:
		"histos/{sample}.histo"
	shell:
		"module laod samtools; module load lumpy; samtools view -r {params.rg} {input} | "
		"tail -n+100000 | "
		" $HPC_LUMPY_DIR/scripts/pairend_distro.pl -r 101 -X 4 -N 10000 -o {output} > params.script_out"

#rule run_lumpy:
#	input:
#		tmp_sample="{sample}",
#		histo="histos/{sample}.histo",
#		discordants="discordants/{sample}.discordants.bam",
#		splitters="split_reads/{sample}.splitters.bam"
#	output:
#		"calls/multi_sample.lumpy.vcf"
#	shell:
#		command_line = "lumpy -mw 4 -tt 0 " 
#		print(config["samples"])
#		#for tmp in config["samples"]:
#			line = "-pe id:{input.tmp_sample},bam_file:{input.discordants},histo_file:{input.histo},mean:500,stdev:50,read_length:101,min_non_overlap:101,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20 "
#		#	line = line + " -sr id:{tmp_sample},bam_file:{input.splitters},back_distance:10,weight:1,min_mapping_threshold:20
#		#	command_line = command_line + line
#		#command_line = command_line + " > {output}"
#		#command_line
#

rule report:
	input:
		T1=expand("histos/{sample}.histo", sample=config["samples"]),
		T2=expand("split_reads/{sample}.splitters.bam", sample=config["samples"]), 
		T3=expand("discordants/{sample}.discordants.unsorted.bam", sample=config["samples']),
		T4=expand("benchmarks/{sample}.bwa.benchmark.txt",sample=config["samples"])
	output:
		"lumpy.report.html"
	run:
		report("""
		A first test of a structural variant calling workflow
		=====================================================
		
		Reads were mapped to the P.nigra
		reference genome and variants were called jointly with
		gatk unified genotyper.

		This resulted in X variants (see Table T1_).
		Benchmark results for BWA can be found in the tables T2_.
		""", output[0], T1=input[0])
