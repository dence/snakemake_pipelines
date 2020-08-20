rule all:
	input:
		"freebayes_calls.report.html",
		"mosdepth.report.html",
		"touched_successfully.txt"

rule pretrim_fastqc:
	input:
		fq1="/orange/kirst/d.ence/pinus_taeda_L/fr1_sequencing/RAPiD-Genomics_HJYNGBBXX_UFL_104301_{sample}_R1_001.fastq.gz",
		fq2="/orange/kirst/d.ence/pinus_taeda_L/fr1_sequencing/RAPiD-Genomics_HJYNGBBXX_UFL_104301_{sample}_R2_001.fastq.gz"
	output:
		out="fastqc_out/pre_process/{sample}/",
		log_file="logs/fastqc_pre_out/{sample}.log"
	log:
		"logs/fastqc_pre_out/{sample}.log"
	benchmark:
		"benchmarks/{sample}.fastqc.benchmark.txt"
	shell:
		"module load fastqc ; fastqc -f fastq -o {output.out} {input} &> {log}"
rule trim_reads:
	input:
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
rule posttrim_fastqc:
	input:
		"trimmed_reads/{sample}/{sample}.trimmed.R1.fastq.gz",
		"trimmed_reads/{sample}/{sample}.trimmed.R2.fastq.gz"
	output:
		"fastqc_out/post_trim/{sample}/",
		"logs/fastqc_post_out/{sample}.log"
	log:
		"logs/fastqc_post_out/{sample}.log"
	benchmark:
		"benchmarks/{sample}.fastqc.benchmark.txt"
	shell:
		"module load fastqc ; fastqc -f fastq -o {output} {input} &> {log}"

rule hisat2_align:
	input:
		fq1="trimmed_reads/{sample}/{sample}.trimmed.R1.fastq.gz",
		fq2="trimmed_reads/{sample}/{sample}.trimmed.R2.fastq.gz"
	output:
		temp("mapped_reads/{sample}.hisat2.bam")
	log:
		"logs/hisat2/{sample}.log"
	params:
		rg="--rg-id={sample} --rg \"SM:{sample}\" --rg \"LB:{sample}\" --rg \"PL:ILLUMINA\"",
		ref=config["reference"]["V1_01"]["custom"]
		#ref=config["reference"]["V2_01"]["custom"]
	benchmark:
		"benchmarks/{sample}.hisat2.benchmark.txt"
	shell:
		"module load hisat2; module load samtools; hisat2 -p 4 --mm --time -x {params.ref} {params.rg} -1 {input.fq1} -2 {input.fq2} | samtools view -bS - > {output} 2> {log}"

rule bwa_mem_align:
	input:
		fq1="trimmed_reads/{sample}/{sample}.trimmed.R1.fastq.gz",
		fq2="trimmed_reads/{sample}/{sample}.trimmed.R2.fastq.gz"
	params:
		sample="{sample}",
		ref=config["reference"]["V1_01"]["custom"]
	output:
		temp("mapped_reads/{sample}.bwa_mem.bam")
	log:
		"logs/bwa_mem/{sample}.bwa_mem.log"
	shell:
		"module load bwa/0.7.17 ; bwa mem {params.ref} " 
		+ "-R '@RG\\tID:{params.sample}\\tSM:{params.sample}\\tPL:ILLUMINA'" 
		+ "{input.fq1} {input.fq2} > {output} 2> {log}"	

rule samtools_sort_bwa_mem:
	input:
		"mapped_reads/{sample}.bwa_mem.bam"
	output:
		temp("sorted_reads/{sample}.bwa_mem.sorted.bam")
	log:
		"logs/samtools_sort/{sample}.log"
	benchmark:
		"benchmarks/{sample}.sort.benchmark.txt"
	threads: 4		
	shell:
		"module load samtools; samtools sort --threads {threads} -T sorted_reads/{wildcards.sample} -O bam {input} > {output}"

rule samtools_sort_hisat2:
	input:
		"mapped_reads/{sample}.hisat2.bam"
	output:
		temp("sorted_reads/{sample}.hisat2.sorted.bam")
	log:
		"logs/samtools_sort/{sample}.hisat2.log"
	benchmark:
		"benchmarks/{sample}.sort.benchmark.txt"
	threads: 4		
	shell:
		"module load samtools; samtools sort --threads {threads} -T sorted_reads/{wildcards.sample} -O bam {input} > {output} 2> {log}"

rule samtools_rmdup_bwa_mem:
	input:
		"sorted_reads/{sample}.bwa_mem.sorted.bam"
	output:
		temp("rmduped_reads/{sample}.bwa_mem.sorted.rmdup.bam")
	log:
		"logs/samtools_rmdup/{sample}.log"
	benchmark:
		"benchmarks/{sample}.rmdup.benchmark.txt"
	shell:
		"module load samtools; samtools rmdup {input} {output}"

rule samtools_rmdup_hisat2:
	input:
		"sorted_reads/{sample}.hisat2.sorted.bam"
	output:
		temp("rmduped_reads/{sample}.hisat2.sorted.rmdup.bam")
	log:
		"logs/samtools_rmdup/{sample}.log"
	benchmark:
		"benchmarks/{sample}.rmdup.benchmark.txt"
	shell:
		"module load samtools; samtools rmdup {input} {output}"

rule samtools_index_rmduped_bwa_mem:
	input:
		"rmduped_reads/{sample}.bwa_mem.sorted.rmdup.bam"
	output:
		temp("rmduped_reads/{sample}.bwa_mem.sorted.rmdup.bam.bai")
	log:
		"logs/samtools_index_sorted.{sample}.log"
	benchmark:
		"benchmarks/{sample}.index_sorted.benchmark.txt"
	shell:
		"module load samtools; samtools index {input}"

rule samtools_index_rmduped_hisat2:
	input:
		"rmduped_reads/{sample}.hisat2.sorted.rmdup.bam"
	output:
		temp("rmduped_reads/{sample}.hisat2.sorted.rmdup.bam.bai")
	log:
		"logs/samtools_index_sorted.{sample}.log"
	benchmark:
		"benchmarks/{sample}.index_sorted.benchmark.txt"
	shell:
		"module load samtools; samtools index {input}"

rule samtools_index_sorted_bwa_mem:
	input:
		"sorted_reads/{sample}.bwa_mem.sorted.rmdup.bam"
	output:
		temp("sorted_reads/{sample}.bwa_mem.sorted.rmdup.bam.bai")
	log:
		"logs/samtools_index_sorted/{sample}.log"
	benchmark:
		"benchmarks/{sample}.index_sorted.benchmark.txt"	
	shell:
		"module load samtools; samtool index {input}"

rule samtools_index_sorted_hisat2:
	input:
		"sorted_reads/{sample}.hisat2.sorted.rmdup.bam"
	output:
		temp("sorted_reads/{sample}.hisat2.sorted.rmdup.bam.bai")
	log:
		"logs/samtools_index_sorted/{sample}.log"
	benchmark:
		"benchmarks/{sample}.index_sorted.benchmark.txt"	
	shell:
		"module load samtools; samtool index {input}"

rule gatk_indel_creator_bwa_mem:
	input:
		bam="rmduped_reads/{sample}.bwa_mem.sorted.rmdup.bam",
		bai="rmduped_reads/{sample}.bwa_mem.sorted.rmdup.bam.bai"
	output:
		temp("realigner_intervals/{sample}.bwa_mem.intervals")
	params:
		ref=config["reference"]["V1_01"]["custom"]
	log:
		"logs/realigner_intervals/{sample}.intervals.log"
	benchmark:
		"benchmarks/{sample}.ineervals.benchmark.txt"
	shell:
		"module load gatk;  java -jar -Xmx9g /apps/gatk/3.7.0/GenomeAnalysisTK.jar -T RealignerTargetCreator --filter_reads_with_N_cigar -I {input.bam} -o {output} -R {params.ref}"

rule gatk_indel_creator_hisat2:
	input:
		bam="rmduped_reads/{sample}.hisat2.sorted.rmdup.bam",
		bai="rmduped_reads/{sample}.hisat2.sorted.rmdup.bam.bai"
	output:
		temp("realigner_intervals/{sample}.hisat2.intervals")
	params:
		ref=config["reference"]["V1_01"]["custom"]
	log:
		"logs/realigner_intervals/{sample}.intervals.hisat2.log"
	benchmark:
		"benchmarks/{sample}.intervals.benchmark.txt"
	shell:
		"module load gatk;  java -jar -Xmx10g /apps/gatk/3.7.0/GenomeAnalysisTK.jar -T RealignerTargetCreator --filter_reads_with_N_cigar -I {input.bam} -o {output} -R {params.ref} &> {log}"

rule gatk_indel_realign_bwa_mem:
	input:
		bam="rmduped_reads/{sample}.bwa_mem.sorted.rmdup.bam",
		bai="rmduped_reads/{sample}.bwa_mem.sorted.rmdup.bam.bai",
		interval="realigner_intervals/{sample}.bwa_mem.intervals"
	output:
		temp("realigned_bams/{sample}.bwa_mem.sorted.rmdup.realigned.bam")
	params:
		ref=config["reference"]["V1_01"]["custom"]
	log:
		"logs/realigner/{sample}.realigner.log"
	benchmark:
		"benchmarks/{sample}.intervals.benchmark.txt"
	shell:
		"module load gatk; java -jar -Xmx4g /apps/gatk/3.7.0/GenomeAnalysisTK.jar -T IndelRealigner --filter_reads_with_N_cigar -R {params.ref} -I {input.bam} -targetIntervals {input.interval} -o {output} &> {log}"

rule gatk_indel_realign_hisat2:
	input:
		bam="rmduped_reads/{sample}.hisat2.sorted.rmdup.bam",
		bai="rmduped_reads/{sample}.hisat2.sorted.rmdup.bam.bai",
		interval="realigner_intervals/{sample}.hisat2.intervals"
	output:
		temp("realigned_bams/{sample}.hisat2.sorted.rmdup.realigned.bam")
	params:
		ref=config["reference"]["V1_01"]["custom"]
	log:
		"logs/realigner/{sample}.realigner.hisat2.log"
	benchmark:
		"benchmarks/{sample}.intervals.benchmark.txt"
	shell:
		"module load gatk; java -jar -Xmx4g /apps/gatk/3.7.0/GenomeAnalysisTK.jar -T IndelRealigner --filter_reads_with_N_cigar -R {params.ref} -I {input.bam} -targetIntervals {input.interval} -o {output}"

rule samtools_index_realigned_bwa_mem:
	input:
		"realigned_bams/{sample}.bwa_mem.sorted.rmdup.realigned.bam"
	output:
		temp("realigned_bams/{sample}.bwa_mem.sorted.rmdup.realigned.bai")
	log:
		"logs/samtools_index_realigned/{sample}.log"
	benchmark:
		"benchmarks/{sample}.index_realigned.benchmark.txt"
	shell:
		"module load samtools; samtool index {input}"

rule samtools_index_realigned_hisat2:
	input:
		"realigned_bams/{sample}.hisat2.sorted.rmdup.realigned.bam"
	output:
		temp("realigned_bams/{sample}.hisat2.sorted.rmdup.realigned.bai")
	log:
		"logs/samtools_index_realigned/{sample}.log"
	benchmark:
		"benchmarks/{sample}.index_realigned.benchmark.txt"
	shell:
		"module load samtools; samtool index {input}"

rule merge_lane_bams_bwa_mem:
	input:
		#making a dumb assumption about the names of the bams to merged. specific to the Fr1 project. DE
		L4_bam="rmduped_reads/{sample}_L004.bwa_mem.sorted.rmdup.bam",
		L5_bam="rmduped_reads/{sample}_L005.bwa_mem.sorted.rmdup.bam"
	output:
		temp("merged_lane_bams/{sample}.bwa_mem.sorted.rmdup.merged.bam")
	log:
		"logs/picard_merge_sam_files/{sample}.log"
	shell:
		"module load picard; java -jar $HPC_PICARD_DIR/picard.jar MergeSamFiles I={input.L4_bam} I={input.L5_bam} O={output} &> {log}"

rule merge_lane_bams_hisat2:
	input:
		#making a dumb assumption about the names of the bams to merged. specific to the Fr1 project. DE
		L4_bam="rmduped_reads/{sample}_L004.hisat2.sorted.rmdup.bam",
		L5_bam="rmduped_reads/{sample}_L005.hisat2.sorted.rmdup.bam"
	output:
		temp("merged_lane_bams/{sample}.hisat2.sorted.rmdup.merged.bam")
	log:
		"logs/picard_merge_sam_files/{sample}.log"
	shell:
		"module load picard; java -jar $HPC_PICARD_DIR/picard.jar MergeSamFiles I={input.L4_bam} I={input.L5_bam} O={output} &> {log}"

rule index_merged_bams_bwa_mem:
	input:
		"merged_lane_bams/{sample}.bwa_mem.sorted.rmdup.merged.bam"
	output:
		temp("merged_lane_bams/{sample}.bwa_mem.sorted.rmdup.merged.bam.bai")
	shell:
		"module load samtools; samtools index {input}"

rule index_merged_bams_hisat2:
	input:
		"merged_lane_bams/{sample}.hisat2.sorted.rmdup.merged.bam"
	output:
		temp("merged_lane_bams/{sample}.hisat2.sorted.rmdup.merged.bam.bai")
	shell:
		"module load samtools; samtools index {input}"

rule samtools_index_replaced_bwa_mem:
	input:
		"RG_replaced_bams/{sample}.bwa_mem.sorted.rmdup.merged.bam"
	output:
		temp("RG_replaced_bams/{sample}.bwa_mem.sorted.rmdup.merged.bam.bai")
	log:
		"logs/samtools_index_realigned/{sample}.log"
	benchmark:
		"benchmarks/{sample}.index_realigned.benchmark.txt"
	shell:
		"module load samtools; samtools index {input}"

rule samtools_index_replaced_hisat2:
	input:
		"RG_replaced_bams/{sample}.hisat2.sorted.rmdup.merged.bam"
	output:
		temp("RG_replaced_bams/{sample}.hisat2.sorted.rmdup.merged.bam.bai")
	log:
		"logs/samtools_index_realigned/{sample}.log"
	benchmark:
		"benchmarks/{sample}.index_realigned.benchmark.txt"
	shell:
		"module load samtools; samtools index {input}"

rule ReplaceRG_merged_bwa_mem:
	input:
		"merged_lane_bams/{sample}.bwa_mem.sorted.rmdup.merged.bam"
	output:
		temp("RG_replaced_bams/{sample}.bwa_mem.sorted.rmdup.merged.bam")
	params:
		RG_fields="RGID={sample} RGLB={sample} RGPL=illumina RGPU={sample} RGSM={sample}"
	log:
		"logs/picard_replaceRG/{sample}.log"
	shell:
		"module load picard; java -jar $HPC_PICARD_DIR/picard.jar AddOrReplaceReadGroups I={input} O={output} {params.RG_fields} &> {log}"

rule ReplaceRG_merged_hisat2:
	input:
		"merged_lane_bams/{sample}.hisat2.sorted.rmdup.merged.bam"
	output:
		temp("RG_replaced_bams/{sample}.hisat2.sorted.rmdup.merged.bam")
	params:
		RG_fields="RGID={sample} RGLB={sample} RGPL=illumina RGPU={sample} RGSM={sample}"
	log:
		"logs/picard_replaceRG/{sample}.log"
	shell:
		"module load picard; java -jar $HPC_PICARD_DIR/picard.jar AddOrReplaceReadGroups I={input} O={output} {params.RG_fields} &> {log}"

rule samtools_index_merged_bwa_mem:
        input:
                "merged_lane_bams/{sample}.bwa_mem.sorted.rmdup.merged.bam"
        output:
                temp("merged_lane_bams/{sample}.bwa_mem.sorted.rmdup.merged.bam.bai")
        log:
                "logs/samtools_index_merged/{sample}.log"
        benchmark:
                "benchmarks/{sample}.index_merged.benchmark.txt"
        shell:
                "module load samtools; samtools index {input}"

rule samtools_index_merged_hisat2:
        input:
                "merged_lane_bams/{sample}.hisat2.sorted.rmdup.merged.bam"
        output:
                temp("merged_lane_bams/{sample}.hisat2.sorted.rmdup.merged.bam.bai")
        log:
                "logs/samtools_index_merged/{sample}.log"
        benchmark:
                "benchmarks/{sample}.index_merged.benchmark.txt"
        shell:
                "module load samtools; samtools index {input}"

rule mosdepth_bwa_mem:
	input:
		bam="RG_replaced_bams/{sample}.bwa_mem.sorted.rmdup.merged.bam",
		bai="RG_replaced_bams/{sample}.bwa_mem.sorted.rmdup.merged.bam.bai"
	params:
		options=" --by /home/d.ence/projects/pinus_taeda_L/Fr1_project/probe_files/RG_0610_merged.fixed.bed --no-per-base ./mosdepth/bwa_mem/{sample}",
		gz="mosdepth/bwa_mem/{sample}.regions.bed.gz"	
	output:
		bed=temp("mosdepth/bwa_mem/{sample}.regions.bed")	
	log:
		"logs/mosdepth/bwa_mem/{sample}.mosdepth.log"
	shell:
		"module load mosdepth; mosdepth {params.options} {input.bam} ; gunzip {params.gz}"	

rule mosdepth_hisat2:
	input:
		bam="RG_replaced_bams/{sample}.hisat2.sorted.rmdup.merged.bam",
		bai="RG_replaced_bams/{sample}.hisat2.sorted.rmdup.merged.bam.bai"
	params:
		options=" --by /home/d.ence/projects/pinus_taeda_L/Fr1_project/probe_files/RG_0610_merged.fixed.bed --no-per-base ./mosdepth/hisat2/{sample}",
		gz="mosdepth/hisat2/{sample}.regions.bed.gz"	
	output:
		bed=temp("mosdepth/hisat2/{sample}.regions.bed")	
	log:
		"logs/mosdepth/hisat2/{sample}.mosdepth.log"
	shell:
		"module load mosdepth; mosdepth {params.options} {input.bam} ; gunzip {params.gz}"	
		
rule compile_non_10_5_indv_mosdepth_report_bwa_mem:
	input:
		non_10_5_indv_mosdepths=expand("mosdepth/bwa_mem/{sample}.regions.bed",sample=config["non_10_5_indv"]),
		non_10_5_indv_list="/home/d.ence/projects/pinus_taeda_L/Fr1_project/aligning_Fr1_samples/aligning_to_custom_R_gene_V1_1_ref/calls/bam_list.non_10_5_indv.txt"
	output:
		filename="mosdepth/bwa_mem/non_10_5_indv.mosdepth.report.txt"
	log:
		"logs/mosdepth/bwa_mem/non_10_5_indv.mosdepth.report.log"
	run:
		open_list = open(input.non_10_5_indv_list,'r')
		header_line = "sample_ID\tfamily\t"
		header_line_finished=0
		sample_lines = []

		for l in open_list:
			short_sample_ID = l.split('/')[10].split('_i')[0]
			long_sample_ID = l.split('/')[10].split('.')[0]
			
			curr_mosdepth = open("mosdepth/bwa_mem/" + long_sample_ID + ".regions.bed")
			tmp = curr_mosdepth.readline() #skip the header line
			
			family = ""
			if(short_sample_ID == "P02_WG07"):
				family = "20-1010"
			else:
				family = "10-5_X_4-6664"
			new_line = short_sample_ID + "\t" + family + "\t"
			for depth_line in curr_mosdepth:
				region = depth_line.strip().split("\t")[0]
				print("region is:\t" + region) 
				mean_covg = depth_line.strip().split("\t")[3]
				print("mean_covg:\t" + mean_covg)
		
				new_line = new_line + mean_covg + "\t" 
				if(header_line_finished == 0):
					header_line = header_line + region + "_mean_covg\t"
			
			header_line_finished=1
			sample_lines.append(new_line)
			curr_mosdepth.close()
		
		all_sample = open(output.filename,'w')
		all_sample.write(header_line + "\n")
		for line in sample_lines:
			all_sample.write(line + "\n")
		all_sample.close()

rule compile_non_10_5_pooled_mosdepth_report_bwa_mem:
	input:
		non_10_5_pooled_mosdepths=expand("mosdepth/bwa_mem/{sample}.regions.bed",sample=config["non_10_5_pooled"]),
		non_10_5_pooled_list="/home/d.ence/projects/pinus_taeda_L/Fr1_project/aligning_Fr1_samples/aligning_to_custom_R_gene_V1_1_ref/calls/bam_list.non_10_5_pooled.txt"
	output:
		filename="mosdepth/bwa_mem/non_10_5_pooled.mosdepth.report.txt"
	log:
		filename="logs/mosdepth/bwa_mem/non_10_5_pooled.report.log"
	run:
		open_list = open(input.non_10_5_pooled_list,'r')
		header_line = "sample_ID\tfamily\t"
		header_line_finished=0
		sample_lines = []

		open_log = open(log.filename, 'w')
		for l in open_list:
			short_sample_ID = l.split('/')[10].split('_i')[0]
			long_sample_ID = l.split('/')[10].split('.')[0]

			curr_mosdepth = open("mosdepth/bwa_mem//" + long_sample_ID + ".regions.bed")
			open_log.write("opened this file:\tmosdepth/" + long_sample_ID + ".regions.bed\n")
			tmp = curr_mosdepth.readline()
			open_log.write(tmp + "\n")
		
			family=""
			open_log.write("short_sample_ID is:\t" + short_sample_ID + "\n")
			if(short_sample_ID == "P04_WA08"):
				family="CL04"
			elif(short_sample_ID == "P04_WB08"):
				family="CL05"
			elif(short_sample_ID == "P04_WC08"):
				family="PD18"
			elif(short_sample_ID == "P04_WD08"):
				family="PD35"
			new_line = short_sample_ID + "\t" + family + "\t"
			open_log.write("family is:\t" + family + "\n")
			
			for depth_line in curr_mosdepth:
				open_log.write("depth_line is:\t" + depth_line + "\n")
				region = depth_line.strip().split("\t")[0]
				mean_covg = depth_line.strip().split("\t")[3]
				
				open_log.write("region is:\t" + region + "\n")
				open_log.write("mean_covg is:\t" + mean_covg + "\n")

				if(header_line_finished == 0):
					header_line = header_line + region + "_mean_covg\t"
				new_line = new_line + mean_covg + "\t"
			sample_lines.append(new_line)	
			curr_mosdepth.close()
		
			header_line_finished=1
		all_sample = open(output.filename,'w')
		all_sample.write(header_line + "\n")
		for line in sample_lines:
			all_sample.write(line + "\n")
		all_sample.close()
		open_log.close() 

rule compile_10_5_megs_report_bwa_mem:
	input:
		megs_mosdepths=expand("mosdepth/bwa_mem/{sample}.regions.bed",sample=config["Fr1_meg_samples"]),
		megs_list="/home/d.ence/projects/pinus_taeda_L/Fr1_project/aligning_Fr1_samples/aligning_to_custom_R_gene_V1_1_ref/calls/bam_list.megs.txt"
	output:
		filename="mosdepth/bwa_mem/10_5_megs.mosdepth.report.txt"
	log:
		filename="logs/mosdepth/bwa_mem/10_5_megs.mosdepth.report.log"
	run:
		open_list = open(input.megs_list,'r')
		open_log = open(log.filename,'w')
		header_line = "sample_ID\tfamily\t"
		header_line_finished=0
		sample_lines = []

		for l in open_list:
			short_sample_ID = l.split('/')[10].split('_i')[0]
			long_sample_ID = l.split('/')[10].split('.')[0]

			curr_filename = "mosdepth/bwa_mem/" + long_sample_ID + ".regions.bed"
			curr_mosdepth = open(curr_filename,'r')
			tmp = curr_mosdepth.readline()
			
			new_line = short_sample_ID + "\t10_5_meg\t"
	
			line_count = 0	
			for depth_line in curr_mosdepth:
				region = depth_line.strip().split("\t")[0]
				mean_covg = depth_line.strip().split("\t")[3]
				
				if(header_line_finished == 0):
					header_line = header_line + region + "_mean_covg\t"
				new_line = new_line + mean_covg + "\t"
				line_count = line_count + 1
			sample_lines.append(new_line)
			curr_mosdepth.close()

			open_log.write("this file:\t" + curr_filename + " had this many lines:\t" + str(line_count) + "\n")
			header_line_finished = 1
		all_sample = open(output.filename,'w')
		all_sample.write(header_line + "\n")
		for line in sample_lines:
			all_sample.write(line + "\n")
		all_sample.close()
		open_log.close()

rule compile_10_5_prog_report_bwa_mem:
	input:
		prog_mosdepths=expand("mosdepth/bwa_mem/{sample}.regions.bed",sample=config["Fr1_prog"]),
		prog_list="/home/d.ence/projects/pinus_taeda_L/Fr1_project/aligning_Fr1_samples/aligning_to_custom_R_gene_V1_1_ref/calls/bam_list.Fr1_prog.txt"
	output:
		filename="mosdepth/bwa_mem/10_5_prog.mosdepth.report.txt"
	log:
		filename="logs/mosdepth/bwa_mem/10_5_prog.mosdepth.report.log"
	run:
		open_list = open(input.prog_list,'r')
		header_line = "sample_ID\tfamily\t"
		header_line_finished = 0
		sample_lines = []
		open_log = open(log.filename,'w')
	
		for l in open_list:
			short_sample_ID = l.split('/')[10].split('_i')[0]
			long_sample_ID = l.split('/')[10].split('.')[0]

			curr_filename = "mosdepth/bwa_mem/" + long_sample_ID + ".regions.bed"
			curr_mosdepth = open(curr_filename,'r')
			tmp = curr_mosdepth.readline()
			
			new_line = short_sample_ID + "\t10_5_prog\t"
		
			line_counter = 0	
			for depth_line in curr_mosdepth:
				region = depth_line.strip().split("\t")[0]
				mean_covg = depth_line.strip().split("\t")[3]
		
				if(header_line_finished == 0):	
					header_line = header_line + region + "_mean_covg\t"
				new_line = new_line + mean_covg + "\t"
				line_counter = line_counter + 1
			sample_lines.append(new_line)
			curr_mosdepth.close()
			open_log.write("this file:\t" + curr_filename + " had this many lines:\t" + str(line_counter) + "\n")
	
			header_line_finished = 1	
		all_sample = open(output.filename,'w')
		all_sample.write(header_line + "\n")
		for line in sample_lines:
			all_sample.write(line + "\n")
		all_sample.close()
		open_log.close()
		 
rule compile_non_10_5_indv_mosdepth_report_hisat2:
	input:
		non_10_5_indv_mosdepths=expand("mosdepth/hisat2/{sample}.regions.bed",sample=config["non_10_5_indv"]),
		non_10_5_indv_list="/home/d.ence/projects/pinus_taeda_L/Fr1_project/aligning_Fr1_samples/aligning_to_custom_R_gene_V1_1_ref/calls/bam_list.non_10_5_indv.txt"
	output:
		filename="mosdepth/hisat2/non_10_5_indv.mosdepth.report.txt"
	log:
		"logs/mosdepth/hisat2/non_10_5_indv.mosdepth.report.log"
	run:
		open_list = open(input.non_10_5_indv_list,'r')
		header_line = "sample_ID\tfamily\t"
		header_line_finished=0
		sample_lines = []

		for l in open_list:
			short_sample_ID = l.split('/')[10].split('_i')[0]
			long_sample_ID = l.split('/')[10].split('.')[0]
			
			curr_mosdepth = open("mosdepth/hisat2/" + long_sample_ID + ".regions.bed")
			tmp = curr_mosdepth.readline() #skip the header line
			
			family = ""
			if(short_sample_ID == "P02_WG07"):
				family = "20-1010"
			else:
				family = "10-5_X_4-6664"
			new_line = short_sample_ID + "\t" + family + "\t"
			for depth_line in curr_mosdepth:
				region = depth_line.strip().split("\t")[0]
				print("region is:\t" + region) 
				mean_covg = depth_line.strip().split("\t")[3]
				print("mean_covg:\t" + mean_covg)
		
				new_line = new_line + mean_covg + "\t" 
				if(header_line_finished == 0):
					header_line = header_line + region + "_mean_covg\t"
			
			header_line_finished=1
			sample_lines.append(new_line)
			curr_mosdepth.close()
		
		all_sample = open(output.filename,'w')
		all_sample.write(header_line + "\n")
		for line in sample_lines:
			all_sample.write(line + "\n")
		all_sample.close()

rule compile_non_10_5_pooled_mosdepth_report_hisat2:
	input:
		non_10_5_pooled_mosdepths=expand("mosdepth/hisat2/{sample}.regions.bed",sample=config["non_10_5_pooled"]),
		non_10_5_pooled_list="/home/d.ence/projects/pinus_taeda_L/Fr1_project/aligning_Fr1_samples/aligning_to_custom_R_gene_V1_1_ref/calls/bam_list.non_10_5_pooled.txt"
	output:
		filename="mosdepth/hisat2/non_10_5_pooled.mosdepth.report.txt"
	log:
		filename="logs/mosdepth/hisat2/non_10_5_pooled.report.log"
	run:
		open_list = open(input.non_10_5_pooled_list,'r')
		header_line = "sample_ID\tfamily\t"
		header_line_finished=0
		sample_lines = []

		open_log = open(log.filename, 'w')
		for l in open_list:
			short_sample_ID = l.split('/')[10].split('_i')[0]
			long_sample_ID = l.split('/')[10].split('.')[0]

			curr_mosdepth = open("mosdepth/hisat2//" + long_sample_ID + ".regions.bed")
			open_log.write("opened this file:\tmosdepth/" + long_sample_ID + ".regions.bed\n")
			tmp = curr_mosdepth.readline()
			open_log.write(tmp + "\n")
		
			family=""
			open_log.write("short_sample_ID is:\t" + short_sample_ID + "\n")
			if(short_sample_ID == "P04_WA08"):
				family="CL04"
			elif(short_sample_ID == "P04_WB08"):
				family="CL05"
			elif(short_sample_ID == "P04_WC08"):
				family="PD18"
			elif(short_sample_ID == "P04_WD08"):
				family="PD35"
			new_line = short_sample_ID + "\t" + family + "\t"
			open_log.write("family is:\t" + family + "\n")
			
			for depth_line in curr_mosdepth:
				open_log.write("depth_line is:\t" + depth_line + "\n")
				region = depth_line.strip().split("\t")[0]
				mean_covg = depth_line.strip().split("\t")[3]
				
				open_log.write("region is:\t" + region + "\n")
				open_log.write("mean_covg is:\t" + mean_covg + "\n")

				if(header_line_finished == 0):
					header_line = header_line + region + "_mean_covg\t"
				new_line = new_line + mean_covg + "\t"
			sample_lines.append(new_line)	
			curr_mosdepth.close()
		
			header_line_finished=1
		all_sample = open(output.filename,'w')
		all_sample.write(header_line + "\n")
		for line in sample_lines:
			all_sample.write(line + "\n")
		all_sample.close()
		open_log.close() 

rule compile_10_5_megs_report_hisat2:
	input:
		megs_mosdepths=expand("mosdepth/hisat2/{sample}.regions.bed",sample=config["Fr1_meg_samples"]),
		megs_list="/home/d.ence/projects/pinus_taeda_L/Fr1_project/aligning_Fr1_samples/aligning_to_custom_R_gene_V1_1_ref/calls/bam_list.megs.txt"
	output:
		filename="mosdepth/hisat2/10_5_megs.mosdepth.report.txt"
	log:
		filename="logs/mosdepth/hisat2/10_5_megs.mosdepth.report.log"
	run:
		open_list = open(input.megs_list,'r')
		open_log = open(log.filename,'w')
		header_line = "sample_ID\tfamily\t"
		header_line_finished=0
		sample_lines = []

		for l in open_list:
			short_sample_ID = l.split('/')[10].split('_i')[0]
			long_sample_ID = l.split('/')[10].split('.')[0]

			curr_filename = "mosdepth/hisat2/" + long_sample_ID + ".regions.bed"
			curr_mosdepth = open(curr_filename,'r')
			tmp = curr_mosdepth.readline()
			
			new_line = short_sample_ID + "\t10_5_meg\t"
	
			line_count = 0	
			for depth_line in curr_mosdepth:
				region = depth_line.strip().split("\t")[0]
				mean_covg = depth_line.strip().split("\t")[3]
				
				if(header_line_finished == 0):
					header_line = header_line + region + "_mean_covg\t"
				new_line = new_line + mean_covg + "\t"
				line_count = line_count + 1
			sample_lines.append(new_line)
			curr_mosdepth.close()

			open_log.write("this file:\t" + curr_filename + " had this many lines:\t" + str(line_count) + "\n")
			header_line_finished = 1
		all_sample = open(output.filename,'w')
		all_sample.write(header_line + "\n")
		for line in sample_lines:
			all_sample.write(line + "\n")
		all_sample.close()
		open_log.close()

rule compile_10_5_prog_report_hisat2:
	input:
		prog_mosdepths=expand("mosdepth/hisat2/{sample}.regions.bed",sample=config["Fr1_prog"]),
		prog_list="/home/d.ence/projects/pinus_taeda_L/Fr1_project/aligning_Fr1_samples/aligning_to_custom_R_gene_V1_1_ref/calls/bam_list.Fr1_prog.txt"
	output:
		filename="mosdepth/hisat2/10_5_prog.mosdepth.report.txt"
	log:
		filename="logs/mosdepth/hisat2/10_5_prog.mosdepth.report.log"
	run:
		open_list = open(input.prog_list,'r')
		header_line = "sample_ID\tfamily\t"
		header_line_finished = 0
		sample_lines = []
		open_log = open(log.filename,'w')
	
		for l in open_list:
			short_sample_ID = l.split('/')[10].split('_i')[0]
			long_sample_ID = l.split('/')[10].split('.')[0]

			curr_filename = "mosdepth/hisat2/" + long_sample_ID + ".regions.bed"
			curr_mosdepth = open(curr_filename,'r')
			tmp = curr_mosdepth.readline()
			
			new_line = short_sample_ID + "\t10_5_prog\t"
		
			line_counter = 0	
			for depth_line in curr_mosdepth:
				region = depth_line.strip().split("\t")[0]
				mean_covg = depth_line.strip().split("\t")[3]
		
				if(header_line_finished == 0):	
					header_line = header_line + region + "_mean_covg\t"
				new_line = new_line + mean_covg + "\t"
				line_counter = line_counter + 1
			sample_lines.append(new_line)
			curr_mosdepth.close()
			open_log.write("this file:\t" + curr_filename + " had this many lines:\t" + str(line_counter) + "\n")
	
			header_line_finished = 1	
		all_sample = open(output.filename,'w')
		all_sample.write(header_line + "\n")
		for line in sample_lines:
			all_sample.write(line + "\n")
		all_sample.close()
rule report_alignment:
	input:
		prog_report_hisat2="mosdepth/hisat2/10_5_prog.mosdepth.report.txt",
		megs_report_hisat2="mosdepth/hisat2/10_5_megs.mosdepth.report.txt",
		non_10_5_indv_hisat2="mosdepth/hisat2/non_10_5_indv.mosdepth.report.txt",
		non_10_5_pooled_hisat2="mosdepth/hisat2/non_10_5_pooled.mosdepth.report.txt",
		prog_report_bwa_mem="mosdepth/bwa_mem/10_5_prog.mosdepth.report.txt",
		megs_report_bwa_mem="mosdepth/bwa_mem/10_5_megs.mosdepth.report.txt",
		non_10_5_indv_bwa_mem="mosdepth/bwa_mem/non_10_5_indv.mosdepth.report.txt",
		non_10_5_pooled_bwa_mem="mosdepth/bwa_mem/non_10_5_pooled.mosdepth.report.txt"
	output:
		"mosdepth.report.html"
	run:
		from snakemake.utils import report

		report("""
		A first test of a structural variant calling workflow
		=====================================================
		
		Reads were mapped to the P.nigra
		reference genome and variants were called jointly with
		gatk unified genotyper.

		This resulted variants (see Table T1_).
		Benchmark results for BWA can be found in the tables.
		""", output[0], T1=input)

rule freebayes_haploid_hisat2:
	input:
		bam=expand("RG_replaced_bams/{sample}.hisat2.sorted.rmdup.merged.bam",sample=config["Fr1_meg_samples"]),
		bai=expand("RG_replaced_bams/{sample}.hisat2.sorted.rmdup.merged.bam.bai",sample=config["Fr1_meg_samples"]),
		bam_list="/home/d.ence/projects/pinus_taeda_L/Fr1_project/aligning_Fr1_samples/aligning_to_custom_R_gene_V1_1_ref/calls/bam_list.hisat2.megs.txt",
		ref=config["reference"]["V1_01"]["custom"]
	log:
		"logs/freebayes/Fr1_megs.freebayes.log"
	output:
		"calls/freebayes/Fr1_megs.hisat2.freebayes.vcf"
	shell:
		"module load freebayes; freebayes-v1.3.1 --ploidy 1 -pvar 0.75 -theta 0.01 -indels -mnps -min-alternate-fraction 0.8 -min-alternate-count 1 --bam-list {input.bam_list} -f {input.ref} --vcf {output} --targets /home/d.ence/projects/pinus_taeda_L/Fr1_project/probe_files/RG_0610_merged.fixed.bed &> {log}"

rule freebayes_progeny_hisat2:
	input:
		bam=expand("RG_replaced_bams/{sample}.hisat2.sorted.rmdup.merged.bam",sample=config["Fr1_prog"]),
		bai=expand("RG_replaced_bams/{sample}.hisat2.sorted.rmdup.merged.bam.bai",sample=config["Fr1_prog"]),
		bam_list="/home/d.ence/projects/pinus_taeda_L/Fr1_project/aligning_Fr1_samples/aligning_to_custom_R_gene_V1_1_ref/calls/bam_list.Fr1_prog.hisat2.txt",
		ref=config["reference"]["V1_01"]["custom"]
	log:
		"logs/freebayes/Fr1_prog.freebayes.log"
	output:
		"calls/freebayes/Fr1_prog.hisat2.freebayes.vcf"
	shell:
		#"module load freebayes; freebayes-v1.3.1 --ploidy 2 --bam-list {input.bam_list} -f {input.ref} --vcf {output} --targets /home/d.ence/projects/pinus_taeda_L/Fr1_project/probe_files/RG_0610_merged.fixed.bed &> {log}"
		#"module load freebayes; freebayes-v1.3.1 --ploidy 2 --theta 0.01 --pvar 0.75 -indels -mnps --min-alternate-fraction 0.8 --min-alternate-count 1 --bam-list {input.bam_list} -f {input.ref} --vcf {output} --targets /home/d.ence/projects/pinus_taeda_L/Fr1_project/probe_files/RG_0610_merged.fixed.bed &> {log}"
		"module load freebayes; freebayes-v1.3.1 --ploidy 2 --theta 0.01 --pvar 0.75 --bam-list {input.bam_list} -f {input.ref} --vcf {output} --targets /home/d.ence/projects/pinus_taeda_L/Fr1_project/probe_files/RG_0610_merged.fixed.bed &> {log}"

rule freebayes_individual_hisat2:
	input:
		bam="RG_replaced_bams/{sample}.hisat2.sorted.rmdup.merged.bam",
		bai="RG_replaced_bams/{sample}.hisat2.sorted.rmdup.merged.bam.bai",
		ref=config["reference"]["V1_01"]["custom"]
	log:
		"logs/freebayes/individual/Fr1.{sample}.freebayes.hisat2.log"
	output:
		"calls/freebayes/individual/Fr1.{sample}.hisat2.freebayes.vcf"
	shell:
		"module load freebayes; freebayes-v1.3.1 --ploidy 2 -f {input.ref} --vcf {output} --targets /home/d.ence/projects/pinus_taeda_L/Fr1_project/probe_files/RG_0610_merged.fixed.bed {input.bam} &> {log}"

rule freebayes_individual_bwa_mem:
	input:
		bam="RG_replaced_bams/{sample}.bwa_mem.sorted.rmdup.merged.bam",
		bai="RG_replaced_bams/{sample}.bwa_mem.sorted.rmdup.merged.bam.bai",
		ref=config["reference"]["V1_01"]["custom"]
	log:
		"logs/freebayes/individual/Fr1.{sample}.freebayes.bwa_mem.log"
	output:
		"calls/freebayes/individual/Fr1.{sample}.bwa_mem.freebayes.vcf"
	shell:
		"module load freebayes; freebayes-v1.3.1 --ploidy 2 -f {input.ref} --vcf {output} --targets /home/d.ence/projects/pinus_taeda_L/Fr1_project/probe_files/RG_0610_merged.fixed.bed {input.bam} &> {log}"

rule freebayes_haploid_bwa_mem:
	input:
		bam=expand("RG_replaced_bams/{sample}.bwa_mem.sorted.rmdup.merged.bam",sample=config["Fr1_meg_samples"]),
		bai=expand("RG_replaced_bams/{sample}.bwa_mem.sorted.rmdup.merged.bam.bai",sample=config["Fr1_meg_samples"]),
		bam_list="/home/d.ence/projects/pinus_taeda_L/Fr1_project/aligning_Fr1_samples/aligning_to_custom_R_gene_V1_1_ref/calls/bam_list.megs.txt",
		ref=config["reference"]["V1_01"]["custom"]
	log:
		"logs/freebayes/Fr1_megs.freebayes.log"
	output:
		"calls/freebayes/Fr1_megs.bwa_mem.freebayes.vcf"
	shell:
		#"module load freebayes; freebayes-v1.3.1 --skip-coverage 1000 --ploidy 2  --min-coverage 50 --bam-list {input.bam_list} -f {input.ref} --vcf {output} --targets /home/d.ence/projects/pinus_taeda_L/Fr1_project/probe_files/RG_0610_merged.fixed.bed &> {log}"
		"module load freebayes; freebayes-v1.3.1 --ploidy 2 --bam-list {input.bam_list} -f {input.ref} --vcf {output} --targets /home/d.ence/projects/pinus_taeda_L/Fr1_project/probe_files/RG_0610_merged.fixed.bed &> {log}"

#rule freebayes_10_5_prog:
#	input:
#		bam=expand("RG_replaced_bams/{sample}.bwa_mem.sorted.rmdup.merged.bam",sample=config["Fr1_prog"]),
#		bai=expand("RG_replaced_bams/{sample}.bwa_mem.sorted.rmdup.merged.bam.bai",sample=config["Fr1_prog"]),
#		bam_list="/home/d.ence/projects/pinus_taeda_L/Fr1_project/aligning_Fr1_samples/aligning_to_custom_R_gene_V1_1_ref/calls/bam_list.Fr1_prog.txt",
#		ref=config["reference"]["V1_01"]["custom"],
#		deconv="calls/freebayes/Fr1_megs.bwa_mem.freebayes.deconvoluted.vcf.gz"
#	log:
#		"logs/freebayes/Fr1_prog.freebayes.log"
#	output:
#		"calls/freebayes/Fr1_prog.bwa_mem.freebayes.vcf"
#	shell:
#                "module load freebayes; freebayes-v1.3.1 --ploidy 2 --bam-list {input.bam_list} -f {input.ref} --vcf {output} --targets /home/d.ence/projects/pinus_taeda_L/Fr1_project/probe_files/RG_0610_merged.fixed.bed &> {log}"
                #"module load freebayes; freebayes-v1.3.1 --ploidy 2 --only-use-input-alleles {input.deconv} --bam-list {input.bam_list} -f {input.ref} --vcf {output} --targets /home/d.ence/projects/pinus_taeda_L/Fr1_project/probe_files/RG_0610_merged.fixed.bed &> {log}"
		#"module load freebayes; freebayes-v1.3.1 --ploidy 2 --variant-input {input.bgzipped} --min-coverage 50 --bam-list {input.bam_list} -f {input.ref} --vcf {output} --targets /home/d.ence/projects/pinus_taeda_L/Fr1_project/probe_files/RG_0610_merged.fixed.bed &> {log}"
		#"module load freebayes; freebayes-v1.3.1 --ploidy 2 --variant-input {input.bgzipped} --min-coverage 50 --bam-list {input.bam_list} -f {input.ref} --vcf {output} &> {log}"

rule freebayes_non_10_5_indv:
        input:
                bam=expand("RG_replaced_bams/{sample}.bwa_mem.sorted.rmdup.merged.bam",sample=config["non_10_5_indv"]),
                bai=expand("RG_replaced_bams/{sample}.bwa_mem.sorted.rmdup.merged.bam.bai",sample=config["non_10_5_indv"]),
                bam_list="/home/d.ence/projects/pinus_taeda_L/Fr1_project/aligning_Fr1_samples/aligning_to_custom_R_gene_V1_1_ref/calls/bam_list.non_10_5_indv.txt",
                ref=config["reference"]["V1_01"]["custom"]
        log:
                "logs/freebayes/non_10_5_indv.freebayes.log"
        output:
                "calls/freebayes/non_10_5_indv.bwa_mem.freebayes.vcf" 
        shell:
                #"module load freebayes; freebayes-v1.3.1 --ploidy 2 --skip-coverage 1000 --min-coverage 50 --bam-list {input.bam_list} -f {input.ref} --vcf {output} --targets /home/d.ence/projects/pinus_taeda_L/Fr1_project/probe_files/RG_0610_merged.fixed.bed &> {log}"
                "module load freebayes; freebayes-v1.3.1 --ploidy 2 --min-coverage 50 --bam-list {input.bam_list} -f {input.ref} --vcf {output} --targets /home/d.ence/projects/pinus_taeda_L/Fr1_project/probe_files/RG_0610_merged.fixed.bed &> {log}"

rule freebayes_non_10_5_indv_hisat2:
        input:
                bam=expand("RG_replaced_bams/{sample}.hisat2.sorted.rmdup.merged.bam",sample=config["non_10_5_indv"]),
                bai=expand("RG_replaced_bams/{sample}.hisat2.sorted.rmdup.merged.bam.bai",sample=config["non_10_5_indv"]),
                bam_list="/home/d.ence/projects/pinus_taeda_L/Fr1_project/aligning_Fr1_samples/aligning_to_custom_R_gene_V1_1_ref/calls/bam_list.non_10_5_indv.txt",
                ref=config["reference"]["V1_01"]["custom"]
        log:
                "logs/freebayes/non_10_5_indv.freebayes.log"
        output:
                "calls/freebayes/non_10_5_indv.hisat2.freebayes.vcf" 
        shell:
                #"module load freebayes; freebayes-v1.3.1 --ploidy 2 --skip-coverage 1000 --min-coverage 50 --bam-list {input.bam_list} -f {input.ref} --vcf {output} --targets /home/d.ence/projects/pinus_taeda_L/Fr1_project/probe_files/RG_0610_merged.fixed.bed &> {log}"
                "module load freebayes; freebayes-v1.3.1 --ploidy 2 --min-coverage 50 --bam-list {input.bam_list} -f {input.ref} --vcf {output} --targets /home/d.ence/projects/pinus_taeda_L/Fr1_project/probe_files/RG_0610_merged.fixed.bed &> {log}"

rule freebayes_non_10_5_pooled_single_sample:
	input:
                bam="RG_replaced_bams/{sample}.bwa_mem.sorted.rmdup.merged.bam",
                bai="RG_replaced_bams/{sample}.bwa_mem.sorted.rmdup.merged.bam.bai",
                ref=config["reference"]["V1_01"]["custom"]
	output:
		vcf="calls/freebayes/{sample}.non_10_5_pooled.bwa_mem.freebayes.vcf"
	log:
                "logs/freebayes/{sample}.non_10_5_pooled.freebayes.log"
	shell:
                #"module load freebayes; freebayes-v1.3.1 --ploidy 20 --skip-coverage 1000 --use-best-n-alleles 2 --pooled-discrete  --min-coverage 50 --bam {input.bam} -f {input.ref} --vcf {output.vcf} --targets /home/d.ence/projects/pinus_taeda_L/Fr1_project/probe_files/RG_0610_merged.fixed.bed &> {log}"
                "module load freebayes; freebayes-v1.3.1 --ploidy 20 --use-best-n-alleles 2 --pooled-discrete  --min-coverage 50 --bam {input.bam} -f {input.ref} --vcf {output.vcf} --targets /home/d.ence/projects/pinus_taeda_L/Fr1_project/probe_files/RG_0610_merged.fixed.bed &> {log}"

rule freebayes_non_10_5_pooled_single_sample_hisat2:
	input:
                bam="RG_replaced_bams/{sample}.hisat2.sorted.rmdup.merged.bam",
                bai="RG_replaced_bams/{sample}.hisat2.sorted.rmdup.merged.bam.bai",
                ref=config["reference"]["V1_01"]["custom"]
	output:
		vcf="calls/freebayes/{sample}.non_10_5_pooled.hisat2.freebayes.vcf"
	log:
                "logs/freebayes/{sample}.non_10_5_pooled.freebayes.log"
	shell:
                #"module load freebayes; freebayes-v1.3.1 --ploidy 20 --skip-coverage 1000 --use-best-n-alleles 2 --pooled-discrete  --min-coverage 50 --bam {input.bam} -f {input.ref} --vcf {output.vcf} --targets /home/d.ence/projects/pinus_taeda_L/Fr1_project/probe_files/RG_0610_merged.fixed.bed &> {log}"
                "module load freebayes; freebayes-v1.3.1 --ploidy 20 --use-best-n-alleles 2 --pooled-discrete  --min-coverage 50 --bam {input.bam} -f {input.ref} --vcf {output.vcf} --targets /home/d.ence/projects/pinus_taeda_L/Fr1_project/probe_files/RG_0610_merged.fixed.bed &> {log}"

rule bgzip_and_tabix_megs_vcf:
	input:
		megs_vcf="calls/freebayes/Fr1_megs.bwa_mem.freebayes.deconvoluted.vcf"
	output:
		bgzipped="calls/freebayes/Fr1_megs.bwa_mem.freebayes.deconvoluted.vcf.gz",
		tabix_indexed="calls/freebayes/Fr1_megs.bwa_mem.freebayes.deconvoluted.vcf.gz.tbi"
	shell:
		"module load vcflib; cat {input.megs_vcf} | bgzip -c > {output.bgzipped} ; tabix -p vcf {output.bgzipped} ; touch {output.tabix_indexed}"

rule bgzip_and_tabix_individual_hisat2_vcf:
	input:
		vcf="calls/freebayes/individual/Fr1.{sample}.hisat2.freebayes.vcf"
	output:
		bgzipped="calls/freebayes/individual/Fr1.{sample}.hisat2.freebayes.vcf.gz",
		tabix_indexed="calls/freebayes/individual/Fr1.{sample}.hisat2.freebayes.vcf.gz.tbi"
	shell:
		"module load vcflib; cat {input.vcf} | bgzip -c > {output.bgzipped} ; tabix -p vcf {output.bgzipped} ; touch {output.tabix_indexed}"
	
rule bgzip_and_tabix_individual_bwa_mem_vcf:
	input:
		vcf="calls/freebayes/individual/Fr1.{sample}.bwa_mem.freebayes.vcf"
	output:
		bgzipped="calls/freebayes/individual/Fr1.{sample}.bwa_mem.freebayes.vcf.gz",
		tabix_indexed="calls/freebayes/individual/Fr1.{sample}.bwa_mem.freebayes.vcf.gz.tbi"
	shell:
		"module load vcflib; cat {input.vcf} | bgzip -c > {output.bgzipped} ; tabix -p vcf {output.bgzipped} ; touch {output.tabix_indexed}"

rule bgzip_and_tabix_vcf:
	input:
		vcf="calls/freebayes/{sample}.non_10_5_pooled.bwa_mem.freebayes.vcf"
	output:
		bgzipped="calls/freebayes/{sample}.non_10_5_pooled.bwa_mem.freebayes.vcf.gz",
		tabix_indexed="calls/freebayes/{sample}.non_10_5_pooled.bwa_mem.freebayes.vcf.gz.tbi"
	shell:
		"module load vcflib; cat {input.vcf} | bgzip -c > {output.bgzipped} ; tabix -p vcf {output.bgzipped} ; touch {output.tabix_indexed}"

rule bgzip_and_tabix_vcf_hisat2:
	input:
		vcf="calls/freebayes/{sample}.non_10_5_pooled.hisat2.freebayes.vcf"
	output:
		bgzipped="calls/freebayes/{sample}.non_10_5_pooled.hisat2.freebayes.vcf.gz",
		tabix_indexed="calls/freebayes/{sample}.non_10_5_pooled.hisat2.freebayes.vcf.gz.tbi"
	shell:
		"module load vcflib; cat {input.vcf} | bgzip -c > {output.bgzipped} ; tabix -p vcf {output.bgzipped} ; touch {output.tabix_indexed}"

rule postprocess_megs:
	input:
		megs_vcf="calls/freebayes/Fr1_megs.bwa_mem.freebayes.vcf",
		ref=config["reference"]["V1_01"]["custom"]
	output:
		megs_deconv_vcf="calls/freebayes/Fr1_megs.bwa_mem.freebayes.deconvoluted.vcf"
	shell:
		"module load bedtools vcflib vt; cat {input.megs_vcf} | vcfallelicprimitives --keep-info --keep-geno | vt normalize -r {input.ref} - | perl -ane 'if(/CHROM/){{ s/_i5\S+//g}} print $_' > {output.megs_deconv_vcf}"

rule postprocess_megs_hisat2:
	input:
		megs_vcf="calls/freebayes/Fr1_megs.hisat2.freebayes.vcf",
		ref=config["reference"]["V1_01"]["custom"]
	output:
		megs_deconv_vcf="calls/freebayes/Fr1_megs.hisat2.freebayes.deconvoluted.vcf"
	shell:
		"module load bedtools vcflib vt; cat {input.megs_vcf} | vcfallelicprimitives --keep-info --keep-geno | vt normalize -r {input.ref} - | perl -ane 'if(/CHROM/){{ s/_i5\S+//g}} print $_' > {output.megs_deconv_vcf}"

rule postprocess_progs_hisat2:
	input:
		prog_vcf="calls/freebayes/Fr1_prog.hisat2.freebayes.vcf",
		ref=config["reference"]["V1_01"]["custom"]
	output:
		prog_deconv_vcf="calls/freebayes/Fr1_prog.hisat2.freebayes.deconvoluted.vcf"
	shell:
		"module load bedtools vcflib vt; cat {input.prog_vcf} | vcfallelicprimitives --keep-info --keep-geno | vt normalize -r {input.ref} - | perl -ane 'if(/CHROM/){{ s/_i5\S+//g}} print $_' > {output.prog_deconv_vcf};" 

rule postprocess_progs:
	input:
		prog_vcf="calls/freebayes/Fr1_prog.bwa_mem.freebayes.vcf",
		ref=config["reference"]["V1_01"]["custom"]
	output:
		prog_deconv_vcf="calls/freebayes/Fr1_prog.bwa_mem.freebayes.deconvoluted.vcf"
	shell:
		"module load bedtools vcflib vt; cat {input.prog_vcf} | vcfallelicprimitives --keep-info --keep-geno | vt normalize -r {input.ref} - | perl -ane 'if(/CHROM/){{ s/_i5\S+//g}} print $_' > {output.prog_deconv_vcf};" 

rule postprocess_other_vcfs:
	input:
		non_10_5_pooled_vcf="calls/freebayes/non_10_5_pooled.bwa_mem.freebayes.vcf",
		non_10_5_indv_vcf="calls/freebayes/non_10_5_indv.bwa_mem.freebayes.vcf",
		ref=config["reference"]["V1_01"]["custom"]
	output:
		non_10_5_pooled_deconv_vcf="calls/freebayes/non_10_5_pooled.bwa_mem.freebayes.deconvoluted.vcf",
		non_10_5_indv_deconv_vcf="calls/freebayes/non_10_5_indv.bwa_mem.freebayes.deconvoluted.vcf"
	shell:
		"module load bedtools vcflib vt; cat {input.non_10_5_pooled_vcf} | vcfallelicprimitives --keep-info --keep-geno | vt normalize -r {input.ref} - | perl -ane 'if(/CHROM/){{ s/_i5\S+//g}} print $_' > {output.non_10_5_pooled_deconv_vcf}; cat {input.non_10_5_indv_vcf} | vcfallelicprimitives --keep-info --keep-geno | vt normalize -r {input.ref} - | perl -ane 'if(/CHROM/){{ s/_i5\S+//g}} print $_' > {output.non_10_5_indv_deconv_vcf}"
		
rule postprocess_other_vcfs_hisat2:
	input:
		non_10_5_pooled_vcf="calls/freebayes/non_10_5_pooled.hisat2.freebayes.vcf",
		non_10_5_indv_vcf="calls/freebayes/non_10_5_indv.hisat2.freebayes.vcf",
		ref=config["reference"]["V1_01"]["custom"]
	output:
		non_10_5_pooled_deconv_vcf="calls/freebayes/non_10_5_pooled.hisat2.freebayes.deconvoluted.vcf",
		non_10_5_indv_deconv_vcf="calls/freebayes/non_10_5_indv.hisat2.freebayes.deconvoluted.vcf"
	shell:
		"module load bedtools vcflib vt; cat {input.non_10_5_pooled_vcf} | vcfallelicprimitives --keep-info --keep-geno | vt normalize -r {input.ref} - | perl -ane 'if(/CHROM/){{ s/_i5\S+//g}} print $_' > {output.non_10_5_pooled_deconv_vcf}; cat {input.non_10_5_indv_vcf} | vcfallelicprimitives --keep-info --keep-geno | vt normalize -r {input.ref} - | perl -ane 'if(/CHROM/){{ s/_i5\S+//g}} print $_' > {output.non_10_5_indv_deconv_vcf}"

rule merge_freebayes_non_10_5_pooled_singles_bwa_mem:
	input:
		bgzipped_vcf=expand("calls/freebayes/{sample}.non_10_5_pooled.bwa_mem.freebayes.vcf.gz",sample=config["non_10_5_pooled"])
	output:
		vcf="calls/freebayes/non_10_5_pooled.bwa_mem.freebayes.vcf"
	log:
		"logs/freebayes/non_10_5_pooled.freebayes.log"
	shell:
		"module load vcftools; vcf-merge -d --ref-for-missing 0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0 {input.bgzipped_vcf} > {output.vcf} 2> {log}"

rule merge_freebayes_non_10_5_pooled_singles_hisat2:
	input:
		bgzipped_vcf=expand("calls/freebayes/{sample}.non_10_5_pooled.hisat2.freebayes.vcf.gz",sample=config["non_10_5_pooled"])
	output:
		vcf="calls/freebayes/non_10_5_pooled.hisat2.freebayes.vcf"
	log:
		"logs/freebayes/non_10_5_pooled.hisat2.freebayes.log"
	shell:
		"module load vcftools; vcf-merge -d --ref-for-missing 0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0 {input.bgzipped_vcf} > {output.vcf} 2> {log}"

#rule merge_freebayes_Fr1_prog_bwa_mem:
#	input:
#		bgzipped_vcf=expand("calls/freebayes/individual/Fr1.{sample}.bwa_mem.freebayes.vcf.gz",sample=config["Fr1_prog"])
#	output:
#		vcf="calls/freebayes/Fr1_prog.bwa_mem.freebayes.vcf"
#	log:
#		"logs/freebayes/Fr1_prog.bwa_mem.freebayes.log"
#	shell:
#		"module load vcftools; vcf-merge -d --ref-for-missing 0/0 {input.bgzipped_vcf} > {output.vcf} 2> {log}"

rule merge_freebayes_Fr1_prog_bwa_mem:
	input:
		bgzipped_vcf=expand("calls/freebayes/individual/Fr1.{sample}.bwa_mem.freebayes.vcf.gz",sample=config["Fr1_prog"])
	output:
		vcf="calls/freebayes/Fr1_prog.bwa_mem.freebayes.vcf"
	log:
		"logs/freebayes/Fr1_prog.bwa_mem.freebayes.log"
	shell:
		"module load bcftools; bcftools merge -0 {input.bgzipped_vcf} --output {output} --output-type v" 

#rule merge_freebayes_Fr1_prog_hisat2:
#	input:
#		bgzipped_vcf=expand("calls/freebayes/individual/Fr1.{sample}.hisat2.freebayes.vcf.gz",sample=config["Fr1_prog"])
#	output:
#		vcf="calls/freebayes/Fr1_prog.hisat2.freebayes.vcf"
#	log:
#		"logs/freebayes/Fr1_prog.hisat2.freebayes.log"
#	shell:
#		"module load vcftools; vcf-merge -d --ref-for-missing 0/0 {input.bgzipped_vcf} > {output.vcf} 2> {log}"

#rule merge_freebayes_Fr1_prog_hisat2:
#	input:
#		bgzipped_vcf=expand("calls/freebayes/individual/Fr1.{sample}.hisat2.freebayes.vcf.gz",sample=config["Fr1_prog"])
#	output:
#		vcf="calls/freebayes/Fr1_prog.hisat2.freebayes.vcf"
#	log:
#		"logs/freebayes/Fr1_prog.hisat2.freebayes.log"
#	shell:
#		"module load bcftools; bcftools merge -0 {input.bgzipped_vcf} --output {output} --output-type v"

#rule freebayes_non_10_5_pooled:
#        input:
#                bam=expand("RG_replaced_bams/{sample}.bwa_mem.sorted.rmdup.realigned.merged.bam",sample=config["non_10_5_pooled"]),
#                bai=expand("RG_replaced_bams/{sample}.bwa_mem.sorted.rmdup.realigned.merged.bam.bai",sample=config["non_10_5_pooled"]),
#                bam_list="/home/d.ence/projects/pinus_taeda_L/Fr1_project/aligning_Fr1_samples/aligning_to_custom_R_gene_V1_1_ref/calls/bam_list.non_10_5_pooled.txt",
#                ref=config["reference"]["V1_01"]["custom"]
#        log:
#                "logs/freebayes/non_10_5_pooled.freebayes.log"
#        output:
#                "calls/freebayes/non_10_5_pooled.bwa_mem.freebayes.vcf"
#        shell:
#                "module load freebayes; freebayes-v1.3.1 --ploidy 20  --use-best-n-alleles 20 --pooled-discrete  --min-coverage 50 --bam-list {input.bam_list} -f {input.ref} --vcf {output} --targets /home/d.ence/projects/pinus_taeda_L/Fr1_project/probe_files/RG_0610_merged.fixed.bed &> {log}"

rule report_vcfs:
	input:
		megs_vcf="calls/freebayes/Fr1_megs.bwa_mem.freebayes.deconvoluted.vcf",
		prog_vcf="calls/freebayes/Fr1_prog.bwa_mem.freebayes.deconvoluted.vcf",
		non_10_5_indv_vcf="calls/freebayes/non_10_5_indv.bwa_mem.freebayes.deconvoluted.vcf",
		non_10_5_pooled_vcf="calls/freebayes/non_10_5_pooled.bwa_mem.freebayes.deconvoluted.vcf",
		megs_hisat2="calls/freebayes/Fr1_megs.hisat2.freebayes.deconvoluted.vcf",
		prog_hisat2="calls/freebayes/Fr1_prog.hisat2.freebayes.deconvoluted.vcf",
		non_10_5_indv_hisat2_vcf="calls/freebayes/non_10_5_indv.hisat2.freebayes.deconvoluted.vcf",
		non_10_5_pooled_hisat2_vcf="calls/freebayes/non_10_5_pooled.hisat2.freebayes.deconvoluted.vcf"
			
	output:
		"freebayes_calls.report.html"
	run:
		from snakemake.utils import report
		with open(input.megs_vcf) as vcf:
			n_calls = sum(1 for l in vcf if not l.startswith("#"))

		report("""
		A first test of a structural variant calling workflow
		=====================================================
		
		Reads were mapped to the P.nigra
		reference genome and variants were called jointly with
		gatk unified genotyper.

		This resulted in {n_calls} variants (see Table T1_).
		Benchmark results for BWA can be found in the tables.
		""", output[0], T1=input[0])

rule break_progeny_into_groups_filter_qual:
	input:
		prog_vcf="calls/freebayes/Fr1_prog.hisat2.freebayes.deconvoluted.vcf",
		rest_list=config["args_lists"]["susc"],
		susc_list=config["args_lists"]["rest"]
	output:
		susc_vcf="filtering_calls_March_2020_hisat2/Fr1_susceptible_samples.hisat2.fixed_sample_names.deconvoluted.recode.vcf",
		rest_vcf="filtering_calls_March_2020_hisat2/Fr1_resistant_samples.hisat2.fixed_sample_names.deconvoluted.recode.vcf",
	params:
		susc_out="filtering_calls_March_2020_hisat2/Fr1_susceptible_samples.hisat2.fixed_sample_names.deconvoluted",
		rest_out="filtering_calls_March_2020_hisat2/Fr1_resistant_samples.hisat2.fixed_sample_names.deconvoluted"
	shell:
		"module load vcftools; vcftools --max-alleles 2 --minQ 20 --vcf {input.prog_vcf} --keep {input.susc_list} --recode --recode-INFO-all --out {params.susc_out}; vcftools --max-alleles 2 --minQ 20 --vcf {input.prog_vcf} --keep {input.rest_list} --recode --recode-INFO-all --out {params.rest_out}; "
		
rule find_rare_susc_variants:
	input:
		"filtering_calls_March_2020_hisat2/Fr1_susceptible_samples.hisat2.fixed_sample_names.deconvoluted.recode.vcf"
	output:
		thresh_zero_susc="filtering_calls_March_2020_hisat2/Fr1_susceptible_samples.hisat2.MAF_MAX_0.00.fixed_sample_names.recode.vcf",
		thresh_one_susc="filtering_calls_March_2020_hisat2/Fr1_susceptible_samples.hisat2.MAF_MAX_0.01.fixed_sample_names.recode.vcf",
		thresh_five_susc="filtering_calls_March_2020_hisat2/Fr1_susceptible_samples.hisat2.MAF_MAX_0.05.fixed_sample_names.recode.vcf",
		thresh_ten_susc="filtering_calls_March_2020_hisat2/Fr1_susceptible_samples.hisat2.MAF_MAX_0.10.fixed_sample_names.recode.vcf"
	params:
		thresh_zero_out="filtering_calls_March_2020_hisat2/Fr1_susceptible_samples.hisat2.MAF_MAX_0.00.fixed_sample_names",
		thresh_one_out="filtering_calls_March_2020_hisat2/Fr1_susceptible_samples.hisat2.MAF_MAX_0.01.fixed_sample_names",
		thresh_five_out="filtering_calls_March_2020_hisat2/Fr1_susceptible_samples.hisat2.MAF_MAX_0.05.fixed_sample_names",
		thresh_ten_out="filtering_calls_March_2020_hisat2/Fr1_susceptible_samples.hisat2.MAF_MAX_0.10.fixed_sample_names"
	shell:
		"module load vcftools; vcftools --recode --recode-INFO-all --vcf {input} --minQ 20 --max-maf 0.00 --out {params.thresh_zero_out}; vcftools --recode --recode-INFO-all --vcf {input} --minQ 20 --max-maf 0.01 --out {params.thresh_one_out}; vcftools --recode --recode-INFO-all --vcf {input} --minQ 20 --max-maf 0.05 --out {params.thresh_five_out};  vcftools --recode --recode-INFO-all --vcf {input} --minQ 20 --max-maf 0.10 --out {params.thresh_ten_out};"

rule filter_megs_het_errors_missing_data:
	input:
		"calls/freebayes/Fr1_megs.hisat2.freebayes.deconvoluted.vcf"
	output:
		"filtering_calls_March_2020_hisat2/Fr1_megs.hisat2.freebayes.filtered_missing.biallelic.recode.vcf"
	log:
		"logs/filtering_calls_March_2020.log"
	params:
		out="filtering_calls_March_2020_hisat2/Fr1_megs.hisat2.freebayes.filtered_missing.biallelic",
		no_het_vcf="filtering_calls_March_2020_hisat2/megs.no_het.vcf",
		imiss_missing_prehet="filtering_calls_March_2020_hisat2/Fr1_megs.hisat2.freebayes.filtered_missing.biallelic.imiss",
		out_no_het="filtering_calls_March_2020_hisat2/megs.no_het",
		imiss_missing_post_het="filtering_calls_March_2020_hisat2/megs.no_het.imiss",
		missing_list="filtering_calls_March_2020_hisat2/megs.indvs_missing_data.txt"
	shell:
		"""
		module load vcftools; vcftools --vcf {input} --missing-indv --out {params.out} 2>> {log}
		cat {params.imiss_missing_prehet} | grep -v "F_MISS" | perl -ane 'chomp; my ($indv, $a, $b, $c, $F_MISS) = split(/\t/); if($F_MISS > 0.2){{print "$indv\n"}}' > {params.missing_list} 2>> {log}
		module load bcftools; bcftools +setGT {input} -- -t q -i 'GT="het"' -n "./." > {params.no_het_vcf} 2>> {log}
		module load vcftools; vcftools --vcf {params.no_het_vcf} --missing-indv --out {params.out_no_het} 2>> {log}
		cat {params.imiss_missing_post_het} | grep -v "F_MISS" | perl -ane 'chomp; my ($indv, $a, $b, $c, $F_MISS) = split(/\t/); if($F_MISS > 0.4){{ print "$indv\n"}}' >> {params.missing_list} 2>> {log}
		module load vcftools; vcftools --vcf {params.no_het_vcf} --recode --recode-INFO-all --max-alleles 2 --minQ 20 --remove {params.missing_list} --out {params.out} 2>> {log}
		"""
rule intersect_megagams_and_rare_susceptible:
	input:
		#megs_vcf=config["args_lists"]["chi_filtered_megs"],
		megs_vcf="filtering_calls_March_2020_hisat2/Fr1_megs.hisat2.freebayes.filtered_missing.biallelic.recode.vcf",
		thresh_zero_susc="filtering_calls_March_2020_hisat2/Fr1_susceptible_samples.hisat2.MAF_MAX_0.00.fixed_sample_names.recode.vcf",
		thresh_one_susc="filtering_calls_March_2020_hisat2/Fr1_susceptible_samples.hisat2.MAF_MAX_0.01.fixed_sample_names.recode.vcf",
		thresh_five_susc="filtering_calls_March_2020_hisat2/Fr1_susceptible_samples.hisat2.MAF_MAX_0.05.fixed_sample_names.recode.vcf",
		thresh_ten_susc="filtering_calls_March_2020_hisat2/Fr1_susceptible_samples.hisat2.MAF_MAX_0.10.fixed_sample_names.recode.vcf"
	output:
		megs_rare_susc_zero="filtering_calls_March_2020_hisat2/Fr1_susceptible_samples.megagam_cands.MAF_MAX_0.00.fixed_sample_names.vcf",
		megs_rare_susc_one="filtering_calls_March_2020_hisat2/Fr1_susceptible_samples.megagam_cands.MAF_MAX_0.01.fixed_sample_names.vcf",
		megs_rare_susc_five="filtering_calls_March_2020_hisat2/Fr1_susceptible_samples.megagam_cands.MAF_MAX_0.05.fixed_sample_names.vcf",
		megs_rare_susc_ten="filtering_calls_March_2020_hisat2/Fr1_susceptible_samples.megagam_cands.MAF_MAX_0.10.fixed_sample_names.vcf"	
	shell:
		"module load bedtools; bedtools intersect -header -wa -a {input.megs_vcf} -b {input.thresh_zero_susc} > {output.megs_rare_susc_zero}; bedtools intersect -header -wa -a {input.megs_vcf} -b {input.thresh_one_susc} > {output.megs_rare_susc_one}; bedtools intersect -header -wa -a {input.megs_vcf} -b {input.thresh_five_susc} > {output.megs_rare_susc_five}; bedtools intersect -header -wa -a {input.megs_vcf} -b {input.thresh_ten_susc} > {output.megs_rare_susc_ten};"

rule find_frequent_rest_variants:
	input:
		"filtering_calls_March_2020_hisat2/Fr1_resistant_samples.hisat2.fixed_sample_names.deconvoluted.recode.vcf"
	output:
		"filtering_calls_March_2020_hisat2/Fr1_resistant_samples.hisat2.MAF_0.10.fixed_sample_names.recode.vcf"
	params:
		"filtering_calls_March_2020_hisat2/Fr1_resistant_samples.hisat2.MAF_0.10.fixed_sample_names"	
	shell:
		"module load vcftools; vcftools --recode --recode-INFO-all --vcf {input} --minQ 20 --maf 0.10 --max-maf 0.5 --out {params}"	

rule intersect_rare_susceptible_and_frequent_rest:
	input:
		frequent_rest="filtering_calls_March_2020_hisat2/Fr1_resistant_samples.hisat2.MAF_0.10.fixed_sample_names.recode.vcf",
		thresh_zero_susc="filtering_calls_March_2020_hisat2/Fr1_susceptible_samples.hisat2.MAF_MAX_0.00.fixed_sample_names.recode.vcf",
		thresh_one_susc="filtering_calls_March_2020_hisat2/Fr1_susceptible_samples.hisat2.MAF_MAX_0.01.fixed_sample_names.recode.vcf",
		thresh_five_susc="filtering_calls_March_2020_hisat2/Fr1_susceptible_samples.hisat2.MAF_MAX_0.05.fixed_sample_names.recode.vcf",
		thresh_ten_susc="filtering_calls_March_2020_hisat2/Fr1_susceptible_samples.hisat2.MAF_MAX_0.10.fixed_sample_names.recode.vcf"
	output:
		both_rare_zero_frequent_rest="filtering_calls_March_2020_hisat2/Fr1_susceptible_samples.MAX_MAF_0.00.MAF_0.10_in_rest.recode.vcf",
		both_rare_one_frequent_rest="filtering_calls_March_2020_hisat2/Fr1_susceptible_samples.MAX_MAF_0.01.MAF_0.10_in_rest.recode.vcf",
		both_rare_five_frequent_rest="filtering_calls_March_2020_hisat2/Fr1_susceptible_samples.MAX_MAF_0.05.MAF_0.10_in_rest.recode.vcf",
		both_rare_ten_frequent_rest="filtering_calls_March_2020_hisat2/Fr1_susceptible_samples.MAX_MAF_0.10.MAF_0.10_in_rest.recode.vcf"
	shell:
		"module load bedtools; bedtools intersect -header -wa -a {input.thresh_zero_susc} -b {input.frequent_rest} > {output.both_rare_zero_frequent_rest}; bedtools intersect -header -wa -a {input.thresh_one_susc} -b {input.frequent_rest} > {output.both_rare_one_frequent_rest}; bedtools intersect -header -wa -a {input.thresh_five_susc} -b {input.frequent_rest} > {output.both_rare_five_frequent_rest}; bedtools intersect -header -wa -a {input.thresh_ten_susc} -b {input.frequent_rest} > {output.both_rare_ten_frequent_rest}; "
rule intersect_megs_rare_susc_frequent_rest:
	input:
		megs_rare_susc_zero="filtering_calls_March_2020_hisat2/Fr1_susceptible_samples.megagam_cands.MAF_MAX_0.00.fixed_sample_names.vcf",
		megs_rare_susc_one="filtering_calls_March_2020_hisat2/Fr1_susceptible_samples.megagam_cands.MAF_MAX_0.01.fixed_sample_names.vcf",
		megs_rare_susc_five="filtering_calls_March_2020_hisat2/Fr1_susceptible_samples.megagam_cands.MAF_MAX_0.05.fixed_sample_names.vcf",
		megs_rare_susc_ten="filtering_calls_March_2020_hisat2/Fr1_susceptible_samples.megagam_cands.MAF_MAX_0.10.fixed_sample_names.vcf",
		frequent_rest="filtering_calls_March_2020_hisat2/Fr1_resistant_samples.hisat2.MAF_0.10.fixed_sample_names.recode.vcf"
	output:
		megs_rare_susc_zero_freq_rest="filtering_calls_March_2020_hisat2/segr_megagam.rare_susc_0.00.freq_rest.vcf",
		megs_rare_susc_one_freq_rest="filtering_calls_March_2020_hisat2/segr_megagam.rare_susc_0.01.freq_rest.vcf",
		megs_rare_susc_five_freq_rest="filtering_calls_March_2020_hisat2/segr_megagam.rare_susc_0.05.freq_rest.vcf",
		megs_rare_susc_ten_freq_rest="filtering_calls_March_2020_hisat2/segr_megagam.rare_susc_0.10.freq_rest.vcf"
	shell:
		"module load bedtools; bedtools intersect -header -wa -a {input.frequent_rest} -b {input.megs_rare_susc_zero} > {output.megs_rare_susc_zero_freq_rest}; bedtools intersect -header -wa -a {input.frequent_rest} -b {input.megs_rare_susc_one} > {output.megs_rare_susc_one_freq_rest}; bedtools intersect -header -wa -a {input.frequent_rest} -b {input.megs_rare_susc_five} > {output.megs_rare_susc_five_freq_rest}; bedtools intersect -header -wa -a {input.frequent_rest} -b {input.megs_rare_susc_ten} > {output.megs_rare_susc_ten_freq_rest}"

rule check_finished:
	input:
		both_rare_zero_frequent_rest="filtering_calls_March_2020_hisat2/Fr1_susceptible_samples.MAX_MAF_0.00.MAF_0.10_in_rest.recode.vcf",
		both_rare_one_frequent_rest="filtering_calls_March_2020_hisat2/Fr1_susceptible_samples.MAX_MAF_0.01.MAF_0.10_in_rest.recode.vcf",
		both_rare_five_frequent_rest="filtering_calls_March_2020_hisat2/Fr1_susceptible_samples.MAX_MAF_0.05.MAF_0.10_in_rest.recode.vcf",
		both_rare_ten_frequent_rest="filtering_calls_March_2020_hisat2/Fr1_susceptible_samples.MAX_MAF_0.10.MAF_0.10_in_rest.recode.vcf",
		megs_rare_susc_zero="filtering_calls_March_2020_hisat2/Fr1_susceptible_samples.megagam_cands.MAF_MAX_0.00.fixed_sample_names.vcf",
		megs_rare_susc_one="filtering_calls_March_2020_hisat2/Fr1_susceptible_samples.megagam_cands.MAF_MAX_0.01.fixed_sample_names.vcf",
		megs_rare_susc_five="filtering_calls_March_2020_hisat2/Fr1_susceptible_samples.megagam_cands.MAF_MAX_0.05.fixed_sample_names.vcf",
		megs_rare_susc_ten="filtering_calls_March_2020_hisat2/Fr1_susceptible_samples.megagam_cands.MAF_MAX_0.10.fixed_sample_names.vcf",
		megs_rare_susc_zero_freq_rest="filtering_calls_March_2020_hisat2/segr_megagam.rare_susc_0.00.freq_rest.vcf",
		megs_rare_susc_one_freq_rest="filtering_calls_March_2020_hisat2/segr_megagam.rare_susc_0.01.freq_rest.vcf",
		megs_rare_susc_five_freq_rest="filtering_calls_March_2020_hisat2/segr_megagam.rare_susc_0.05.freq_rest.vcf",
		megs_rare_susc_ten_freq_rest="filtering_calls_March_2020_hisat2/segr_megagam.rare_susc_0.10.freq_rest.vcf"
		
	output:
		"touched_successfully.txt"
	shell:
		"touch {input.both_rare_zero_frequent_rest} > {output}; touch {input.both_rare_one_frequent_rest} > {output}; touch {input.both_rare_five_frequent_rest} > {output}; touch {input.both_rare_ten_frequent_rest} > {output}; touch {input.megs_rare_susc_zero} > {output}; touch {input.megs_rare_susc_one} > {output}; touch {input.megs_rare_susc_five} > {output}; touch {input.megs_rare_susc_ten} > {output}; touch {input.megs_rare_susc_zero_freq_rest} > {output}; touch {input.megs_rare_susc_one_freq_rest} > {output}; touch {input.megs_rare_susc_five_freq_rest} > {output}; touch {input.megs_rare_susc_ten_freq_rest} > {output}"



 	
