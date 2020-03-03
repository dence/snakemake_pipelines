rule all:
	input:
		"freebayes_calls.report.html"

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

rule hisat2_align:
	input:
		fq1="trimmed_reads/{sample}/{sample}.trimmed.R1.fastq.gz",
		fq2="trimmed_reads/{sample}/{sample}.trimmed.R2.fastq.gz"
	params:
		sample="{sample}",
		ref=config["reference"]["V1_01"]["custom"]
	output:
		temp("mapped_reads/{sample}.hisat2.bam")
	log:
		"logs/hisat2/{sample}.bwa_mem.log"
	shell:
		"module load bwa/0.7.17 ; bwa mem {params.ref} " 
		+ "-R '@RG\\tID:{params.sample}\\tSM:{params.sample}\\tPL:ILLUMINA'" 
		+ "{input.fq1} {input.fq2} > {output} 2> {log}"	

rule samtools_sort:
	input:
		"mapped_reads/{sample}.hisat2.bam"
	output:
		temp("sorted_reads/{sample}.hisat2.sorted.bam")
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

rule samtools_rmdup:
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

rule samtools_index_rmduped:
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

rule samtools_index_sorted:
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

rule samtools_index_sorted_hisat2:
	input:
		"sorted_reads/{sample}.hisat2.sorted.rmdup.bam"
	output:
		"sorted_reads/{sample}.hisat2.sorted.rmdup.bam.bai"
	log:
		"logs/samtools_index_sorted/{sample}.log"
	benchmark:
		"benchmarks/{sample}.index_sorted.benchmark.txt"	
	shell:
		"module load samtools; samtool index {input}"
rule gatk_indel_creator:
	input:
		bam="rmduped_reads/{sample}.hisat2.sorted.rmdup.bam",
		bai="rmduped_reads/{sample}.hisat2.sorted.rmdup.bam.bai"
	output:
		temp("realigner_intervals/{sample}.hisat2.intervals")
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
		"realigner_intervals/{sample}.hisat2.intervals"
	params:
		ref=config["reference"]["V1_01"]["custom"]
	log:
		"logs/realigner_intervals/{sample}.intervals.hisat2.log"
	benchmark:
		"benchmarks/{sample}.intervals.benchmark.txt"
	shell:
		"module load gatk;  java -jar -Xmx10g /apps/gatk/3.7.0/GenomeAnalysisTK.jar -T RealignerTargetCreator --filter_reads_with_N_cigar -I {input.bam} -o {output} -R {params.ref} &> {log}"

rule gatk_indel_realign:
	input:
		bam="rmduped_reads/{sample}.hisat2.sorted.rmdup.bam",
		bai="rmduped_reads/{sample}.hisat2.sorted.rmdup.bam.bai",
		interval="realigner_intervals/{sample}.hisat2.intervals"
	output:
		"realigned_bams/{sample}.hisat2.sorted.rmdup.realigned.bam"
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
		"realigned_bams/{sample}.hisat2.sorted.rmdup.realigned.bam"
	params:
		ref=config["reference"]["V1_01"]["custom"]
	log:
		"logs/realigner/{sample}.realigner.hisat2.log"
	benchmark:
		"benchmarks/{sample}.intervals.benchmark.txt"
	shell:
		"module load gatk; java -jar -Xmx4g /apps/gatk/3.7.0/GenomeAnalysisTK.jar -T IndelRealigner --filter_reads_with_N_cigar -R {params.ref} -I {input.bam} -targetIntervals {input.interval} -o {output}"

rule samtools_index_realigned:
	input:
		"realigned_bams/{sample}.hisat2.sorted.rmdup.realigned.bam"
	output:
		"realigned_bams/{sample}.hisat2.sorted.rmdup.realigned.bai"
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
		"realigned_bams/{sample}.hisat2.sorted.rmdup.realigned.bai"
	log:
		"logs/samtools_index_realigned/{sample}.log"
	benchmark:
		"benchmarks/{sample}.index_realigned.benchmark.txt"
	shell:
		"module load samtools; samtool index {input}"

rule merge_lane_bams:
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

rule index_merged_bams:
	input:
		"merged_lane_bams/{sample}.hisat2.sorted.rmdup.merged.bam"
	output:
		temp("merged_lane_bams/{sample}.hisat2.sorted.rmdup.merged.bam.bai")
	shell:
		"module load samtools; samtools index {input}"

rule index_merged_bams_hisat2:
	input:
		"merged_lane_bams/{sample}.hisat2.sorted.rmdup.merged.bam"
	output:
		temp("merged_lane_bams/{sample}.hisat2.sorted.rmdup.merged.bam.bai")
	shell:
		"module load samtools; samtools index {input}"

rule samtools_index_replaced:
	input:
		"RG_replaced_bams/{sample}.hisat2.sorted.rmdup.merged.bam"
	output:
		"RG_replaced_bams/{sample}.hisat2.sorted.rmdup.merged.bam.bai"
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
		"RG_replaced_bams/{sample}.hisat2.sorted.rmdup.merged.bam.bai"
	log:
		"logs/samtools_index_realigned/{sample}.log"
	benchmark:
		"benchmarks/{sample}.index_realigned.benchmark.txt"
	shell:
		"module load samtools; samtools index {input}"

rule ReplaceRG_merged:
	input:
		"merged_lane_bams/{sample}.hisat2.sorted.rmdup.merged.bam"
	output:
		"RG_replaced_bams/{sample}.hisat2.sorted.rmdup.merged.bam"
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
		"RG_replaced_bams/{sample}.hisat2.sorted.rmdup.merged.bam"
	params:
		RG_fields="RGID={sample} RGLB={sample} RGPL=illumina RGPU={sample} RGSM={sample}"
	log:
		"logs/picard_replaceRG/{sample}.log"
	shell:
		"module load picard; java -jar $HPC_PICARD_DIR/picard.jar AddOrReplaceReadGroups I={input} O={output} {params.RG_fields} &> {log}"

rule samtools_index_merged:
        input:
                "merged_lane_bams/{sample}.hisat2.sorted.rmdup.merged.bam"
        output:
                "merged_lane_bams/{sample}.hisat2.sorted.rmdup.merged.bam.bai"
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
                "merged_lane_bams/{sample}.hisat2.sorted.rmdup.merged.bam.bai"
        log:
                "logs/samtools_index_merged/{sample}.log"
        benchmark:
                "benchmarks/{sample}.index_merged.benchmark.txt"
        shell:
                "module load samtools; samtools index {input}"

rule mosdepth_hisat2:
	input:
		bam="RG_replaced_bams/{sample}.hisat2.sorted.rmdup.merged.bam",
		bai="RG_replaced_bams/{sample}.hisat2.sorted.rmdup.merged.bam.bai"
	params:
		options=" --by /home/d.ence/projects/pinus_taeda_L/Fr1_project/probe_files/RG_0610_merged.fixed.bed --no-per-base ./mosdepth/{sample}",
		gz="mosdepth/hisat2/{sample}.regions.bed.gz"	
	output:
		temp(bed="mosdepth/hisat2/{sample}.regions.bed")	
	log:
		"logs/mosdepth/hisat2/{sample}.mosdepth.log"
	shell:
		"module load mosdepth; mosdepth {params.options} {input.bam} ; gunzip {params.gz}"	

rule mosdepth_hisat2:
	input:
		bam="RG_replaced_bams/{sample}.hisat2.sorted.rmdup.merged.bam",
		bai="RG_replaced_bams/{sample}.hisat2.sorted.rmdup.merged.bam.bai"
	params:
		options=" --by /home/d.ence/projects/pinus_taeda_L/Fr1_project/probe_files/RG_0610_merged.fixed.bed --no-per-base ./mosdepth/{sample}",
		gz="mosdepth/hisat2/{sample}.regions.bed.gz"	
	output:
		temp(bed="mosdepth/hisat2/{sample}.regions.bed")	
	log:
		"logs/mosdepth/hisat2/{sample}.mosdepth.log"
	shell:
		"module load mosdepth; mosdepth {params.options} {input.bam} ; gunzip {params.gz}"	
		
rule compile_non_10_5_indv_mosdepth_report:
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

rule compile_non_10_5_pooled_mosdepth_report:
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

rule compile_10_5_megs_report:
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

rule compile_10_5_prog_report:
	input:
		prog_mosdepths=expand("mosdepth/hisat2/{sample}.regions.bed",sample=config["Fr1_prog_samples"]),
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
		prog_mosdepths=expand("mosdepth/hisat2/{sample}.regions.bed",sample=config["Fr1_prog_samples"]),
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
rule report:
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

		This resulted in {n_calls} variants (see Table T1_).
		Benchmark results for BWA can be found in the tables.
		""", output[0], T1=input)



