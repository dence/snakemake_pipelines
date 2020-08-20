rule all:
	input:
		"freebayes_calls.report.html"

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
		"module load freebayes; freebayes-v1.3.1 --ploidy 2 --bam-list {input.bam_list} -f {input.ref} --vcf {output} --targets /home/d.ence/projects/pinus_taeda_L/Fr1_project/probe_files/RG_0610_merged.fixed.bed &> {log}"

#rule freebayes_progeny_hisat2:
#	input:
#		bam=expand("RG_replaced_bams/{sample}.hisat2.sorted.rmdup.merged.bam",sample=config["Fr1_prog_samples"]),
#		bai=expand("RG_replaced_bams/{sample}.hisat2.sorted.rmdup.merged.bam.bai",sample=config["Fr1_prog_samples"]),
#		bam_list="/home/d.ence/projects/pinus_taeda_L/Fr1_project/aligning_Fr1_samples/aligning_to_custom_R_gene_V1_1_ref/calls/bam_list.Fr1_prog.hisat2.txt",
#		ref=config["reference"]["V1_01"]["custom"]
#	log:
#		"logs/freebayes/Fr1_prog.freebayes.log"
#	output:
#		"calls/freebayes/Fr1_prog.hisat2.freebayes.vcf"
#	shell:
#		"module load freebayes; freebayes-v1.3.1 --ploidy 2 --bam-list {input.bam_list} -f {input.ref} --vcf {output} --targets /home/d.ence/projects/pinus_taeda_L/Fr1_project/probe_files/RG_0610_merged.fixed.bed &> {log}"

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
#		bam=expand("RG_replaced_bams/{sample}.bwa_mem.sorted.rmdup.merged.bam",sample=config["Fr1_prog_samples"]),
#		bai=expand("RG_replaced_bams/{sample}.bwa_mem.sorted.rmdup.merged.bam.bai",sample=config["Fr1_prog_samples"]),
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
		"module load bedtools vcflib vt; cat {input.megs_vcf} | vcfallelicprimitives --keep-info --keep-geno | vt normalize -r {input.ref} - > {output.megs_deconv_vcf}"

rule postprocess_megs_hisat2:
	input:
		megs_vcf="calls/freebayes/Fr1_megs.hisat2.freebayes.vcf",
		ref=config["reference"]["V1_01"]["custom"]
	output:
		megs_deconv_vcf="calls/freebayes/Fr1_megs.hisat2.freebayes.deconvoluted.vcf"
	shell:
		"module load bedtools vcflib vt; cat {input.megs_vcf} | vcfallelicprimitives --keep-info --keep-geno | vt normalize -r {input.ref} - > {output.megs_deconv_vcf}"

rule postprocess_progs_hisat2:
	input:
		prog_vcf="calls/freebayes/Fr1_prog.hisat2.freebayes.vcf",
		ref=config["reference"]["V1_01"]["custom"]
	output:
		prog_deconv_vcf="calls/freebayes/Fr1_prog.hisat2.freebayes.deconvoluted.vcf"
	shell:
		"module load bedtools vcflib vt; cat {input.prog_vcf} | vcfallelicprimitives --keep-info --keep-geno | vt normalize -r {input.ref} - > {output.prog_deconv_vcf};" 

rule postprocess_progs:
	input:
		prog_vcf="calls/freebayes/Fr1_prog.bwa_mem.freebayes.vcf",
		ref=config["reference"]["V1_01"]["custom"]
	output:
		prog_deconv_vcf="calls/freebayes/Fr1_prog.bwa_mem.freebayes.deconvoluted.vcf"
	shell:
		"module load bedtools vcflib vt; cat {input.prog_vcf} | vcfallelicprimitives --keep-info --keep-geno | vt normalize -r {input.ref} - > {output.prog_deconv_vcf};" 

rule postprocess_other_vcfs:
	input:
		non_10_5_pooled_vcf="calls/freebayes/non_10_5_pooled.bwa_mem.freebayes.vcf",
		non_10_5_indv_vcf="calls/freebayes/non_10_5_indv.bwa_mem.freebayes.vcf",
		ref=config["reference"]["V1_01"]["custom"]
	output:
		non_10_5_pooled_deconv_vcf="calls/freebayes/non_10_5_pooled.bwa_mem.freebayes.deconvoluted.vcf",
		non_10_5_indv_deconv_vcf="calls/freebayes/non_10_5_indv.bwa_mem.freebayes.deconvoluted.vcf"
	shell:
		"module load bedtools vcflib vt; cat {input.non_10_5_pooled_vcf} | vcfallelicprimitives --keep-info --keep-geno | vt normalize -r {input.ref} - > {output.non_10_5_pooled_deconv_vcf}; cat {input.non_10_5_indv_vcf} | vcfallelicprimitives --keep-info --keep-geno | vt normalize -r {input.ref} - > {output.non_10_5_indv_deconv_vcf}"
		
rule postprocess_other_vcfs_hisat2:
	input:
		non_10_5_pooled_vcf="calls/freebayes/non_10_5_pooled.hisat2.freebayes.vcf",
		non_10_5_indv_vcf="calls/freebayes/non_10_5_indv.hisat2.freebayes.vcf",
		ref=config["reference"]["V1_01"]["custom"]
	output:
		non_10_5_pooled_deconv_vcf="calls/freebayes/non_10_5_pooled.hisat2.freebayes.deconvoluted.vcf",
		non_10_5_indv_deconv_vcf="calls/freebayes/non_10_5_indv.hisat2.freebayes.deconvoluted.vcf"
	shell:
		"module load bedtools vcflib vt; cat {input.non_10_5_pooled_vcf} | vcfallelicprimitives --keep-info --keep-geno | vt normalize -r {input.ref} - > {output.non_10_5_pooled_deconv_vcf}; cat {input.non_10_5_indv_vcf} | vcfallelicprimitives --keep-info --keep-geno | vt normalize -r {input.ref} - > {output.non_10_5_indv_deconv_vcf}"

rule merge_freebayes_non_10_5_pooled_singles:
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

rule merge_freebayes_Fr1_prog:
	input:
		bgzipped_vcf=expand("calls/freebayes/individual/Fr1.{sample}.bwa_mem.freebayes.vcf.gz",sample=config["Fr1_prog"])
	output:
		vcf="calls/freebayes/Fr1_prog.bwa_mem.freebayes.vcf"
	log:
		"logs/freebayes/Fr1_prog.bwa_mem.freebayes.log"
	shell:
		"module load vcftools; vcf-merge -d --ref-for-missing 0/0 {input.bgzipped_vcf} > {output.vcf} 2> {log}"

rule merge_freebayes_Fr1_prog_hisat2:
	input:
		bgzipped_vcf=expand("calls/freebayes/individual/Fr1.{sample}.hisat2.freebayes.vcf.gz",sample=config["Fr1_prog"])
	output:
		vcf="calls/freebayes/Fr1_prog.hisat2.freebayes.vcf"
	log:
		"logs/freebayes/Fr1_prog.hisat2.freebayes.log"
	shell:
		"module load vcftools; vcf-merge -d --ref-for-missing 0/0 {input.bgzipped_vcf} > {output.vcf} 2> {log}"

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

rule report:
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



