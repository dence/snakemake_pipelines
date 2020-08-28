import os

def get_bam_list(sample_type_list):
    bam_list = []
    for sample_type in sample_type_list:
        bam_list.extend(expand("results/realigned/{sample}.realigned.bam", /
        sample=get_sample_subset("sample_type")))
    return bam_list

def get_bai_list(sample_type_list):
    bam_list = []
    for sample_type in sample_type_list:
        bam_list.extend(expand("results/realigned/{sample}.realigned.bam.bai", /
        sample=get_sample_subset("sample_type")))
    return bam_list



def make_bam_list_file(prefix,bam_list_obj):
    #print("in make_bam_list")
    #print(bam_list_obj)
    #print every bam from bam_list_obj to a list in the results dir
    #use prefix to name the list files
    #return the list filename

    #need to add a line to check for "results" dir and
    #make the dir if not there
    if not os.path.exists("./results/"):
        os.makedirs("./results/")

    list_filename="results/" + prefix + ".list.txt"
    list_file = open(list_filename,'w')

    for bam_file in bam_list_obj:
        list_file.write(bam_file)

    list_file.close()
    return list_filename

rule freebayes_haploid:
    input:
        bam=get_bam_list(["10_5_megagametophyte"]),
        bai=get_bai_list(["10_5_megagametophyte"]),
        bam_list=make_bam_list_file("10_5_megagametophyte", get_bam_list(["10_5_megagametophyte"])),
        ref=get_reference

    params:
        targets="--targets " + config["resources"]["intervals"],
        settings="-pvar 0.75 -theta 0.01 -indels -mnps -min-alternate-fraction 0.8 -min-alternate-count 1 "
    log:
        "logs/freebayes/Fr1_megs.freebayes.log"
    output:
        "results/calls/freebayes/Fr1_megs.freebayes.vcf"
    shell:
        "unset TMPDIR; module load freebayes/1.3.1; freebayes-v1.3.1. --ploidy 1 {params.settings} {params.targets}  --bam-list {input.bam_list} -f {input.ref} --vcf {output}  &> {log}"

rule freebayes_diploid:
	input:
        bam=get_bam_list(["nongalled_10_5_prog","galled_10_5_prog","20_1010","elite_pine_family","10_5_x_4_6664"]),
        bai=get_bai_list(["nongalled_10_5_prog","galled_10_5_prog","20_1010","elite_pine_family","10_5_x_4_6664"]),
        bam_list=make_bam_list_file("10_5_megagametophyte", get_bam_list(["10_5_megagametophyte"])),
        ref=get_reference

        bam_list=make_bam_list("diploid", /
        get_bam_list(["nongalled_10_5_prog","galled_10_5_prog","20_1010","elite_pine_family","10_5_x_4_6664"])),
        ref=get_reference
	log:
		"logs/freebayes/Fr1_prog.freebayes.log"
	output:
		"calls/freebayes/Fr1_prog.hisat2.freebayes.vcf"
	shell:
		#"module load freebayes; freebayes-v1.3.1 --ploidy 2 --bam-list {input.bam_list} -f {input.ref} --vcf {output} --targets /home/d.ence/projects/pinus_taeda_L/Fr1_project/probe_files/RG_0610_merged.fixed.bed &> {log}"
		#"module load freebayes; freebayes-v1.3.1 --ploidy 2 --theta 0.01 --pvar 0.75 -indels -mnps --min-alternate-fraction 0.8 --min-alternate-count 1 --bam-list {input.bam_list} -f {input.ref} --vcf {output} --targets /home/d.ence/projects/pinus_taeda_L/Fr1_project/probe_files/RG_0610_merged.fixed.bed &> {log}"
		"module load freebayes; freebayes-v1.3.1 --ploidy 2 --theta 0.01 --pvar 0.75 --bam-list {input.bam_list} -f {input.ref} --vcf {output} --targets /home/d.ence/projects/pinus_taeda_L/Fr1_project/probe_files/RG_0610_merged.fixed.bed &> {log}"
