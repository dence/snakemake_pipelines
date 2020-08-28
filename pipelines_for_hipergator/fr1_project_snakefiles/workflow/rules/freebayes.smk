def make_bam_list(prefix,bam_list_obj):
    #print("in make_bam_list")
    #print(bam_list_obj)
    #print every bam from bam_list_obj to a list in the results dir
    #use prefix to name the list files
    #return the list filename
    list_filename="results/" + prefix + ".list.txt"
    list_file = open(list_filename,'w')

    for bam_file in bam_list_obj:
        list_file.write(bam_file)

    list_file.close()
    return list_filename

rule freebayes_haploid:
    input:
        bam=expand("results/realigned/{sample}.realigned.bam", sample=get_sample_subset("10_5_megagametophyte")),
        bai=expand("results/realigned/{sample}.realigned.bam.bai", sample=get_sample_subset("10_5_megagametophyte")),
        bam_list=make_bam_list("10_5_megagametophyte", expand("results/realigned/{sample}.realigned.bam", sample=get_sample_subset("10_5_megagametophyte"))),
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
