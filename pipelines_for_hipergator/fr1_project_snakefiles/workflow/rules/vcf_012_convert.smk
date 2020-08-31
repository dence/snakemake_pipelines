rule vcf_012_convert:
#read in a vcf file
#output a 012 file
    input:
        "results/calls/freebayes/Fr1_diploid.freebayes.vcf"
    output:
        "results/vcf_012/gall_vs_nogall.012",
        "results/vcf_012/gall_vs_nogall.012.pos",
        "results/vcf_012/gall_vs_nogall.012.indv"
    params:
        "results/vcf_012/gall_vs_nogall",
        make_sample_type_list_file(["nongalled_10_5_prog","galled_10_5_prog"])
    log:
    shell:
        "module load vcftools; vcftools --vcf {input[0]} --012 --keep {params[1]} --out {params[0]}"
