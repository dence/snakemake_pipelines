rule deconvolute_haploid_vcf:
    input:
        "results/calls/freebayes/{call_group}.freebayes.vcf",
        get_reference
    output:
        "results/calls/postprocessed/{call_group}.freebayes.postprocess.vcf"
    shell:
        "module load bedtools; module load vcflib; module load vt; cat {input[0]} | vcfallelicprimitives --keep-info --keep-geno | vt normalize -r {input[1]} - > {output}"
