
rule freebayes_haploid:
    input:
        bam=get_bam_list(["10_5_megagametophyte"]),
        bai=get_bai_list(["10_5_megagametophyte"]),
        bam_list=make_bam_list_file("10_5_megagametophyte", get_bam_list(["10_5_megagametophyte"])),
        ref=get_reference
    params:
        targets="--targets " + config["resources"]["intervals"],
        settings="--pvar 0.75 --theta 0.01 --min-alternate-fraction 0.8 --min-alternate-count 1 "
    log:
        "logs/freebayes/Fr1_megs.freebayes.log"
    output:
        "results/calls/freebayes/Fr1_megs.freebayes.vcf"
    threads:
        10
    shell:
        "unset TMPDIR; module load freebayes/1.3.1; freebayes-parallel /blue/kirst/d.ence/ref_genomes/loblolly_pine/V1_1/ref_for_Fr1_probes/ptaeda.v1.01.masked.and_elite_RNAseq.probes_bed_scaffolds.fasta.regions {threads} {params.settings} {params.targets}  --bam-list {input.bam_list} -f {input.ref} --vcf {output}  &> {log}"

rule freebayes_diploid:
    input:
        bam=get_bam_list(["nonqgalled_10_5_prog","galled_10_5_prog","20_1010","elite_pine_family","10_5_x_4_6664"]),
        bai=get_bai_list(["nongalled_10_5_prog","galled_10_5_prog","20_1010","elite_pine_family","10_5_x_4_6664"]),
        bam_list=make_bam_list_file("diploid", get_bam_list(["nongalled_10_5_prog","galled_10_5_prog","20_1010","elite_pine_family","10_5_x_4_6664"])),
        ref=get_reference
    params:
        targets="--targets " + config["resources"]["intervals"],
        settings="--pvar 0.75 --theta 0.01 --min-alternate-fraction 0.8 --min-alternate-count 1 "
    log:
        "logs/freebayes/Fr1_prog.freebayes.log"
    output:
        "results/calls/freebayes/Fr1_diploid.freebayes.vcf"
    threads:
        10
    shell:
        "unset TMPDIR; module load freebayes/1.3.1; freebayes-parallel /blue/kirst/d.ence/ref_genomes/loblolly_pine/V1_1/ref_for_Fr1_probes/ptaeda.v1.01.masked.and_elite_RNAseq.probes_bed_scaffolds.fasta.regions {threads} {params.settings} {params.targets}  --bam-list {input.bam_list} -f {input.ref} --vcf {output}  &> {log}"
