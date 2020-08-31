#Daniel Ence
#August 31, 2020

rule gall_vs_nogall_rrBLUP_GWAS:
    input:
        data_012="results/vcf_012/gall_vs_nogall.012",
        pos_012="results/vcf_012/gall_vs_nogall.012.pos",
        indv="results/vcf_012/gall_vs_nogall.012.indv"
    output:
        manhattanplot="results/GWAS/rrBLUP/rrBLUP.manhattanplot.pdf",
        kinship="results/GWAS/rrBLUP/kinship.rrblup.pdf",
        qqplot="results/GWAS/rrBLUP/qqplot.pdf"
    params:
        data_plate_layout="config/UFL_104301_Plate_Layout_Dec2016.txt",
        data_plate_extractions="config/DNA_plate_layout.txt"
    log:
        "logs/R/rrBLUP_GWAS_gall_vs_nogall.log"
    script:
        "../scripts/gall_vs_nogall_GWAS.R"
