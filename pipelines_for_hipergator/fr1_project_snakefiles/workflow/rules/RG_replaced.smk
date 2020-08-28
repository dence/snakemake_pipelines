rule ReplaceRG_merged:
    input:
        "results/merged_lane_bams/{sample}.merged.bam"
    output:
        "results/RG_replaced_bams/{sample}.merged.bam"
    params:
        RG_fields="RGID={sample} RGLB={sample} RGPL=illumina RGPU={sample} RGSM={sample}"
    log:
        "logs/picard_replaceRG/{sample}.log"
    shell:
        "unset TMPDIR; module load picard; java -jar $HPC_PICARD_DIR/picard.jar AddOrReplaceReadGroups I={input} O={output} {params.RG_fields} &> {log}"
