rule ReplaceRG_merged_bwa_mem:
    input:
        "results/merged_lane_bams/{sample}.merged.bam"
    output:
        temp("results/RG_replaced_bams/{sample}.bwa_mem.sorted.rmdup.merged.bam")
    params:
        RG_fields="RGID={sample} RGLB={sample} RGPL=illumina RGPU={sample} RGSM={sample}"
    log:
        "logs/picard_replaceRG/{sample}.log"
    shell:
        "module load picard; java -jar $HPC_PICARD_DIR/picard.jar AddOrReplaceReadGroups I={input} O={output} {params.RG_fields} &> {log}"
