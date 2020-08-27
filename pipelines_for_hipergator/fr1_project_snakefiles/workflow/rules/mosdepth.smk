rule mosdepth_bwa_mem:
    input:
        bam="results/realigned/{sample}.sorted.rmdup.realigned.bam"
        bai="results/realigned/{sample}.sorted.rmdup.realigned.bam.bai"
    params:
        options=" --by config[\"resources\"][\"intervals\"] --no-per-base ./results/mosdepth/{sample}",
        gz="results/mosdepth/{sample}.regions.bed.gz"
    output:
        bed=temp("results/mosdepth/{sample}.regions.bed")
    log:
        "logs/mosdepth//{sample}.mosdepth.log"
    shell:
        "module load mosdepth; mosdepth {params.options} {input.bam} ; gunzip {params.gz}"
