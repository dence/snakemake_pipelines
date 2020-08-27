rule merge_lane_bams:
	input:
        lambda w: expand("results/mapped/{sample}-{unit}",
            sample = w.sample,
            unit = units.loc[units['sample'] == w.sample].unit.to_list()
        )
		#making a dumb assumption about the names of the bams to merged. specific to the Fr1 project. DE
        #L4_bam="rmduped_reads/{sample}_{unit}.bam",
		#L5_bam="rmduped_reads/{sample}_{unit}.bam"
	output:
		temp("results/merged_lane_bams/{sample}.merged.bam")
	log:
		"logs/picard_merge_sam_files/{sample}.log"
	shell:
		"module load picard; java -jar $HPC_PICARD_DIR/picard.jar MergeSamFiles I={input[0]} I={input[1]} O={output} &> {log}"
