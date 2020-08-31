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
        "unset TMPDIR; module load mosdepth; mosdepth {params.options} {input.bam} ; gunzip {params.gz}"

rule compile_mosdepth_report_bwa_mem:
	input:
		mosdepths="results/mosdepth/{sample}.regions.bed"
	output:
		filename=expand("mosdepth/mosdepth{aligner}.report.txt",aligner=config["settings"]["aligner"])
	log:
	   expand("logs/mosdepth/mosdepth.{aligner}.report.log",aligner=config["settings"]["aligner"])
	run:
        header_line = "sample_ID\tfamily\t"
        header_line_finished=0
        sample_lines = []

        for l in open_list:
            short_sample_ID = l.split('/')[10].split('_i')[0]
            long_sample_ID = l.split('/')[10].split('.')[0]

            curr_mosdepth = open("mosdepth/bwa_mem/" + long_sample_ID + ".regions.bed")
            tmp = curr_mosdepth.readline() #skip the header line

            family = ""
            if(short_sample_ID == "P02_WG07"):
                family = "20-1010"
            else:
                family = "10-5_X_4-6664"
            new_line = short_sample_ID + "\t" + family + "\t"
            for depth_line in curr_mosdepth:
                region = depth_line.strip().split("\t")[0]
                print("region is:\t" + region)
                mean_covg = depth_line.strip().split("\t")[3]
                print("mean_covg:\t" + mean_covg)

                new_line = new_line + mean_covg + "\t"
                if(header_line_finished == 0):
                    header_line = header_line + region + "_mean_covg\t"
            header_line_finished=1
            sample_lines.append(new_line)
            curr_mosdepth.close()

        all_sample = open(output.filename,'w')
        all_sample.write(header_line + "\n")
        for line in sample_lines:
            all_sample.write(line + "\n")
        all_sample.close()
