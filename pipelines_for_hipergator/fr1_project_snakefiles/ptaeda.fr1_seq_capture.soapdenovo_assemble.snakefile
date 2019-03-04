rule all:
	input:
		"longest_scaffold_report.txt"

#adapted from https://github.com/eclarke/sisyphus/blob/master/Snakefile
rule soap_config:
	input:
		L4_r1 = "/ufrc/kirst/d.ence/pinus_taeda_L/Fr1_project/aligning_Fr1_samples/aligning_to_V1_1/running_all_files/trimmed_reads/{sample}_L004/{sample}_L004.trimmed.R1.fastq.gz",
		L4_r2 = "/ufrc/kirst/d.ence/pinus_taeda_L/Fr1_project/aligning_Fr1_samples/aligning_to_V1_1/running_all_files/trimmed_reads/{sample}_L004/{sample}_L004.trimmed.R2.fastq.gz",
		L5_r1 = "/ufrc/kirst/d.ence/pinus_taeda_L/Fr1_project/aligning_Fr1_samples/aligning_to_V1_1/running_all_files/trimmed_reads/{sample}_L005/{sample}_L005.trimmed.R1.fastq.gz",
		L5_r2 = "/ufrc/kirst/d.ence/pinus_taeda_L/Fr1_project/aligning_Fr1_samples/aligning_to_V1_1/running_all_files/trimmed_reads/{sample}_L005/{sample}_L005.trimmed.R2.fastq.gz"
	log:
		'logs/soapdenovo/{sample}.config.log'	
	output:
		config=temp('configs/{sample}_soap_config.txt')
	run:
		import os
		tmp_log = open(log[0],'w')
		tmp_log.write(os.getcwd() + "\n")
		tmp_log.close()
		if(not os.path.isdir( os.path.join(os.getcwd(), "configs"))):
			os.mkdir(os.getcwd() + "/" + "configs")
		with open(os.getcwd() + "/" +  output[0],'w') as out:
			out.write((
				"max_rd_len=100\n"
				"[LIB]\n"
				"avg_ins=200\n"
				"reverse_seq=0\nasm_flags=3\nrank=1\n"
				"q1=" + input[0] + "\n"
				"q2=" + input[1] + "\n"
				"q1=" + input[2] + "\n"
				"q2=" + input[3] + "\n"))
					
rule soapdenovo:
	input:
		'configs/{sample}_soap_config.txt'
	output:
		scaffolds='{sample}_assembly/{sample}.scafSeq',
		log_file='logs/soapdenovo/{sample}.log'
	log:
		'logs/soapdenovo/{sample}.log'
	params: 
                prefix = '{sample}_assembly/{sample}'
	shell:
		"module load soapdenovo ; SOAPdenovo all -p 1 -s {input} -K 63 -R -o {params.prefix} >& {log} "

rule report:
	input:
		T1=expand("logs/soapdenovo/{sample}.log",sample=config["samples"])
	output:
		"longest_scaffold_report.txt"
	run:
		report = open(output[0],'w')
		for curr_soap in input.T1:
			with open(curr_soap) as soap_out:
				for l in soap_out:
					if(l.startswith("Longest scaffold")):
						longest_scaf = l
						report.write(soap_out.name + ":\t" + longest_scaf + "\n")
				
