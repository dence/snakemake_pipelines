
rule gatk_haplotypecaller_haploid:
	input:
		bam="results/realigned/{haploid_samples}.realigned.bam",
		bai="results/realigned/{haploid_samples}.realigned.bam.bai",
		ref=get_reference
	params:
        #targets="--targets " + config["resources"]["intervals"],
		settings="-ploidy 1 "
	log:
		"logs/gatk_haplotypecaller/Fr1_{haploid_samples}.gatk_haplotypecaller.log"
	output:
		temp("results/calls/gatk_haplotypecaller/Fr1_{haploid_samples}.haplotypecaller.g.vcf.gz")
	shell:
        #"unset TMPDIR; module load freebayes/1.3.2; freebayes {params.settings} {params.targets}  --bam-list {input.bam_list} -f {input.ref} --vcf {output}  &> {log}"
		"unset TMPDIR; module load gatk/3.7.0; java -jar -Xmx9g /apps/gatk/3.7.0/GenomeAnalysisTK.jar -T HaplotypeCaller  -R {input.ref} -I {input.bam} -o {output} -ERC GVCF {params.settings} &> {log}"

rule gatk_CombineGVCFS_haploid:
	input:
		gvcfs=get_gvcf_list(["10_5_megagametophyte"]),
		ref=get_reference
	log:
		"logs/gatk_combinegvcfs/haploid_samples.log"
	output:
		temp("results/calls/gatk_combinegvcfs/Fr1_megs.haplotypecaller.combined.g.vcf.gz")
	shell: """
			unset TMPDIR; module load gatk/3.7.0;
			GVCFS=$(echo {input.gvcfs});
			echo {{$GVCFS}} &>> {log}
			VARIANTS=""
			for file in ${{GVCFS}}
			do
				VARIANTS="--variant ${{file}} ${{VARIANTS}}"
				echo ${{VARIANTS}} &>> {log}
			done
			echo "checking VARIANTS variable"
			echo ${{VARIANTS}} &>> {log}
			java -jar -Xmx9g /apps/gatk/3.7.0/GenomeAnalysisTK.jar -T CombineGVCFs  -R {input.ref} -o {output} $VARIANTS &>> {log}
			"""

rule gatk_GenotypeGVCFs_haploid:
	input:
		combined_gvcf="results/calls/gatk_combinegvcfs/Fr1_megs.haplotypecaller.combined.g.vcf.gz",
		ref=get_reference
	log:
		"logs/gatk_genotype_gvcfs/haploid_samples.log"
	output:
		"results/calls/gatk_genotype_gvcfs/Fr1_megs.haplotypecaller.genotyped.vcf.gz"
	shell:
		"unset TMPDIR; module load gatk/3.7.0; java -jar -Xmx9g /apps/gatk/3.7.0/GenomeAnalysisTK.jar -T GenotypeGVCFs  -R {input.ref} -V {input.combined_gvcf} -o {output} &> {log}"

rule gatk_haplotypecaller_diploid:
	input:
		bam="results/realigned/{diploid_samples}.realigned.bam",
		bai="results/realigned/{diploid_samples}.realigned.bam.bai",
		ref=get_reference
	params:
        #targets="--targets " + config["resources"]["intervals"],
		settings="-ploidy 2 "
	log:
		"logs/gatk_haplotypecaller/Fr1_{diploid_samples}.gatk_haplotypecaller.log"
	output:
		temp("results/calls/gatk_haplotypecaller/Fr1_{diploid_samples}.haplotypecaller.g.vcf.gz")
	shell:
        #"unset TMPDIR; module load freebayes/1.3.2; freebayes {params.settings} {params.targets}  --bam-list {input.bam_list} -f {input.ref} --vcf {output}  &> {log}"
		"unset TMPDIR; module load gatk/3.7.0; java -jar -Xmx9g /apps/gatk/3.7.0/GenomeAnalysisTK.jar -T HaplotypeCaller  -R {input.ref} -I {input.bam} -o {output} -ERC GVCF {params.settings} &> {log}"

rule gatk_CombineGVCFS_diploid:
	input:
		gvcfs=get_gvcf_list(["nongalled_10_5_prog","galled_10_5_prog","20_1010","elite_pine_family","10_5_x_4_6664"]),
		ref=get_reference
	log:
		"logs/gatk_combinegvcfs/diploid_samples.log"
	output:
		temp("results/calls/gatk_combinegvcfs/diploid_samples.haplotypecaller.combined.g.vcf.gz")
	shell: """
		unset TMPDIR; module load gatk/3.7.0;
		GVCFS=$(echo {input.gvcfs});
		echo {{$GVCFS}} &>> {log}
		VARIANTS=""
		for file in ${{GVCFS}}
		do
			VARIANTS="--variant ${{file}} ${{VARIANTS}}"
			echo ${{VARIANTS}} &>> {log}
		done
		echo "checking VARIANTS variable"
		echo ${{VARIANTS}} &>> {log}
		java -jar -Xmx9g /apps/gatk/3.7.0/GenomeAnalysisTK.jar -T CombineGVCFs  -R {input.ref} -o {output} $VARIANTS &>> {log}
		"""

rule gatk_GenotypeGVCFs_diploid:
	input:
		combined_gvcf="results/calls/gatk_combinegvcfs/diploid_samples.haplotypecaller.combined.g.vcf.gz",
		ref=get_reference
	log:
		"logs/gatk_genotype_gvcfs/diploid_samples.log"
	output:
		"results/calls/gatk_genotype_gvcfs/diploid_samples.haplotypecaller.genotyped.vcf.gz"
	shell:
		"unset TMPDIR; module load gatk/3.7.0; java -jar -Xmx9g /apps/gatk/3.7.0/GenomeAnalysisTK.jar -T GenotypeGVCFs  -R {input.ref} -V {input.combined_gvcf} -o {output} &> {log}"
