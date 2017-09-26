import re
import argparse

parser = argparse.ArgumentParser(description='replace template fields with specific fastq files')
parser.add_argument("fq1", type=str, help="first fastq file")
parser.add_argument("fq2", type=str, help="second fastq file")
parser.add_argument("template", type=str, help="template config file for masurca")

args = parser.parse_args()

template_file = open(args.template,'r')
for line in template_file:
	line = line.strip()
	if(re.search("FQ1",line)):
		line  = re.sub("FQ1",args.fq1,line)
	
	if(re.search("FQ2",line)):
		line = re.sub("FQ2",args.fq2,line)
	print line
		



