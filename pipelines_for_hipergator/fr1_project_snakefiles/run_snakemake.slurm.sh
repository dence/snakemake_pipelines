#!/bin/sh
#SBATCH --job-name=fr1_snakemake # Job name
#SBATCH --mail-type=ALL              # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=d.ence@mail.ufl.edu  # Where to send mail
#SBATCH --nodes=1                    # Use one node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=1            # Number of CPU cores per task
#SBATCH --mem=2gb                  # Total memory limit
#SBATCH --time=96:00:00              # Time limit hrs:min:sec
#SBATCH --output=fr1_snakemake_%j.out     # Standard output and error log
#SBATCH --qos=peter
#SBATCH --account=peter
#SBATCH --partition=hpg1-compute



#Daniel Ence
#August 26, 2020


unset TMPDIR
module load python3 
mkdir results
snakemake --configfile config/config.yaml --snakefile ./workflow/Snakefile -c 50 --jobs 50 --directory . --cluster-config ../hipergator.cluster.json --cluster "sbatch --qos={cluster.qos} -p {cluster.partition} -c {cluster.c} -n {cluster.N} --mail-type=FAIL --mail-user=d.ence@ufl.edu -t {cluster.time} --mem={cluster.mem} -J "fr1_align" -o fr1_align_%j.out -D /home/d.ence/projects/pinus_taeda_L/Fr1_project/test_pipelines/snakemake_pipelines/pipelines_for_hipergator/fr1_project_snakefiles"
