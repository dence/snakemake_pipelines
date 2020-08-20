#Daniel Ence, 08/20/2020
#based heavily (probably wholly) on common rule from
#https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/blob/master/workflow/rules/common.smk

from snakemake.utils import validate
import pandas as pd

configfile: "config/config.yaml"
validate(config, schema="../schemas/confi.schema.yaml")

samples = pd.read_csv(config["samples"], sep="\t", dtype=str).set_index("sample", drop=False)
samples.index.names = ["sample_id"]
validate(samples, schema="../schemas/samples.schema.yaml")

#my samples spread across two lanes,
#so need to account for this in my workflow
units = pd.read_csv(
    config["units"],dtype=str,sep="\t").set_index(["sample"], "unit"], drop=False)
units.index.names = ["sample_id", "unit_id"]
units.index = units.index.set_levels(
    [i.astype(str) for i in units.index.levels])
validate(units, schema="../schemas/units.schema.yaml")

report: "../report/workflow.rst"

##### wildcard constraints ####

wildcard_constraints:
    sample="|".join(samples.index),
    units="|".join(units["units"])


###### helpers #######

def is_single_end(sample, unit):
        """Determine whether unit is single-ended."""
        fq_present = pd.isnull(units.loc[(sample, unit), "fq2"])
        if isinstance(fq2_present, pd.core.series.Series):

            #if this is the case, get_fastqs cannot work properly
            raise ValuleError(
                f"Multiple fq2 entries found for sample-unit combination {sample}-{unit}.\n"
                "This is most likely due to a faulty units.tsv file, e.g. "
                "a unit name is used twice for the same sample.\n"
                "Try checking your units.tsv for duplicates."
            )
        return fq2_present

def get_fastqs(wildcards):
    """Get raw FASTQ files from unit sheet."""
    if is_single_end(wildcards.sample, wildcards.unit):
        return units.loc[ (wildcards.sample, wildcards.unit), "fq1")]
    else:
        u = units.loc[ (wildcards.sample, wildcards.unit), ["fq1","fq2"] ].dropna()
        return [ f"{u.fq1}", f"{u.fq2}" ]

def get_trimmed(wildcards):
    if not is_single_end(**wildcards):
        #paired-end sample
        return expand("results/trimmed/{sample}-{unit}.{group}.fastq.gz",
                        group=[1,2], **wildcards)
    #single end sample
    return expand("results/trimmed/{sample}-{unit}.{group}.fastq.gz", **wildcards)

#I don't know if I'll need these rules, so comment out for now
#def get_bioc_species_pkg(wildcards)
#   """Get the package bioconductor package name for the species in config.yaml"""
#   species_letters = config["resources"]["ref"]["species"][0:2].capitalize()

#def get_bioc_pkg_path(wildcards):
#   return "resources/bioconductor/lib/R/library/{pkg}.format(pkg=get_bioc_species_pkg(wildcards))"

#def is_activated(config_element):
#   return config_element['activate'] in {"true","True"}
