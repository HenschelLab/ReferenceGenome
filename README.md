# UAE reference construction and testing
Python and shell script collection.
The pipeline is assumed to run on a high performance computing (HPC) facility. Ours runs CentOS on all nodes. Parallelization is per genome, i.e the command line argument is some identifier that either explicitly or implicitely serves as sample id.
central classes

## Dependencies:
GATK requires a number of additional support files, including a reference genome (indexed).
We use hg19, and locate it as per configuration in pipeline.py.
We use GATK version 4.0.6.0 and Python 3.6+. In addition, Anaconda environments are used.


## Assumptions
Various assumptions regarding the naming conventions of input data are made. See configuration part in `pipeline.py`.

## Pipeline for Variant Calling

Pipeline to construct single sample (g)vcf:
`pipeline.py`
* Variant calling
* QC: FastQC, possibly before and certainly after trimming
* Trimmomatic
* Pipeline - HPC
* Quality control

* GVCF 
* Joint genotype calling with Combined VCF


### Methods Details

Raw data provided is expected in some rawdata dir as per configuration in pipeline.py

The ipeline version that deals with Whole Exomes is `pipeline_WES.py`.
Requires exonic regions mapping file as per configuration in that script (see global variable TR).
WES follows same coordinate system as WGS, thus efforts are combinable.

### Filtering VQSR 
genomes + exomes integrate 

Use GenomicsDB, Joint genotype calling, requires parallelization along genome regions due 
to its computational expenses - split by 10Mbp


### HPC usage
We deploy IBM's LSF queuing system, using bsub for job submisssion. 
bsub scripts are provided in the individual directories and mightt require adaptation to the specific HPC at hand. Note that some bsub scripts


#Annotation 
variant stats
novel variants


# Results
alignment statistics coverage
Joint genotyping
29M variant loci discovered.

### Variant reduction (reference hg19 vs UAE):
We compare the number of called variants with respect to two reference genome: 1. hg19 and 2. our own (UAERG).
The call stats are produced with a simple shell (bash) script, `gatherNrVariantsInVCFs.sh`.

Lists the number of variants for selected genomes, once for reference 
`run concatVariantStats.py`
which generates a single summary spreadsheet (variantResults.csv), what eventually is reported as Table 3 in the manuscript.
`     
