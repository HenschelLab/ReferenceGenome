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

## Pipeline

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

`         12723R    12723U    12753R    12753U    12801R    12801U    12811R    12811U   12723D  12723Dpct   12753D  12753Dpct 
chr1   402961.0  323321.0  405401.0  321405.0  392219.0  308402.0  402864.0  320065.0  79640.0  19.763699  83996.0  20.719238 
chr2   412000.0  324350.0  403261.0  318001.0  399028.0  313695.0  211792.0  329532.0  87650.0  21.274272  85260.0  21.142635 
chr3   333239.0  260852.0  343642.0  277707.0  325838.0  259316.0       NaN  280774.0  72387.0  21.722247  65935.0  19.187119 
chr4   358300.0  281472.0  364340.0  287789.0  354895.0  278450.0       NaN  306839.0  76828.0  21.442367  76551.0  21.010869 
chr5   295294.0  241912.0  300598.0  242585.0  293180.0  236954.0       NaN  232489.0  53382.0  18.077577  58013.0  19.299197 `

