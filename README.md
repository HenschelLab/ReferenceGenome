# UAE reference construction and testing
This repository contains the scripts for data analysis and reference genome construction.
It is a Python and shell script collection.
The variant calling pipeline is assumed to run on a high performance computing (HPC) facility. Ours runs CentOS on all nodes. Parallelization is per genome, i.e the command line argument is some identifier that either explicitly or implicitely serves as sample id.
central classes

## Dependencies:
GATK requires a number of additional support files, including a reference genome (indexed).
We use hg19, and locate it as per configuration in pipeline.py.
We use GATK version 4.0.6.0 and Python 3.6+. In addition, Anaconda environments are used.


## Assumptions
Various assumptions regarding the naming conventions of input data are made. See configuration part in `pipeline.py`.

## Pipeline for Variant Calling

Pipeline summary to construct single sample (g)vcf:
`pipeline.py`

#### Part 1:
Tools for data preprocessing:
* FastQC v0.11.7 tool for quality control.
* Trimmomatic v0.36 tool, to remove low quality and short reads.
* BWA-MEM v0.7.12 for mapping the samples against the reference genome.
* Qualimap v2.2.1 for checking mapping quality and mean coverage per sample.
* Picard v2.9.4 for marking and removing duplicate reads and sorting the resulting BAM file.

#### Part 2:
The joint variant calling workflow is designed to run on our in-house High Performance Computing (HPC), using the following tools:

* Variants were called using the Genome Analysis Toolkit (GATK) v4.0.6.0 GVCF workflow that included
* BQSR
* HaplotypeCaller (-ERC GVCF)

For joint genotyping  on the HPC:
The use of GenomicsDB and Joint genotype calling, requires parallelization along genome regions due to its computational expenses  so we split by 10Mbp 
* GenomicsDBImport: bsub files for region-wise import to genomicsDB are generated, 'genomicsDBscripts.py'
* GenotypeGVCFs: bsub files for region-wise joint gentotyping bsub files for region-wise are generated,  'genotypeGVCFsScripts.py'
* Gather all vcfs genetared from the region split GenotypeGVCFs to one vcf file,  'gatherVCFs.sh'
* VQSR is then applied to the the gathered VCF file 'pipeline_VQSR.py'

#### Part 3:
Variant annotaions
* SnpEFF v4.3t for functional annotaion of the VCF.
* GEMINI v0.30.2 tool was used to annotate each variant by integrating several clinical and functional genome annotations.

For Structural Variants(SV) calling :
* Manta v1.6.0-0 and Delly v0.8.2 joint genotyping germline Structural Variants calling workflows parallelized on our in-house High-Performance Computer (HPC) are used.
For consensus SV call sets from the results of Manta and Delly:
* SURVIVOR v1.0.6 tool was used to merge across SV callers and across individuals and generate a union call set and an intersection call set, for which the Structural Variants frequency was calculated
* The tool AnnotSV v2.1, an integrated tool for structural variations annotation was used to annotate the SV calls.

* Visual representation of the spatial variability of SNVs, and SVs across the UAE genomes has been generated using Circos v-0.69-8

### Methods Details

* Raw data is expected in  rawdata directory as per configuration in pipeline.py
* The pipeline version that deals with Whole Exomes is `pipeline_WES.py`.
* WES requires exonic regions mapping file as per configuration in that script (see global variable TR).
* WES follows same coordinate system as WGS, thus efforts are combinable.


### HPC usage

* We deploy IBM's LSF queuing system, using bsub for job submisssion. 
* bsub scripts are provided in the individual directories and might require adaptation to the specific HPC at hand.

### Variant reduction (reference hg19 vs UAE):
* We compare the number of called variants with respect to two reference genome: 1. hg19 and 2. our own (UAERG).
* The call stats are produced with a simple shell (bash) script, `gatherNrVariantsInVCFs.sh`.

* Lists the number of variants for selected genomes, once for reference, `run concatVariantStats.py`
which generates a single summary spreadsheet (variantResults.csv), what eventually is reported as Table 3 in the manuscript.
     
