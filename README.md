# Results
alignment statistics coverage
Joint genotyping
29M 


# Methods

Reference Genome construction for the UAE

Raw data provided in some raw data dir

pipeline to construct single sample (g)vcf:



Variant calling
Trimmomatic
Pipeline - HPC
Quality control

GVCF - advantages
Combined VCF


exonic regions mapping file - pipeline_WES.py
same coordinate system thus combinable

# Filtering VQSR 
genomes + exomes integrate 

Use GenomicsDB, Joint genotype calling, requires parallelization along genome regions due 
to its computational expenses - split by 10Mbp


#Annotation 
variant stats
novel variants

#UAE reference construction and testing
GATK
Joint Genotype calling

### Variant reduction (reference hg19 vs UAE):
We compare the number of called variants with respect to two reference genome: 1. hg19 and 2. our own (UAERG).
The call stats are produced with a simple shell (bash) script, `gatherNrVariantsInVCFs.sh`.

Lists the number of variants for selected genomes, once for reference 
`run concatVariantStats.py`
which generates a single summary spreadsheet (variantResults.csv)

`         12723R    12723U    12753R    12753U    12801R    12801U    12811R    12811U   12723D  12723Dpct   12753D  12753Dpct 
chr1   402961.0  323321.0  405401.0  321405.0  392219.0  308402.0  402864.0  320065.0  79640.0  19.763699  83996.0  20.719238 
chr2   412000.0  324350.0  403261.0  318001.0  399028.0  313695.0  211792.0  329532.0  87650.0  21.274272  85260.0  21.142635 
chr3   333239.0  260852.0  343642.0  277707.0  325838.0  259316.0       NaN  280774.0  72387.0  21.722247  65935.0  19.187119 
chr4   358300.0  281472.0  364340.0  287789.0  354895.0  278450.0       NaN  306839.0  76828.0  21.442367  76551.0  21.010869 
chr5   295294.0  241912.0  300598.0  242585.0  293180.0  236954.0       NaN  232489.0  53382.0  18.077577  58013.0  19.299197 `

