#!/bin/bash

#BSUB -n 12
#BSUB -q general
#BSUB -W 48:00

#BSUB -o JobOutputs/annotate.%I.%J.out
#BSUB -e JobOutputs/annotate.%I.%J.err

module load java


##Reference_dir="/research/btc_bioinformatic/operations/data_masdar/"
##Picard_dir="/usr/share/BTSSRV/GenomeBrowse/broadinstitute-picard-b7996a6"
##Gatk_dir="/bmshare/gihan/gatk4/gatk-4.0.6.0"
##RESULTS_dir="/research/btc_bioinformatic/results"
##snpEff="/usr/share/BTSSRV/snpEff"
##snpEff= "/research/btc_bioinformatic/operations/snpEff"



# ##annotation with snpEff
java -Xmx32G -jar /research/btc_bioinformatic/operations/snpEff/snpEff.jar \
   -c /research/btc_bioinformatic/operations/snpEff/snpEff.config \
     -stats /research/btc_bioinformatic/results/gathered_vcfs_recal_snp_indel_snpeff.html \
     -i vcf -o gatk  hg19 /research/btc_bioinformatic/results/VQSR/gathered_vcfs_recal_snp_indel.vcf.gz > /research/btc_bioinformatic/results/gathered_vcfs_recal_snp_indel_snpeff.vcf.gz

## annotation dbsnp

java -Xmx32G -jar /research/btc_bioinformatic/operations/snpEff/SnpSift.jar annotate \
    /research/btc_bioinformatic/operations/AuxData/dbSNP/All_20180423.vcf.gz \
     /research/btc_bioinformatic/results/gathered_vcfs_recal_snp_indel_snpeff.vcf > \
   /research/btc_bioinformatic/results/gathered_vcfs_recal_snp_indel_snpeff_dbsnp.vcf

java -Xmx32G -jar /research/btc_bioinformatic/operations/snpEff/SnpSift.jar annotate \
     /research/btc_bioinformatic/operations/AuxData/ClinVAr/clinvar_20190513.vcf.gz \
     /research/btc_bioinformatic/results/gathered_vcfs_recal_snp_indel_snpeff_dbsnp.vcf.gz > \
     /research/btc_bioinformatic/results/gathered_vcfs_recal_snp_indel_snpeff_dbsnp_clinvar.vcf.gz 
## anotation clinvar
##java -Xmx32G -jar /research/btc_bioinformatic/operations/snpEff/SnpSift.jar \
##     annotate /research/btc_bioinformatic/operations/AuxData/ClinVar/clinvar_20190513.vcf.gz \

##     /research/btc_bioinformatic/results/gathered_vcfs_recal_snp_indel_snpeff_dbsnp.vcf.gz > /research/btc_bioinformatic/results/gathered_vcfs_recal_snp_indel_snpeff_dbsnp_clinvar.vcf.gz
## annotation dbNSFP
##java -Xmx32G -jar /research/btc_bioinformatic/operations/snpEff/SnpSift.jar annotate \
##     /research/btc_bioinformatic/operations/AuxData/dbNSFP3.5_hg19/dbNSFP3.5a_hg19.txt.gz \
##     /research/btc_bioinformatic/results/gathered_vcfs_recal_snp_indel_snpeff_dbsnp_clinvar.vcf.gz > \
##     /research/btc_bioinformatic/results/gathered_vcfs_recal_snp_indel_snpeff_dbsnp_clinvar_dbNSFP.vcf.gz



     
     
