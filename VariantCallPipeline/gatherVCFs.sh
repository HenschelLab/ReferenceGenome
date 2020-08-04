#!/bin/sh

module load picard
module load java

java -Xmx32G -jar /research/btc_bioinformatic/operations/picard/build/libs/picard.jar  GatherVcfs \
     I= /research/btc_bioinformatic/results/GenotypeGVCFs/genotypeVcfs.list \
     O= /research/btc_bioinformatic/results/gathered_vcfs.vcf.gz 
