#BSUB -n 12
#BSUB -q general
#BSUB -W 12:00

#BSUB -o BackupJobOutputs/combineGVCFs.%J.out
#BSUB -e BackupJobOutputs/combineGVCFs.%J.err

module load gatk

time gatk --java-options "-Xmx32G" CombineGVCFs\
     -R /research/gutsybugs/applications/RD/data_masdar/ucsc.hg19.fasta \
     --variant /scratch/GenomePipe/120067/GATK/120067_HaplotypeCaller.gvcf.gz\
     --variant /scratch/GenomePipe/120165/GATK/120165_HaplotypeCaller.gvcf.gz\
     --variant /scratch/GenomePipe/13024/GATK/13024_HaplotypeCaller.gvcf.gz\
     --variant /scratch/GenomePipe/120060/GATK/120060_HaplotypeCaller.gvcf.gz\
     --variant /scratch/GenomePipe/120129/GATK/120129_HaplotypeCaller.gvcf.gz\
     --variant /scratch/GenomePipe/120287/GATK/120287_HaplotypeCaller.gvcf.gz\
     --variant /scratch/GenomePipe/120136/GATK/120136_HaplotypeCaller.gvcf.gz\
     --variant /scratch/GenomePipe/120069/GATK/120069_HaplotypeCaller.gvcf.gz\
     --variant /scratch/GenomePipe/120195/GATK/120195_HaplotypeCaller.gvcf.gz\
     --variant /scratch/GenomePipe/12964/GATK/12964_HaplotypeCaller.gvcf.gz\
     --variant /scratch/GenomePipe/120131/GATK/120131_HaplotypeCaller.gvcf.gz\
     --variant /scratch/GenomePipe/120164/GATK/120164_HaplotypeCaller.gvcf.gz\
     --variant /scratch/GenomePipe/10885/GATK/10885_HaplotypeCaller.gvcf.gz\
     --variant /scratch/GenomePipe/13116/GATK/13116_HaplotypeCaller.gvcf.gz\
     --variant /scratch/GenomePipe/120133/GATK/120133_HaplotypeCaller.gvcf.gz\
     --variant /scratch/GenomePipe/120190/GATK/120190_HaplotypeCaller.gvcf.gz\
     --variant /scratch/GenomePipe/120075/GATK/120075_HaplotypeCaller.gvcf.gz\
     --variant /scratch/GenomePipe/13049/GATK/13049_HaplotypeCaller.gvcf.gz\
     --variant /scratch/GenomePipe/12960/GATK/12960_HaplotypeCaller.gvcf.gz\
     --variant /scratch/GenomePipe/13119/GATK/13119_HaplotypeCaller.gvcf.gz\
     --variant /scratch/GenomePipe/12657/GATK/12657_HaplotypeCaller.gvcf.gz\
     -O /scratch/GenomePipe/combinedGVCFs_1.batch.vcf
