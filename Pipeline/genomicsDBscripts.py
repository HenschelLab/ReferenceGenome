import os, time

bsub = """#BSUB -n 24
#BSUB -q general                           
#BSUB -W 48:00      
#BSUB -J gDB%s
#BSUB -o GenomicsDB_Jobout/gDB%s.out
#BSUB -e GenomicsDB_Jobout/gDB%s.err

module load gatk/4.0.6.0

time gatk --java-options "-Xmx128G" GenomicsDBImport \
 --variant /scratch/GenomePipe/Batch1_GVCFs/120067_HaplotypeCaller.gvcf.gz\
 --variant /scratch/GenomePipe/Batch1_GVCFs/120165_HaplotypeCaller.gvcf.gz\
 --variant /scratch/GenomePipe/Batch1_GVCFs/13024_HaplotypeCaller.gvcf.gz\
 --variant /scratch/GenomePipe/Batch1_GVCFs/120060_HaplotypeCaller.gvcf.gz\
 --variant /scratch/GenomePipe/Batch1_GVCFs/12608_HaplotypeCaller.gvcf.gz\
 --variant /scratch/GenomePipe/Batch1_GVCFs/120129_HaplotypeCaller.gvcf.gz\
 --variant /scratch/GenomePipe/Batch1_GVCFs/120287_HaplotypeCaller.gvcf.gz\
 --variant /scratch/GenomePipe/Batch1_GVCFs/120136_HaplotypeCaller.gvcf.gz\
 --variant /scratch/GenomePipe/Batch1_GVCFs/120069_HaplotypeCaller.gvcf.gz\
 --variant /scratch/GenomePipe/Batch1_GVCFs/120195_HaplotypeCaller.gvcf.gz\
 --variant /scratch/GenomePipe/Batch1_GVCFs/12964_HaplotypeCaller.gvcf.gz\
 --variant /scratch/GenomePipe/Batch1_GVCFs/120197_HaplotypeCaller.gvcf.gz\
 --variant /scratch/GenomePipe/Batch1_GVCFs/120131_HaplotypeCaller.gvcf.gz\
 --variant /scratch/GenomePipe/Batch1_GVCFs/120164_HaplotypeCaller.gvcf.gz\
 --variant /scratch/GenomePipe/Batch1_GVCFs/10885_HaplotypeCaller.gvcf.gz\
 --variant /scratch/GenomePipe/Batch1_GVCFs/13116_HaplotypeCaller.gvcf.gz\
 --variant /scratch/GenomePipe/Batch1_GVCFs/120133_HaplotypeCaller.gvcf.gz\
 --variant /scratch/GenomePipe/Batch1_GVCFs/120190_HaplotypeCaller.gvcf.gz\
 --variant /scratch/GenomePipe/Batch1_GVCFs/120075_HaplotypeCaller.gvcf.gz\
 --variant /scratch/GenomePipe/Batch1_GVCFs/13049_HaplotypeCaller.gvcf.gz\
 --variant /scratch/GenomePipe/Batch1_GVCFs/120213_HaplotypeCaller.gvcf.gz\
 --variant /scratch/GenomePipe/Batch1_GVCFs/12960_HaplotypeCaller.gvcf.gz\
 --variant /scratch/GenomePipe/Batch1_GVCFs/12553_HaplotypeCaller.gvcf.gz\
 --variant /scratch/GenomePipe/Batch1_GVCFs/13119_HaplotypeCaller.gvcf.gz\
 --variant /scratch/GenomePipe/Batch1_GVCFs/12657_HaplotypeCaller.gvcf.gz\
 --genomicsdb-workspace-path /scratch/GenomePipe/GenomicsDB/%s \
 --intervals %s
"""
for chrm in list(range(1,21)) + ['X', 'Y']:
    chrom = 'chr%s' %chrm
    scriptfile = '/scratch/GenomePipe/GenomicsDB/bsubA_%s.sh' % chrom
    with open(scriptfile, 'w') as w:
        w.write(bsub % (chrm, chrom, chrom, chrom, chrom))

    time.sleep(1)
    cmd = 'bsub < %s' % scriptfile
    print(cmd)
    os.system(cmd)
    
