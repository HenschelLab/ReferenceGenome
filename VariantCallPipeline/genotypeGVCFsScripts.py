import os, time

bsub = """#BSUB -n 6
#BSUB -q general                           
#BSUB -W 48:00 
#BSUB -J gtyp_%s
#BSUB -o GenomicsDB_Jobout/genotypeGVCFs_%s.%%J.out
#BSUB -e GenomicsDB_Jobout/genotypeGVCFs_%s.%%J.err

module load gatk/4.0.6.0

time gatk --java-options "-Xmx12G" GenotypeGVCFs \
  -R /research/btc_bioinformatic/operations/data_masdar/ucsc.hg19.fasta\
  -V gendb:///research/btc_bioinformatic/results/GenomicsDB/%s \
  -G StandardAnnotation -new-qual \
  -O /research/btc_bioinformatic/results/GenotypeGVCFs/%s.vcf 
"""


chromLen = {'20': 63025520, '21': 48129895, '22': 51304566, '1': 249250621, '3': 198022430, '2': 243199373, '5': 180915260, '4': 191154276, '7': 159138663, '6': 171115067, '9': 141213431, '8': 146364022, 'Y': 59373566, 'X': 155270560, '11': 135006516, '10': 135534747, '13': 115169878, '12': 133851895, '15': 102531392, '14': 107349540, '17': 81195210, '16': 90354753, '19': 59128983, '18': 78077248}

regionSize = 10000000

genotypeGVCFdir = '/research/btc_bioinformatic/results/GenotypeGVCFs/'


for chrm in list(range(1,23)) + ['X', 'Y']: 
    chromL = chromLen[str(chrm)]
    for regidx, start in enumerate(range(chromL)[::regionSize]):
        #if not(chrm==1 and regidx==26): continue ## CHANGE!!!
        interval = 'chr%s:%s-%s' % (chrm, start, min(start+regionSize, chromL))
        
        chrom = 'chr%s' %chrm
        region = '%s_%03d' % (chrom, regidx)
        if os.path.exists('/scratch/GenomePipe/GenotypeGVCFs/%s.vcf.idx' % region):
            

            scriptfile = f'{genotypeGVCFdir}/Scripts/bsub_{region}.sh' #% region
            with open(scriptfile, 'w') as w:
                w.write(bsub % (region, region, region, region, region))

            time.sleep(0.3)
            cmd = 'bsub < %s' % scriptfile
            print(cmd)
            os.system(cmd)

