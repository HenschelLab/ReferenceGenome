"""
GenomicsDBImport bsub script preparation and execution

Generates bsub files for region/chromosome-wise import to genomicsDB
Splitting by region brings the job size below any given time limit (48h).

Note, some chromosome numberings start far into the chromosome. this will yield bogus jobs...
Also, would be nice to skip centromeric regions.

The follow-up script (actually contains a bsub completion dependency) is:
genotypeGVCFsScripts.py

"""


import os, time
from subprocess import call
import glob

bsub = """#BSUB -n 4
#BSUB -q general                           
#BSUB -W 24:00      
#BSUB -J gDB%s
#BSUB -o /research/btc_bioinformatic/results/GenomicsDB_JoboutFamilies/gDB%s.out
#BSUB -e /research/btc_bioinformatic/results/GenomicsDB_JoboutFamilies/gDB%s.err

module load gatk/4.0.6.0

time gatk --java-options "-Xmx12G" GenomicsDBImport %s\
 --genomicsdb-workspace-path /%s/%s \
 --intervals %s
"""
#553636 - 553957
genomicsDBdir = '/research/btc_bioinformatic/results/GenomicsDB'

## collecting all gvcfs from various batches (adjusting to various folder naming conventions)

exomebatch = glob.glob('/research/btc_bioinformatic/results/scratch/*/GATK/*.gvcf.gz')
variantOptions = ' '.join(['--variant %s'%gvcf for gvcf in exomebatch])

chromLen = {'20': 63025520, '21': 48129895, '22': 51304566, '1': 249250621, '3': 198022430, '2': 243199373, '5': 180915260, '4': 191154276, '7': 159138663, '6': 171115067, '9': 141213431, '8': 146364022, 'Y': 59373566, 'X': 155270560, '11': 135006516, '10': 135534747, '13': 115169878, '12': 133851895, '15': 102531392, '14': 107349540, '17': 81195210, '16': 90354753, '19': 59128983, '18': 78077248}

regionSize = 10000000
for chrm in list(range(1,23)) + ['X', 'Y']:
    chromL = chromLen[str(chrm)]
    for regidx, start in enumerate(range(4, chromL)[::regionSize]):
        interval = 'chr%s:%s-%s' % (chrm, start, min(start+regionSize, chromL))
                              
        chrom = 'chr%s' %chrm
        region = '%s_%03d' %(chrom, regidx)
#        if os.path.exists('/scratch/GenomePipe/GenotypeGVCFs/%s.vcf.idx' % region):
#            print ("Job %s previously succeeded" % region)
#            continue
#        print ("Job %s previously failed" % region)
#        call(["rm", "-r", '%s/%s/'%(genomicsDBdir, region)])
        
        scriptfile = '%s/Scripts/bsubB_%s.sh' % (genomicsDBdir, region)
        with open(scriptfile, 'w') as w:
            w.write(bsub % (region, region, region, variantOptions, genomicsDBdir, region, interval))

        time.sleep(0.3)
        cmd = 'bsub < %s' % scriptfile
        print(cmd)
        os.system(cmd)

    

