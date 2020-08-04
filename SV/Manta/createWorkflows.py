import subprocess
import os

REFGENOME="/research/btc_bioinformatic/operations/data_masdar/ucsc.hg19.fasta"
mantaDir = "/research/genomicds1/Manta"
bsub = """#BSUB -n 24
#BSUB -q general
#BSUB -W 3:00

#BSUB -o /research/genomicds1/Manta/JobOutputs/manta.%%I.%%J.out
#BSUB -e /research/genomicds1/Manta/JobOutputs/manta.%%I.%%J.err
#BSUB -J manta[%s]                        
echo "Job: ${LSB_JOBINDEX}"
conda activate manta
cd /research/genomicds1/Manta/${LSB_JOBINDEX};
./runWorkflow.py
"""

samples = []
for line in open("bamfiles.txt"):
    bamfile = line.strip()
    sample = int(bamfile.split('/')[-1].split('_')[0])

    samples.append(sample)
    sampleDir = '%s/%s' %(mantaDir, sample)
    if not os.path.exists(sampleDir):
        os.mkdir(sampleDir)
    mantaCmd = "configManta.py --bam %s --referenceFasta %s  --runDir %s" %(bamfile, REFGENOME, sampleDir)
    
    print (mantaCmd)
    subprocess.call(mantaCmd, shell=True)
    

print(bsub % ",".join(map(str,sorted(samples))))
