#BSUB -n 24
#BSUB -q general
#BSUB -W 3:00

#BSUB -o /research/genomicds1/Manta/JobOutputs/manta.%I.%J.out
#BSUB -e /research/genomicds1/Manta/JobOutputs/manta.%I.%J.err
#BSUB -J manta[1-120]

echo "Job: ${LSB_JOBINDEX}"
conda activate manta
cd /research/genomicds1/Manta/${LSB_JOBINDEX};
./runWorkflow.py
