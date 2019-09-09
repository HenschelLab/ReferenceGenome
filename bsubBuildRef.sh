#BSUB -n 12
#BSUB -q general
#BSUB -W 48:00

#BSUB -o JobOutputs/flipIndels.%I.%J.out
#BSUB -e JobOutputs/flipIndels.%I.%J.err
#BSUB -J flipIndels[1-24]

#echo "Job: ${LSB_JOBINDEX}"
source activate bio3
python buildRefParallel.py ${LSB_JOBINDEX}
