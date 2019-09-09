# Checking coverage on modified loci only, apparently that's what Fakhro et al (Qatar paper) did
# by maintaining an offset, we can check the corresponding positions between hg19 and hg19uae
# run this to check the first 10 positions
# !head /research/btc_bioinformatic/results/AlleleFrequencies/mafUAE_chr1.txt
# and then run the script
from Bio import SeqIO
from subprocess import check_output

alleleFreqDir = '/research/btc_bioinformatic/results/AlleleFrequencies'
applicationDir = '/research/btc_bioinformatic/operations/'
hg19 = '%s/data_masdar/ucsc.hg19.fasta' % applicationDir
refbam = '/research/btc_bioinformatic/operations/scratch/13164/Alignment/13164_001_FR_hg19.paired_aligned_Sorted_Merged_dedup.bam'
ourbam = '/research/btc_bioinformatic/operations/UaeRef/13164/Alignment/13164_001_FR_hg19uae.paired_aligned_Sorted_Merged_dedup.bam'
samtools = '/apps/samtools/samtools-1.2-gcc-4.9.2/bin/samtools'

selChrom = 'chr20'
## Just look at chr1 for the time being:
for refChrom in SeqIO.parse(hg19, 'fasta'):
    if refChrom.description != selChrom: continue
    break
for uaeChrom in SeqIO.parse(f'{alleleFreqDir}/hg19uae.fasta', 'fasta'):
    if uaeChrom.description != selChrom: continue
    break
#samtools view 13164_001_FR_hg19uae.paired_aligned_Sorted_Merged_dedup.bam chr1:723798-723800|wc -l 


chromFile = f'{alleleFreqDir}/mafUAE_{uaeChrom.description}.txt'
counter = 0
offset = 0
depthRef = []
depthUAE = []
for line in open(chromFile):
    pos1, ref, frq, alt = line.split()
    pos = int(pos1)-1 # python counting
    posoff = pos + offset    
    if float(frq)>.97 and len(ref) != 0 and len(alt) != 0: 
        samviewRef = check_output(f'{samtools} view {refbam} {selChrom}:{pos+1}-{pos+1+len(ref)}', shell=True)
        samviewUAE = check_output(f'{samtools} view {ourbam} {selChrom}:{posoff+1}-{posoff+1+len(alt)}', shell=True)
        depthReflocus = len(samviewRef.decode().split('\n'))
        depthUAElocus = len(samviewRef.decode().split('\n'))
        print (int(samviewRef), refChrom.seq[pos : pos + len(ref)], int(samviewUAE), uaeChrom.seq[posoff: posoff + len(alt)], frq)
        depthRef.append(depthReflocus)
        depthUAE.append(depthUAElocus)
        for ind, seq in zip(indent, seqs):
            print(' '*ind + seq) 
        reads =samviewRef.decode().rstrip().split('\n')
        indent = [int(read.split()[3]) - 80512 for read in reads]
        for ind, seq in zip(indent, seqs):
            print(' '*ind + seq) 
        counter += 1
        if counter == 1: break
        
    offset += (len(alt) - len(ref)) ## always update the offset!!
    
