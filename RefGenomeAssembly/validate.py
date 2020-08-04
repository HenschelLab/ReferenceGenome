# by maintaining an offset, we can check the corresponding positions between hg19 and hg19uae

from Bio import SeqIO

alleleFreqDir = '/research/btc_bioinformatic/results/AlleleFrequencies'
applicationDir = '/research/btc_bioinformatic/operations/'
hg19='%s/data_masdar/ucsc.hg19.fasta' % applicationDir

for refChrom in SeqIO.parse(hg19, 'fasta'):
    if refChrom.description != 'chr1': continue
    break
for uaeChrom in SeqIO.parse(f'{alleleFreqDir}/hg19uae.fasta', 'fasta'):
    if uaeChrom.description != 'chr1': continue
    break



chromFile = f'{alleleFreqDir}/mafUAE_{uaeChrom.description}.txt'
counter = 0
offset = 0
for line in open(chromFile):
    pos1, ref, frq, alt = line.split()
    pos = int(pos1)-1 # python counting
    posoff = pos + offset
    print (pos, refChrom.seq[pos : pos + len(ref)], posoff, uaeChrom.seq[posoff: posoff + len(alt)])
    offset += (len(alt) - len(ref))
    
    counter += 1
    if counter == 10: break
