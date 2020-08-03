from Bio import SeqIO
from Bio.Seq import MutableSeq

alleleFreqDir = '/research/btc_bioinformatic/results/AlleleFrequencies'
applicationDir = '/research/btc_bioinformatic/operations/'
hg19='%s/data_masdar/ucsc.hg19.fasta' % applicationDir
warning = 0
mods = 0
for chrom in SeqIO.parse(hg19,'fasta'):
    if chrom.description != 'chrY': continue ## CHANGE!!!                                                              
    chrom1 = chrom.seq
    majorAllelesChrom = 'majorAllelesUAE_chr10.txt'
    if True:
        _,pos, altAllele = row.split()[:3]
        pos = int(pos)
        altAllele, frq = altAllele.split(':')
        if altAllele in list('ACGT'):
            if chrom1[pos] == altAllele:
                warning += 1
            else:
                chrom1[pos] = altAllele
                mods += 1

    with open(f'{alleleFreqDir}/hg19arab.fasta', 'w+') as w:
        w.write(f'>{chrDesc}\n{chrom1.seq}')