'''
Generating the actual reference genome, deriving it from the reference hg19

Flipping those alleles where the max. alternative call is >50% in the vcf file from 
joint genotype calling.

Frequencies are taken from freq files, that are split by chromosome
(as the algorithm operates chromosome-wise)
'''

import pdb
from Bio import SeqIO
#from Bio.Alphabet import generic_dna

alleleFreqDir = '/research/btc_bioinformatic/results/AlleleFrequencies'
applicationDir = '/research/btc_bioinformatic/operations/'
hg19='%s/data_masdar/ucsc.hg19.fasta' % applicationDir
hg19Y = 'hg19_chrY.fasta' # for testing purposes

ourref = ''
lastPosRef = 0
verbose = True
for chrom in SeqIO.parse(hg19, 'fasta'): ## CHANGE!!! 
    if chrom.description == 'chrM': continue ## no mitoch.
    chromFile = f'{alleleFreqDir}/mafUAE_{chrom.description}.txt'
    for line in open(chromFile):
    #ref, pos1, alt in indels:
        pos1, ref, frq, alt = line.split()
        pos = int(pos1)-1 # python counting

        if verbose and ref.upper() != str(chrom.seq[pos:pos+len(ref)]).upper():
            pdb.set_trace()

        ourref += chrom[lastPosRef: pos] + alt
        #if verbose:
        #    before = chrom[pos-10 : pos]
        #    after = chrom[pos+len(ref) : pos+len(ref)+10]
        #    print (f'{before} > {ref} < {after}')
        #    print (f'{ourref[-10:]} > {alt}')
        #    print (line.rstrip())
        #    pdb.set_trace()

        lastPosRef = pos + len(ref)
    ourref += chrom[lastPosRef:]
    ## writing/appending chromosome sequence to disc    
    with open(f'{alleleFreqDir}/hg19arab.fasta', 'w+') as w:
        w.write(f'>{chrom.description}\n{ourref}')
    print(f'Appended {chrom.description}')
    