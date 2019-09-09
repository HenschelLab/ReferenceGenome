'''
Generating the actual reference genome, deriving it from the reference hg19

Flipping those alleles where the max. alternative call is >50% in the vcf file from 
joint genotype calling.

Frequencies are taken from freq files, that are split by chromosome (see file splitChromsEtc.py)
(as the algorithm operates chromosome-wise)

Uses Biopython for fasta seq handling

Paralellized per chromosome (chr-id comes as command line, 1..24 -> 1..22,X,Y

TODO: create a liftover file while creating the new reference genome 
'''

import pdb
from Bio import SeqIO
#from Bio.Alphabet import generic_dna
import sys
alleleFreqDir = '/research/btc_bioinformatic/results/AlleleFrequencies'
applicationDir = '/research/btc_bioinformatic/operations/'
hg19='%s/data_masdar/ucsc.hg19.fasta' % applicationDir
hg19Y = 'hg19_chrY.fasta' # for testing purposes

uaeref = ''
lastPosRef = 0
verbose = True
job = int(sys.argv[-1]) - 1 
chromJob = [f'chr{id}' for id in list(range(1,23)) + ['X','Y']][job]
for chrom in SeqIO.parse(hg19, 'fasta'): 
    if chrom.description != chromJob: continue #== 'chrM': continue ## no mitoch.
    print(f'Extracted {chrom.description}: Length {len(chrom.seq)}')
    chromFile = f'{alleleFreqDir}/mafUAE_{chrom.description}.txt'
    for line in open(chromFile):
    #ref, pos1, alt in indels:
        pos1, ref, frq, alt = line.split()
        pos = int(pos1)-1 # python counting
        ## making sure that the REF entry from our VCF actually is the same as in the hg19 fasta
        if verbose and ref.upper() != str(chrom.seq[pos:pos+len(ref)]).upper():
            pdb.set_trace()

        ## ourref: Bio.Seq dtype, nice: takes description from chrom
        uaeref += chrom[lastPosRef: pos].seq + alt 
        if verbose:
            before = chrom[pos-10 : pos]
            after = chrom[pos+len(ref) : pos+len(ref)+10]
            print (f'{before} > {ref} < {after}')
            print (f'{ourref[-10:]} > {alt}')
            print (line.rstrip())
            #pdb.set_trace()

        lastPosRef = pos + len(ref)
    uaeref += chrom[lastPosRef:]
    ## writing/appending chromosome sequence to disc    
    with open(f'{alleleFreqDir}/hg19uae_{chrom.description}.fasta', 'w') as w:
        SeqIO.write(uaeref, w, "fasta")
        #w.write(f'>{chrom.description}\n{ourref}')

    