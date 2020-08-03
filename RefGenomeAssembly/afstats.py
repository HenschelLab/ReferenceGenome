"""
### Allele frequency statistics
Reads out frq file (produced from VCF using vcftools or so) 
(simply contains chr, pos, ref:freq0, alt1:freq1, alt2:freq2...)
Also includes info for multiallelic loci

The script identifies the alt_i with max freq_i
final result: freqCounter - a histogram counter of 5% sized bins,
accounting how frequent the (max) alt. freqeuncy for a locus is
"""

import numpy as np
from collections import Counter

wdir = '/research/btc_bioinformatic/results'
afdir = wdir + '/AlleleFrequencies'
frqFile = 'ALL_PASS_only_AllelFreq.frq'
prevChrom = None
w = None
freqs = []
freqCounter = Counter()

for line in open(f'{wdir}/{frqFile}'):    
    fields = line.rstrip().split('\t')
    chrom, pos, allelic, count = fields[:4]
    if chrom == 'CHROM': continue
    ref = fields[4].split(':')[0]
    alts0 = [f.split(':') for f in fields[5:]]
    maxalt = max([(float(frq), alt) for (alt, frq) in alts0])
    #freqs.append(maxalt[0])
    category = int(maxalt[0]*20)
    freqCounter[category] += 1
#hist = np.histogram(list(zip(*freqs))[0], bins = np.arange(0, 1.1,0.05))
#print(hist)
#    if maxalt[0] > 0.5 and int(count) > 100:
#        pass




