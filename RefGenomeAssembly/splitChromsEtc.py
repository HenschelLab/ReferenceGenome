#import re
#singleDigitChr = re.compile('chr(\d)$')
"""
This script takes a frequency file as produced by vcftools.
Assumes that entries are sorted by (chromosome, position)

creates a new tab-sep. file per chromosome, containing pos, ref, max-alt, max-altfreq, if max-alt > 50%
"""
wdir = '/research/btc_bioinformatic/results'
afdir = wdir + '/AlleleFrequencies'
frqFile = 'ALL_PASS_only_AllelFreq.frq'
prevChrom = None
w = None
for line in open(f'{wdir}/{frqFile}'):    
    fields = line.rstrip().split('\t')
    chrom, pos, allelic, count = fields[:4]
    if chrom == 'CHROM': continue
    ref = fields[4].split(':')[0]
    alts0 = [f.split(':') for f in fields[5:]]
    maxalt = max([(float(frq), alt) for (alt, frq) in alts0])
    if chrom != prevChrom: ## encountering a new chromosome
        if not w is None:
            w.close()
        #match = singleDigitChr.match(chrom)
        #chromf = f'chr0{match.group(1)}' if match else chrom
        chromfile = f'{afdir}/mafUAE_{chrom}.txt'
        print(f'Opening {chromfile}')
        w = open(chromfile, 'w')
    if maxalt[0] > 0.5 and int(count) > 100:
        w.write(f'{pos}\t{ref}\t{maxalt[0]}\t{maxalt[1]}\n')

    prevChrom = chrom
w.close()





