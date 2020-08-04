from collections import defaultdict, Counter
import gzip
from Levenshtein import distance
import networkx as nx
import numpy as np
import tabix

addedLength = 0

svDict = defaultdict(list)

def pdist(x, dist=distance, threshold=7):
    
    #D = np.zeros(len(x), len(x))
    G = nx.graph()
    for i in range(len(x)):
        for j in range(i+1, len(x)):
            d = dist(x[i], x[j])
            if d<threshold:
                G.add_edge(i,j,weight)
    return G
            
class SeqCluster:
    def __init__(seqs):
        G = pdist(seqs)
        
                       
def stats(v):
    import pdb
    lens = np.array([len(sv.alt) for sv in v])
    seqFreq = Counter([sv.alt for sv in v])
    seq, freq = seqFreq.most_common(1)[0]

    print(seq)
    if str(seq)=="<INS>": ##ideally, we could get something out of this too
        pdb.set_trace()
    s = sorted(list(seqFreq.values()), reverse=True)[:3]
    return len(v), len(set([sv.alt for sv in v])), s, lens.mean(), lens.std()

def stats2vcf(v):
    seqFreq = Counter([sv.alt for sv in v])
    seq, freq = seqFreq.most_common(1)[0]
    if str(seq)=="<INS>": ##ideally, we could get something out of this too
        return
    rep = None
    if freq == len(v):
        rep = [sv for sv in v if sv.alt == seq][0]
        rep.freq = freq
        rep.freqCat = 'shared'
    elif freq > len(v)/2:
        rep = [sv for sv in v if sv.alt == seq][0] ## pick a represent. with the most frequent seq (not efficient)
        rep.freq = freq
        rep.freqCat = 'common'
    elif freq > 1:
        rep = [sv for sv in v if sv.alt == seq][0] ## pick a represent. with the most frequent seq (not efficient)
        rep.freq = freq
        rep.freqCat = 'polymorph'
    elif freq == 1:
        rep = [sv for sv in v if sv.alt == seq][0] ## pick a represent. with the most frequent seq (not efficient)
        rep.freq = freq
        rep.freqCat = 'singleton'
    return rep

def processReps(reps):
    for rep in reps:
        vcf = [rep.chrom, rep.start, rep.id, rep.ref, rep.alt, rep.qual, rep.filter, rep.info, rep.format, rep.sampleInfo]
        #print('\t'.join(vcf))
        
        #return len(rep.alt) - len(rep.ref)

class SV:
    def __init__(self, line, sample):
        self.chrom, self.start, self.id, self.ref, self.alt, self.qual, self.filter, self.info, self.format, self.sampleInfo = line.split()
        self.sample = sample
        self.svType = self.id.split(':')[0].split("Manta")[-1]

samples = 0
maxlen = 0
## produce diploidSV.txt with
## find /research/genomicds1/Manta/1*/results/variants -name "diploidSV.vcf.gz" > diploidSV.txt
for dfile in open("diploidSV.txt"):
    diploidfile = dfile.strip()
    samples += 1
    for line in gzip.open(diploidfile):
        line=line.decode().strip()
        if line.startswith('#'): continue
        sv = SV(line, diploidfile.split('/')[1])
        svDict[(sv.chrom, sv.start, sv.svType)].append(sv)
        ## printout for clustering
        #if sv.svType == 'INS' and len(sv.alt)>50:
        #    print(f'>{sv.chrom}|{sv.start}|{sv.id}|{sv.info}\n{sv.alt}')
        maxlen = max(maxlen, len(sv.alt))
if False:
    ## frequent Insertions (present in at least half the samples)
    freqIns = {k:v for k,v in svDict.items() if k[2]=='INS' and len(v)>samples/2}
    reps0 = [stats2vcf(v) for k,v in freqIns.items()]
    reps = [rep for rep in reps0 if not rep is None]
    #reps = sorted(reps, key=lambda rep: (int(rep.chrom[3:]), int(rep.start)))
    for rep in reps:
        if len(rep.alt) - len(rep.ref) >= 100:
            print(f'>{rep.chrom}|{rep.start}|{rep.id}|{rep.freq}|{rep.freqCat}\n{rep.alt}')
