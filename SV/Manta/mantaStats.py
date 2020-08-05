import numpy as np
from collections import defaultdict
l = defaultdict(list)
for line in open("mantaStats.txt"):
    if line.startswith("/research"): continue
    nr, svtype = line.strip().split()[:2]
    l[svtype].append(int(nr))
for k,v in l.items():
    a = np.array(v)
    print(k, len(a), a.mean(), a.std(), a.sum())

total = [sum([l[key][i] for key in l.keys()]) for i in range(118)]
print(np.array(total).mean())

    """
In [19]: run mantaStats.py
             avg               std                sum
MantaBND 118 1638.677966101695 296.50680620295515 193364
MantaDEL 118 4824.872881355932 490.5836293642221 569335
MantaDUP 118 477.9406779661017 94.70241458565837 56397
MantaINS 118 3095.7796610169494 653.5084341342741 365302
"""

#avg SV:
#Out[24]: 10037.271186440677
    
