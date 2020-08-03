import pdb
#            ....5....0....5....0....5....0
hg19 =      "AAAAABBBBBCCCCCDDDDDEEEEEFFFFF"
newgenome = "AAAacctABBBCCCCDDDDDEEEEEffF"
#           'AAAacctABBCCCCCDDDDDEEEEEff'
#            ....5....0....5....0....5....0
segments = []
# 1 - 1
# 2 - 2
# 3 - 3
# 4 - 4,5,6,7
# 5 - 8 A
# 6 - 9 B
# 7 - 10 B
# 8 - 11 B
# 9 - del
# 10 - del
# 11 - del
# 12 - 12  ... 
indels = [('A', 4, 'acct'), ('BBC', 9, ''), ('FFFF', 26, 'ff')]
our = ''
lastPosRef = 0
for ref, pos1, alt in indels:
    pos = pos1-1 # python counting
    if ref != hg19[pos:pos+len(ref)]:
        pdb.set_trace()
    else:
        print(ref, hg19[pos:pos+len(ref)])
    our += hg19[lastPosRef: pos] + alt
    lastPosRef = pos + len(ref)
our += hg19[lastPosRef:]

print (newgenome)
print ("".join(our))

#    segments.append((lastPos, pos, 'ref')) # (0,4) AAA
#    segments.append(())
    