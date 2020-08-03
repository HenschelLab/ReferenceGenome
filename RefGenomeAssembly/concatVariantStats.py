## This script simply appends individual calls/chromosome stats from samples, against ref
## Result is written to variantReduction.csv
## individual tables were created using the shell script gatherNrVariantsInVCFs.sh

import pandas as pd
import glob

refFiles = sorted(glob.glob('variantStats_*_scratch.csv'))
uaeFiles = sorted(glob.glob('variantStats_*_UaeRef.csv'))
#refFiles = sorted(glob.glob('variantStats_*_scratch.csv'))

dfs = [pd.read_csv(refFile, sep='\s+', index_col=1, header=None) for refFile in uaeFiles+refFiles]
mega=pd.concat(dfs, axis=1, sort=True)

refn = [rf.split('_')[1]+'R' for rf in refFiles]
uaen = [rf.split('_')[1]+'U' for rf in uaeFiles]
diffn = [rf.split('_')[1]+'D' for rf in refFiles]

mega.columns = uaen + refn

chroms = [f'chr{id}' for id in list(range(1,23)) + ['X','Y']]
df = mega[sorted(uaen + refn)].loc[chroms]
df.to_csv('variantReduction.csv')
