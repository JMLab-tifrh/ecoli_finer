import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pandas import DataFrame as df

# defining rRNA operon locations
res = 500
rrnA = int(0.5*(4033262 + 4033379)/res)
rrnB = int(0.5*(4164390 + 4164507)/res)
rrnC = int(0.5*(3939539 + 3939656)/res)
rrnD = int(0.5*(3427069 + 3426962)/res)
rrnE = int(0.5*(4205886 + 4205994)/res)
rrnG = int(0.5*(2729470 + 2729354)/res)
rrnH = int(0.5*(223485 + 223593)/res)
oriC = int(0.5*(3925744 + 3925975)/res)

rrn_operons = {'rrnA':rrnA, 'rrnB':rrnB, 'rrnC':rrnC,
               'rrnD':rrnD, 'rrnE':rrnE, 'rrnG':rrnG,
               'rrnH':rrnH, 'oriC':oriC}

rnaSeq = pd.read_csv('rnaByDNA_scholes.csv', index_col=0)

genomicPos = np.array(rnaSeq['pos'])
propensity = np.array(rnaSeq['raw_propensity1'])

binsize = 500 # bp
prop_cg = np.zeros(int(4.64e6/binsize))

# coarse graining RNA-seq signal
for i in range(int(4.64e6/binsize)):
    _prop = propensity[np.where(np.logical_and(i*binsize<genomicPos, genomicPos<=(i+1)*binsize)==1)[0]]
    if(_prop.shape[0] > 1):
        prop_cg[i] = _prop.mean()
    else:
        prop_cg[i] = 0
prop_cg /= prop_cg.max() # scaling between 0 to 1

cutoff = 0
dcut = 0.01
err = 5
mean = 43
N = prop_cg.shape[0]

nRegions = 0

while(not mean-err<=nRegions<=mean+err):
    nRegions = 0
    
    labels = np.ones(N)
    labels[list(rrn_operons.values())] = 0
    labels[np.where(prop_cg > cutoff)[0]] = 0
    
    for i in range(N):
        nxt = i+1
        if(nxt+1>=N):
            nxt = nxt-N

        if(labels[nxt]==0 and labels[i]==1):
            nRegions += 1
            
    cutoff = round(cutoff+dcut, 3)
    
print(f"Cutoff={cutoff:.3f}, regions={nRegions:.3f}")
labels = labels.astype('int')

labels_rnaSeq = np.copy(labels)
np.savetxt('labels_rnaSeq.txt', labels_rnaSeq, fmt='%d')
print("labels saved in labels_rnaSeq.txt")
