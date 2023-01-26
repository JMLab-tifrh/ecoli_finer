import numpy as np
import matplotlib.pyplot as plt
from pandas import DataFrame as df
from numpy.linalg import norm
from numpy.random import exponential as exp
from numpy.random import random, randint
import sys

# to write gro files
#====================================================================================
def write_gro(out, pos, resids=[], resnames=[], at_names=[], box=[], mesg="Homebrewed GRO writer"):
    R"""Usage:
write_gro(out, pos, resids, resname, at_names, box)
---------------------------------------------------------------------
out              :  output filename
pos              :  positions of particles. ndarray of shape (n, 3)
                where n is the number of particles.
resids           :  residue id. integer. 
resnames         :  residue names for each particle. list of strings/characters of length n. 
at_names         :  atom names for each particle. list of strings/characters of length n. 
box              :  box dimensions (default = (0, 0, 0))
"""
    usage = """Usage:
write_gro(out, pos, resids, resname, at_names, charges, box)
---------------------------------------------------------------------
out              :  output filename
pos              :  positions of particles. ndarray of shape (n, 3)
                    where n is the number of particles.
resids         :  particle types. list of strings/characters of length n. 
resnames         :  residue names for each particle. list of strings/characters of length n. 
at_names         :  atom names for each particle. list of strings/characters of length n. 
box              :  box dimensions (default = (0, 0, 0))
"""

    if(type(resids) != type([])):
        resids = list([resids])
    if(type(resnames) != type([])):
        resnames = list([resnames])
    if(type(at_names) != type([])):
        at_names = list([at_names])

    n = pos.shape[0]

    if(len(resids) == 0):
        resids = [i for i in range(1, n+1)]
    elif(len(resids) == 1):
        resids = resids*n
    elif(len(resids) < n):
        print(f"Mismatch in size of resids and shape[0] of pos. Got {len(resids)}, expected {n}.")
        print(usage)

    if(len(resnames) == 0):
        resnames = ['UNK']*n
    elif(len(resnames) == 1):
        resnames = resnames*n
    elif(len(resnames) < n):
        print(f"Mismatch in size of resnames and shape[0] of pos. Got {len(resnames)}, expected {n}.")
        print(usage)

    if(len(at_names) == 0):
        at_names = ['X']*n
    elif(len(at_names) == 1):
        at_names = at_names*n
    elif(len(at_names) < n):
        print(f"Mismatch in size of at_names and shape[0] of pos. Got {len(at_names)}, expected {n}.")
        print(usage)

    if(len(box) == 0):
        box = [0, 0, 0]
    elif(len(box) != 3):
        print(f"Mismatch in size of box vectors. Got {len(box)}, expected 3.")
        print(usage)

    # for i in range(n):
    #   print(pos[i, 0])

    # num,type,residue number, residue name, atom name, x, y, z
    with open(out, "w") as w:
        w.write(f"{mesg}\n{n}\n")
        for i in range(n):
            w.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" % (resids[i], resnames[i],
                    at_names[i], (i+1)%100000, pos[i, 0], pos[i, 1], pos[i, 2]))
        w.write("   %.5lf %.5lf %.5lf\n" % (box[0], box[1], box[2]))

    print(f"{out} written successfully.")
#====================================================================================

usage = """usage: python build.py <labels from RNA-seq> <G>[default G=1.0]"""

try:
    labels = np.loadtxt(sys.argv[1], dtype='int')
except:
    print(usage)
    exit(0)

try:
    G = float(sys.argv[2])
except:
    try:
        G = 1.0
    except:
        print("Failed to read value of G.")
        print(usage)
        exit(0)

# parameters
res = 500  # bp
resOld = 5e3  # bp
oriC = int(0.5*(3925744+3925975)/res)
rep_regions = np.array([oriC-int(0.5*4.64e6*(G-1)/res), oriC+int(0.5*4.64e6*(G-1)/res)])
N = int(4.64e6/res)
for i in range(len(rep_regions)):
    if(rep_regions[i]>=N):
        rep_regions[i] = rep_regions[i] - N
labels[rep_regions[0]-1:rep_regions[0]+2] = 0
labels[rep_regions[1]-1:rep_regions[1]+2] = 0

# calculating PAR regions
nRegions = 0
for i in range(labels.shape[0]-1):
    if(labels[i]==0 and labels[i+1]==1):
        nRegions += 1
print(f'Number of PAR domains={nRegions}')

sig = 0.464
oldSig = 1
meanLength = int(10e3/res)
totalBeads = int(4640e3/res)
box = np.array([0.84, 0.84, 3.05])/68.210e-3
resRatio = int(resOld/res)
sig *= 1.1

while True:
    # generating plectoneme lengths
    expLengths = 151*np.ones(1000)
    while(expLengths.max() > 150):
        expLengths = (exp(scale=meanLength-8, size=1000)+8).astype('int')

    # generating coordinates

    indices = []
    for ctr, idx in enumerate(labels):
        indices.append(ctr)
    indices = np.array(indices, dtype='int')

    idxPARs = indices[labels == 1]
    ctr = 0
    idx = 0
    totalLength = idxPARs.shape[0]
    while idx < totalLength:
        labels[idxPARs[idx]] = -1
        idx += (expLengths[ctr])
        ctr += 1

    for i in range(labels.shape[0]-1):
        if(labels[i] == 1 and (labels[i+1] == 0 or labels[i+1] == -1)):
            labels[i] = 2
    labels[-1] = 2

    for i in range(labels.shape[0]-1):
        if(labels[i] == 0 and labels[i+1] == 1):
            labels[i+1] = -1

    for i in range(labels.shape[0]-1):
        if(labels[i] == 0 and labels[i+1] == 2):
            labels[i+1] = 0
    for i in range(labels.shape[0]-1):
        if(labels[i] == -1 and labels[i+1] == 0):
            labels[i] = 0
    try:
        plecLengths = (indices[labels == 2] - indices[labels == -1])+1
        break
    except:
        None

branchLengths = ((2147/res) + (plecLengths - (3500/res))
                 * 0.066775).astype('int')

pos = np.zeros((totalBeads, 3))

nPFR = labels[np.logical_or(labels == 0, labels == -1)].shape[0]
theta = np.linspace(0, 2*np.pi*0.9975, nPFR)
radius = 0.5*nPFR*sig/np.pi
x = radius*np.cos(theta)
y = radius*np.sin(theta)

for idx, X, Y in zip(indices[np.logical_or(labels == 0, labels == -1)], x, y):
    pos[idx, 0] = X
    pos[idx, 1] = Y

for start, end in zip(indices[labels == -1], indices[labels == 2]):
    init = pos[start]
    vec = init/norm(init)
    ctr = 1
    for i in range(start+1, end+1):
        pos[i] = init + ctr*sig*vec
        ctr += 1

trunkLength = []
nBranches = []
for ctr, (start, end) in enumerate(zip(indices[labels == -1], indices[labels == 2])):
    if(plecLengths[ctr] < 12):
        trunkLength.append(plecLengths[ctr])
        nBranches.append(0)
        continue

    rand = random()*0.1 + 0.5
    branchBeads = int(plecLengths[ctr]*(1-rand))
    trunkLength.append(int(plecLengths[ctr]*rand))
    nBranch = np.floor(branchBeads/branchLengths[ctr]).astype('int')
    nBranches.append(nBranch)

    if(nBranch == 0):
        continue

    difference = int((trunkLength[ctr]-3)/nBranch)
    if(difference < branchLengths[ctr]):
        difference = branchLengths[ctr]+1

    trunkBeads = [i for i in range(start, end)]
    locations = []

    for i in range(nBranch):
        loc = start+2 + i*(difference+1)

        if(loc >= end):
            loc = end - branchLengths[ctr]

        #print(ctr, loc, loc+branchLengths[ctr], start, end, nBranch, difference, branchLengths[ctr])

        init = pos[loc]
        vec = np.array([0, 0, (-1)**randint(0, 2)])

        for j in range(1, branchLengths[ctr]+1):
            labels[loc + j] = 3
            pos[loc + j] = init + j*vec*sig
            trunkBeads.remove(loc+j)
        last = loc + branchLengths[ctr]

        init = pos[trunkBeads[0]]
        vec = init/norm(init)
        for ctr1, j in enumerate(trunkBeads[1:]):
            pos[j] = init + sig*(ctr1+1)*vec

# writing gro file
nsr = int(603e3/res)
rt = int(1206e3/res)
ter = int(2041e3/res)
lt = int(2877e3/res)
nsl = int(3758e3/res)
ori = int(4640e3/res)

restypes = []
resnames = []
for ctr, l in enumerate(labels):
    if(l == 0):
        resnames.append('F')
    elif(l == -1):
        resnames.append('F')
    elif(l == 1):
        resnames.append('T')
    elif(l == 2):
        resnames.append('T')
    elif(l == 3):
        resnames.append('B')

    if(ctr+1 <= nsr):
        restypes.append(1)
    elif(ctr+1 > nsr and ctr+1 <= rt):
        restypes.append(2)
    elif(ctr+1 > rt and ctr+1 <= ter):
        restypes.append(3)
    elif(ctr+1 > ter and ctr+1 <= lt):
        restypes.append(4)
    elif(ctr+1 > lt and ctr+1 <= nsl):
        restypes.append(5)
    elif(ctr+1 > nsl and ctr+1 <= ori):
        restypes.append(6)
    else:
        restypes.append(0)

write_gro('pol.gro', pos, resids=restypes, resnames='POL1',
          at_names=resnames, box=box, mesg=f"Chromsome at 500 bp with G={G:.2f}")

# generating bonds
bonds = []

# backbone
backbone = indices[np.logical_or(labels == 0, labels == -1)]+1
for i in range(backbone.shape[0]-1):
    bonds.append([backbone[i], backbone[i+1]])
bonds.append([backbone[0], backbone[-1]])

# trunks
trunkBeads = indices[labels == 1]
branchBeads = indices[labels == 3]
for start, end in zip(indices[labels == -1], indices[labels == 2]):
    beads = [i for i in range(start, end+1)]
    for b in branchBeads:
        if(b in beads):
            beads.remove(b)
    for i in range(len(beads)-1):
        bonds.append([beads[i]+1, beads[i+1]+1])

for i in range(labels.shape[0]-1):
    if((labels[i] == 1 or labels[i] == 3) and labels[i+1] == 3):
        bonds.append([i+1, i+2])

bonds = np.array(bonds, dtype='int')
np.savetxt('bonds.txt', bonds, fmt='%d')
np.savetxt('labels.txt', labels, fmt='%d')

