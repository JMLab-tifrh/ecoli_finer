import numpy as np
from MDAnalysis import Universe
from groio import write_gro
import sys
from glob import glob

#==========================================================================
# function to write gro files
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
#==========================================================================


usage = '''python genRepFork.py <unreplicated chromosome coordinates> <w>
Example:
python genRepFork pol.gro 1.0
'''

if(len(sys.argv) < 2):
    print(usage)
    exit(0)

try:
    G = float(open(sys.argv[1], 'r').readline().split()[-1][2:])
    pos = Universe(sys.argv[1]).trajectory[0].positions*0.1
    top = Universe(sys.argv[1])
except:
    print("Could not read "+sys.argv[1]+".")
    print(usage)
    exit(0)

try:
    labels0 = np.loadtxt('labels.txt', dtype='int')
except:
    print('labels.txt not found. Exiting.')
    print(usage)
    exit(0)

try:
    bonds0 = np.loadtxt('bonds.txt', dtype='int')
except:
    print('bonds.txt not found. Exiting.')
    print(usage)
    exit(0)

print('Replicating chromosome......')
# generating rep fork coordinates
names0 = top.select_atoms('all').names.tolist()
resids0 = top.select_atoms('all').resids.tolist()
index0 = top.select_atoms('all').resids


res = 500 # bp
oriC = int(0.5*(3925744+3925975)/res)
rep_regions = np.array([oriC-int(0.5*4.64e6*(G-1)/res), oriC+int(0.5*4.64e6*(G-1)/res)])
N = int(4.64e6/res)

_pos = np.copy(pos)

rep_regions[rep_regions>=N] -= N 
    
names = np.copy(names0).tolist()
resids = np.copy(resids0).tolist()
bonds = np.copy(bonds0)
labels = np.copy(labels0)    
    
indices = np.array([i for i in range(_pos.shape[0])])

if(G > 1):
    _temp = 0  
    
    displacement = np.array([0, 0, 10])
    if(oriC+int(0.5*4.64e6*(G-1)/res) > N):
        new_pos = np.copy(pos[rep_regions[0]-1:]+displacement)#@Rx
        replicated = np.copy(index0[rep_regions[0]-1:])
        replicated = np.append(replicated, index0[:rep_regions[1]])
        new_pos = np.append(new_pos, pos[:rep_regions[1]]+displacement, axis=0)
        _pos = np.append(_pos, new_pos, axis=0)
        
        labels = np.append(labels, labels[rep_regions[0]-1:])
        _temp= np.where(labels==-1)[0][-1]
        labels = np.append(labels, labels[:rep_regions[1]])
        
        names.extend(names0[rep_regions[0]-1:])
        names.extend(names0[:rep_regions[1]])
        
        resids.extend(resids0[rep_regions[0]-1:])
        resids.extend(resids0[:rep_regions[1]])
        
    else:
        new_pos = np.copy(pos[rep_regions[0]-1:rep_regions[1]])
        _pos = np.append(_pos, new_pos+displacement, axis=0)
        
        replicated = np.copy(index0[rep_regions[0]-1:rep_regions[1]])
        
        names.extend(names0[rep_regions[0]-1:rep_regions[1]])
        
        resids.extend(resids0[rep_regions[0]-1:rep_regions[1]])
        
        labels = np.append(labels, labels[rep_regions[0]-1:rep_regions[1]])
        
    write_gro('pol_fork.gro', _pos, resids=resids, resnames='POL1',
              at_names=names, box=top.trajectory[0]._unitcell[:3]*0.1,
              mesg=f"Chromsome at 500 bp with G={G:.2f}")
    
    print('Joining chromosome beads by bonds......')
    # generating rep fork backbone bonds
    bonded = [f'{min(b)}-{max(b)}' for b in bonds]
    indices = np.array([i+1 for i in range(pos.shape[0])])
    
    if(oriC+int(0.5*4.64e6*(G-1)/res) > N):
        indices = np.append(indices, indices[rep_regions[0]-1:])
        indices = np.append(indices, indices[:rep_regions[1]])
        for i in range(rep_regions[0]-1, N):
            for j in range(1, 100):
                if(f'{i}-{i+j}' in bonded):
                    bonds = np.append(bonds,
                                      np.array([[i,i+j]])-rep_regions[0]+N+1,
                                      axis=0)
                    
        shift = N + (N-rep_regions[0])
        bonds = np.append(bonds, np.array([[_temp-1, shift+1]]), axis=0)
        for i in range(rep_regions[1]):
            for j in range(1,  100):
                if(f'{i}-{i+j}' in bonded):
                    bonds = np.append(bonds,
                                      np.array([[i,i+j]])+shift+1,
                                      axis=0)
    else:
        indices = np.append(indices, indices[rep_regions[0]-1:rep_regions[1]])
        
        for i in range(rep_regions[0]-1, rep_regions[1]):
            for j in range(1, 100):
                if(f'{i}-{i+j}' in bonded):
                    bonds = np.append(bonds,
                                      np.array([[i,i+j]])-rep_regions[0]+N+1,
                                      axis=0)
    
    
    
    bonds = np.append(bonds, np.array([[rep_regions[0], N+1]]), axis=0)
    bonds = np.append(bonds, np.array([[rep_regions[1], _pos.shape[0]]]), axis=0)
    
    bonded = [f'{min(b)}-{max(b)}' for b in bonds]
    for i, j in zip(np.where(np.array(names)=='F')[0], np.where(np.array(names)=='F')[0][1:]):
        if(not f'{i+1}-{j+1}' in bonded and i+1 > N and j+1 > N):
            bonds = np.append(bonds, np.array([[i+1, j+1]]), axis=0)
            
    np.savetxt('bonds_fork.txt',  bonds, fmt='%d')
    np.savetxt('labels_fork.txt', labels.T, fmt="%d")
    np.savetxt('indices_fork.txt', indices.T, fmt="%d")
    
print('Introducing Hi-C bonds......')
# generating mother chromosome HiC bonds
mats = glob('expt_pij*')
print('Autodetecting experimental matrices present in the current directory.')
if(len(mats)==0):
    print('No matrix has been autodetected. Kindly enter path to a matrix. Press ENTER to quit.')
    path = input('path: ')
    if(path==""):
        exit(0)
    try:
        mats = np.loadtxt(path)
    except:
        print('Could not load matrix from provided path. Exiting.')
        exit(0)
elif(len(mats)>1):
    print('Multiple experimental Hi-C matrices found.')
    for i in range(1, len(mats)+1):
        print(f'{i}. {mats[i-1]}')
    choice = int(input('Please select one: '))
    try:
        mat = np.loadtxt(mats[choice])
        print(f'Loaded {mats[choice]} successfully')
    except:
        print(f'Could not load matrix from {mats[choice]}. Exiting.')
        exit(0)
elif(len(mats)==1):
    try:
        mat = np.loadtxt(mats[0])
        print(f'Loaded {mats[0]} successfully')
    except:
        print(f'Could not load matrix from {mats[0]}. Exiting.')
        exit(0)

oldSig = 1.0
w = 1.0
try:
    w = float(sys.argv[2])
except:
    None
k0 = 10

hicBonds = []
mother_ = []
serial = np.array([i for i in range(_pos.shape[0])])
motherChr = indices[serial<N][::10] + 5
motherSerial = serial[serial<N][::10] + 5
for ctr, (i, iSerial) in enumerate(zip(motherChr, motherSerial)):
    for (j, jSerial) in zip(motherChr[ctr+1:], motherSerial[ctr+1:]):
        I = int((i-1)/10)
        J = int((j-1)/10)
        
        if(mat[I, J] > 0):
            dij = oldSig/mat[I, J]
            if(w > 0):
                kij = np.exp(-((dij/oldSig-1)**2)/w)
            else:
                continue            
            if(kij > 1e-3):# and abs(i-j) > 1):
                hicBonds.append([i+1, j+1, dij, k0*kij])
                mother_.append(f'{min([i, j])}-{max([i, j])}')
                

# generating rep fork HiC bonds
if(G > 1.0):
    serial = np.array([i for i in range(_pos.shape[0])])
    repChr = indices[serial>N][np.where(indices[serial>N]%5==0)[0]][1::2]+1
    repSerial = serial[serial>N][np.where(indices[serial>N]%5==0)[0]][1::2]+1
    rep_ = []
    
    motherHiCBonds = len(hicBonds)

    for b in mother_:
        i = int(b.split('-')[0])
        j = int(b.split('-')[1])
        dij = oldSig/mat[int((i-1)/10), int((j-1)/10)]
        kij = np.exp(-((dij/oldSig-1)**2)/w)
        if(np.where(repChr==i)[0].shape[0] > 0 and np.where(repChr==j)[0].shape[0] > 0):
            hicBonds.append([repSerial[np.where(repChr==i)[0]][0]+1,
                             repSerial[np.where(repChr==j)[0]][0]+1,
                             dij, k0*kij])

hicBonds = np.array(hicBonds)
with open('hic.itp', 'w') as file:
    for ctr, (i, j, d, k) in enumerate(hicBonds):
        if(G>1):
            if(ctr == motherHiCBonds):
                file.write(';replication fork bonds\n')
        file.write(f'{int(i-1)}\t{int(j-1)}\t1\t{d:.6f}\t{k:.6f}\n')

print(f'''Repication complete. Here are the parameters and system details.
G = {G:.2f}
w = {w:.3f}
Number of beads = {_pos.shape[0]}''')
