import numpy as np
from MDAnalysis import Universe
from MDAnalysis.analysis.distances import distance_array as get_dist
import sys
from groio import write_gro
from glob import glob
import sys

usage = '''python gentopol.py <gro file> <w> <PBC>[default:0] <k0>[default:10] <HiC matrix>[default: search in PWD]
Example:
python gentopol.py pol1.gro 0.2 1'''

if(len(sys.argv) < 2):
    print(usage)
    exit(0)

try:
    top = Universe(sys.argv[1], sys.argv[1])
except:
    print(f'ERROR: Could not load {sys.argv[1]}.')
    exit(0)

try:
    w = float(sys.argv[2])
except:
    print("ERROR: Could not read w value.")
    exit(0)

try:
    pbc = int(float(sys.argv[3]))
except:
    pbc = 0

try:
    k0 = float(sys.argv[4])
except:
    k0 = 10.0

try:
    mat = np.loadtxt(sys.argv[5])
    print(f'Loaded {sys.argv[5]} successfully.')
except:
    if(len(sys.argv) > 5):
       print(f'Could not load {sys.argv[5]}. Autodetection ON.')
    mats = glob("expt*.mat")
    if(len(mats) > 1):
        ch = 0
        print('Multiple experimental matrices found. Please select one.')
        print('0.\tExit')
        for i in range(len(mats)):
            print(f'{i+1}.\t{mats[i]}')
        ch = int(input('Enter your choice: '))
        if(ch == 0):
            print('Exiting...')
            exit(0)
        elif(ch > len(mats)):
            print('ERROR: Incorrect choice. Exiting.')
            exit(0)
    else:
        ch = 0

    try:
        mat = np.loadtxt(mats[0])
        print(f'Loaded {mats[ch]} successfully.')
    except:
        print(f'ERROR: Could not load {mats[ch]}. Exiting.')

try:
    labels = np.loadtxt('labels.txt')
except:
    print('ERROR: Could not load labels.txt. Exiting.')
    exit(0)

try:
    nBB = np.loadtxt('bonds.txt').shape[0]
except:
    print('ERROR: Could not load bonds.txt. Exiting.')
    exit(0)

if(pbc):
    print('Doubling box size and placing system at the center to enable PBC.')

## some parameters
nbeads = top.select_atoms('all').n_atoms
pos = top.select_atoms('all').positions*0.1
box = top.trajectory[0]._unitcell[:3]*0.1
newpos = []
res = 0.5 #kbp
newRes = int(5/res) 

sig = 0.464
factor = (res/newRes)**(1.0/3.0)
oldSig = 1.0

## generating topology and PSF
indices = np.linspace(1, nbeads, nbeads).astype('int')
topname = 'topol.top' #f"{(sys.argv[1]).split('.')[0]}" + '.top'
top = open(topname, 'w')
pol = open('pol.itp', 'w')

psf = open('str.psf', 'w')
psf.write("%d !NATOM\n"%(pos.shape[0]))

top.write("#include \"param.itp\"\n#include \"pol.itp\"\n#include \"polysome.itp\"\n")
top.write('#include "monomer30S.itp"\n#include "monomer50S.itp"\n\n')
pol.write("""[ moleculetype ]
POL1    0\n\n
[ atoms ]\n""")

resnames = Universe(sys.argv[1]).select_atoms('all').names
resids = Universe(sys.argv[1]).select_atoms('all').resids

for i in range(0, pos.shape[0]):
    pol.write(f'{i+1}\t{resnames[i]}\t{resids[i]}\tPOL1\t{resnames[i]}\t{i+1}\t{0:.4f}\n')
    psf.write("   %5d %4s 1    %3s  %4s   1     0.000000        0.0000           0\n"%(i+1,f"POL1",f"{resnames[i]}",f"{resnames[i]}"))

G = float(open(sys.argv[1], 'r').readline().split(' ')[-1][2:-1])
if(G > 1):
    I, J = np.loadtxt('bonds_fork.txt', dtype='int', comments='#', unpack=True)
else:
    I, J = np.loadtxt('bonds.txt', dtype='int', comments='#', unpack=True)

psf.write("%d !NBOND\n"%(int(I.shape[0])))
ctr = 0
pol.write('\n[ bonds ]\n')
# backbone + plectoneme bonds
k_adj = 300.0
for i, j in zip(I, J):
    pol.write(f'{i}\t{j}\t1\t{sig:.6f}\t{k_adj:.6f}\n')
    ctr+=1
    psf.write("%8d%8d"%(i,j))
    if ctr%4==0:
        psf.write("\n")

hicBeads = indices[indices<int(4.64e6/500)][::10] + 5

# HiC bonds
pol.write('\n#ifdef HIC\n')
pol.write(f"; HiC bonds at w={w:.2f}, k0={k0:.2f}\n")
pol.write('#include "hic.itp"\n')
pol.write('#endif\n')

# Position restraints
pol.write('\n#ifdef POSRES\n')
pol.write('; spherocylindrical position restraints\n')
pol.write('[ position_restraints ]\n')
for i in range(1, pos.shape[0]+1):
    if(pbc):
        pol.write("%d\t%d\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n" %
                (i, 2, 6, 0.5*box[0], 600.0, box[2]-box[0], box[0], box[1], box[2]))
    else:
        pol.write("%d\t%d\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n" %
                (i, 2, 6, 0.5*box[0], 600.0, box[2]-box[0], 0.5*box[0], 0.5*box[1], 0.5*box[2]))

pol.write('#endif\n')

first_line = open(sys.argv[1], 'r').readline()

# defining molecules
top.write(f'''\n[ system ]
{first_line}

[ molecules ]
POL1    1
;R70S    1292
;R30S    2102
;R50S    2102''')
top.close()
print(f'Succesfully written topology to {topname}.')
