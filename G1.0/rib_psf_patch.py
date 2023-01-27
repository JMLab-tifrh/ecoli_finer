import numpy as np
import sys
from MDAnalysis import Universe

usage = '''usage:
python rib_psf_patch.py <polymer PSF> <gro> [N polysomes] [N ribosomes] [polydispersity]
<polymer PSF>    :    PSF file containing information about the chromosome only.
<gro>            :    GROMOS 96 format file containing coordinates and atom information.
N polysomes      :    number of polysomes [default = 1538]
N monomers       :    number of monomers [default = 0]
polydispersity   :    number of monomers in each polysome [default = 13]
example:
python rib_psf_patch.py str.psf'''

if(len(sys.argv) < 3):
    print(usage)
    exit(0)

try:
    psfInp = open(sys.argv[1], 'r')
    print(f'Successfully loaded {sys.argv[1]}.')
except:
    print(f'Could not load {sys.argv[1]}. Exiting.')
    print(usage)
    exit(0)

try:
    top = Universe(sys.argv[2])
    print(f'Successfully loaded {sys.argv[2]}.')
except:
    print('Coulnd not load {sys.argv[2]}. Exiting.')
    print(usage)
    exit(0)

try:
    N_polysomes = int(sys.argv[3])
except:
    N_polysomes = 1538

try:
    N_monomers = int(sys.argv[4])
except:
    N_monomers = 0

try:
    pdi = int(sys.argv[5])
except:
    pdi = 13

print(f'Number of polysomes = {N_polysomes}')
print(f'Number of beads per polysome (polydispersity) = {pdi}')
print(f'Number of monomers = {N_monomers}')

out = open('full.psf', 'w')

N_particles = top.select_atoms('all').n_atoms
resnames = top.select_atoms('all').resnames
at_names =  top.select_atoms('all').names
resids = top.select_atoms('all').resids
resids[resids > 6] = 0

# writing atom information
out.write(f'{N_particles} !NATOM\n')
for i in range(N_particles):
    out.write("   %5d %4s %d    %3s  %4s   1     0.000000        0.0000           0\n"%(i+1,f"{resnames[i]}", resids[i],f"{at_names[i]}",f"{at_names[i]}"))

# writing bond information
line = psfInp.readline()
while('!NBOND' not in line):
    line = psfInp.readline()[:-1]
    line = line.split(' ')
    while('' in line):
        line.remove('')
chr_bonds = int(line[0])
total_bonds = chr_bonds + (pdi-1)*N_polysomes

out.write(f'\n{total_bonds} !NBOND\n')

ctr = 0
lines = psfInp.readlines()
for l in lines:
    l = l[:-1]
    pairs = l.split(' ')
    while '' in pairs:
        pairs.remove('')
    ctr += int(len(pairs)/2)
    out.write(l)
    if(ctr%4 == 0):
        out.write('\n')

N_chromosome = top.select_atoms('resname POL1').n_atoms

for i in range(N_polysomes):
    for j in range(1, pdi):
        out.write("%8d%8d"%((pdi*i)+j+N_chromosome,(pdi*i)+j+1+N_chromosome))
        ctr += 1
        if(ctr % 4 == 0):
            out.write('\n')
out.close()
