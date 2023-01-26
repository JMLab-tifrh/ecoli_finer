from math import sqrt
import sys

## Masses
m0   = 649*500 # Da
m30S = 0.85e6/m0
m50S = 1.36e6/m0
m70S = 2.31e6/m0
m0 = 1

## Sizes
sig0   = 6.821e-8 # m
sig30S = 14e-9/sig0
sig50S = 17e-9/sig0
sig70S = 21e-9/sig0
sig = (0.1**(1/3.0))


## Epsilons
ep0 = 1
ep_rib = 1
ep6 = 0
ep12 = 1
try:
    ep6, ep12 = float(sys.argv[1]),  float(sys.argv[2])
except:
    print(f'Using default values ep6={ep6:.2f},  ep12={ep12:.2f}. Pass an argument to change it.')

out = open('param.itp', 'w')
out.write(f'''[ defaults ]
1       1       no      1       1

[ atomtypes ]
;atomtype       mass    charge  particle        c6      c12
;-------        ----    -----   -------         ---     ---
F               1.00    0.00    A               0.0     {4*ep0*(sig**12):.3e}
T               1.00    0.00    A               0.0     {4*ep0*(sig**12):.3e}
B               1.00    0.00    A               0.0     {4*ep0*(sig**12):.3e}
R3              {m30S:.2f}    0.00    A               0.0     {4*ep_rib*(sig30S**12):.3e}
R5              {m50S:.2f}    0.00    A               0.0     {4*ep_rib*(sig50S**12):.3e}
R7              {m70S:.2f}    0.00    A               0.0     {4*ep_rib*(sig70S**12):.3e}
[ nonbond_params ]
;i      j       func    c6      c12
;--     --      ----    --      ---
;intra chromosome
F       F       1       0.00    {4*ep0*(sig**12):.3e}
T       T       1       0.00    {4*ep0*(sig**12):.3e}
B       B       1       0.00    {4*ep0*(sig**12):.3e}
F       T       1       0.00    {4*ep0*(sig**12):.3e}
F       B       1       0.00    {4*ep0*(sig**12):.3e}
B       T       1       0.00    {4*ep0*(sig**12):.3e}
;ribosomal interactions
R3      R3      1       0.00    {4*ep_rib*(sig30S**12):.3e}
R5      R5      1       0.00    {4*ep_rib*(sig50S**12):.3e}
R7      R7      1       0.00    {4*ep_rib*(sig70S**12):.3e}
R3      R5      1       0.00    {sqrt(4*ep_rib*(sig30S**12)*4*ep_rib*(sig50S**12)):.3e}
R5      R7      1       0.00    {sqrt(4*ep_rib*(sig50S**12)*4*ep_rib*(sig70S**12)):.3e}
R7      R3      1       0.00    {sqrt(4*ep_rib*(sig70S**12)*4*ep_rib*(sig30S**12)):.3e}
;chromosome - ribosome interactions 
F       R3      1       {sqrt(4*ep6*(sig30S**6)*4*ep6*(sig**6)):.3e}    {sqrt(4*ep12*(sig30S**12)*4*ep12*(sig**12)):.3e} 
T       R3      1       {sqrt(4*ep6*(sig30S**6)*4*ep6*(sig**6)):.3e}    {sqrt(4*ep12*(sig30S**12)*4*ep12*(sig**12)):.3e}
B       R3      1       {sqrt(4*ep6*(sig30S**6)*4*ep6*(sig**6)):.3e}    {sqrt(4*ep12*(sig30S**12)*4*ep12*(sig**12)):.3e}
F       R5      1       {sqrt(4*ep6*(sig30S**6)*4*ep6*(sig**6)):.3e}    {sqrt(4*ep12*(sig50S**12)*4*ep12*(sig**12)):.3e}
T       R5      1       {sqrt(4*ep6*(sig30S**6)*4*ep6*(sig**6)):.3e}    {sqrt(4*ep12*(sig50S**12)*4*ep12*(sig**12)):.3e}
B       R5      1       {sqrt(4*ep6*(sig30S**6)*4*ep6*(sig**6)):.3e}    {sqrt(4*ep12*(sig50S**12)*4*ep12*(sig**12)):.3e}
F       R7      1       {sqrt(4*ep6*(sig30S**6)*4*ep6*(sig**6)):.3e}    {sqrt(4*ep12*(sig70S**12)*4*ep12*(sig**12)):.3e}
T       R7      1       {sqrt(4*ep6*(sig30S**6)*4*ep6*(sig**6)):.3e}    {sqrt(4*ep12*(sig70S**12)*4*ep12*(sig**12)):.3e}
B       R7      1       {sqrt(4*ep6*(sig30S**6)*4*ep6*(sig**6)):.3e}    {sqrt(4*ep12*(sig70S**12)*4*ep12*(sig**12)):.3e}
 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                     actual values                      ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;ep6 = {ep6:.2e}, ep12 = {ep12:.2e}
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;sigma = {sig0} (length scales)
;chromosome(1 bead = 500bp)
;size of 1 chromosome bead = {sig*sig0*1e9:.2f} nm
;mass of 1 chromosome bead = {649*500/1e3:.2e} KDa [1]
;epsilon for all chromosome beads = {ep0:.2f}
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;ribosome parameters
;size of 30S subunit = 14 nm [2]
;mass of 30S subunit = {0.85e3:.2e} KDa [3]
;epsilon of 30S subunit = {ep_rib:.2e}
;
;size of 50S subunit = 17 nm [2]
;mass of 50S subunit = {1.39e3:.2e} KDa [3]
;epsilon of 50S subunit = {ep_rib:.2e}
;
;size of 70S subunit = 21 nm [2]
;mass of 70S subunit = {2.31e3:.2e} KDa [3]
;epsilon of 70S subunit = {ep_rib:.2e}
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                         References                     ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;[1] https://bionumbers.hms.harvard.edu/files/Nucleic%20Acids_Sizes_and_Molecular_Weights_2pgs.pdf
;[2] https://onlinelibrary.wiley.com/doi/full/10.1111/mmi.12805
;[3] https://www.nature.com/articles/nmeth.4147/
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
''')

out.close()
