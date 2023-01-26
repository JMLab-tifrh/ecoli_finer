# this script automates the simulation protocol
# these are the steps required for simulation

# defining executable paths. replace with your local paths to executables
python="python"
gmx="/apps/gmx-5.0.6_sph/bin/gmx"
gmx22="/apps/gmx-2022/bin/gmx"

# 1. classification of chromosome regions into PARs and PFRs
echo "========================================================"
echo "classifying chromosome into PARs and PFRs"
#$python rnaSeq_to_labels.py
echo "========================================================"
# 2. generation of the hyper-branched polymer topology
echo "========================================================"
echo "generating polymer topology"
#$python build.py labels_rnaSeq.txt 1.0 # just run python build.py for usage.
echo "========================================================"

# 3. generation of force-field files
echo "========================================================"
echo "generating force field files"
#$python gentopol.py pol.gro 1.0 
echo "========================================================"

# 4. confining generated initial configuration in the capsule
# for this we will use nrg_min.mdp without Hi-C bonds.
echo "========================================================"
#$gmx grompp -f nrg_min.mdp -c pol.gro -p topol.top -o em1.tpr
#$gmx mdrun -v -deffnm em1 -nt 1
echo "========================================================"

# 5. inserting ribosomes and polysomes. GROMACS versions 2019 and above are required.
echo "========================================================"
$gmx22 insert-molecules -f em1.gro -ci polysome.gro -nmol 1292 -o temp1.gro
$gmx22 insert-molecules -f temp1.gro -ci monomer30S.gro -nmol 2102 -o temp2.gro
$gmx22 insert-molecules -f temp2.gro -ci monomer50S.gro -nmol 2102 -o conf.gro
echo "========================================================"

# 6. energy minimization to remove any overlaps and confine
# any ribosome that falls outside the capsule
echo "========================================================"
$gmx grompp -f nrg_min.mdp -c conf.gro -p topol.top -o em2.tpr
$gmx mdrun -v -deffnm em2 -nt 1
echo "========================================================"

# 7. turn on Hi-C bonds and perform another minimization
echo "========================================================"
$gmx grompp -f nrg_min.mdp -c em2.gro -p topol.top -o em3.tpr
$gmx mdrun -v -deffnm em3 -nt 1
echo "========================================================"

# 8. perform a short equiliration run
echo "========================================================"
$gmx grompp -f prod_md_init.mdp -c em3.gro -p topol.top -o md.tpr -maxwarn 1
$gmx mdrun -v -deffnm md -nt 1 -nsteps 5000 #100000
echo "========================================================"

# 9. center system at the centre of an enlarged box so that PBC
# can be turned on without the particles ever crossing the PBC.
# Turning on PBC allows us to use parallelization of the simulation
# hence a huge boost in the simulation speed. This is only a techicality
# if simulations are being performed on the modified GROMACS we use
echo "========================================================"
$gmx editconf -f md.gro -o start.gro -box 24.62982 24.62982 89.4297 -c
echo "========================================================"

# 10. regenerate the topology so as to incorporate the PBC
# and adjust the position restrain centres. Again this is specific to 
# the GROMACS we are using
echo "========================================================"
$python gentopol.py pol.gro 1.0 1.0
echo "========================================================"

# 11. remove ';' from the last three lines of topol.top. ';' means commented out.
echo "========================================================"
sed -i 's/;//g' topol.top
echo "========================================================"

# 12. perform another minimization so that no particle crosses PBC
echo "========================================================"
$gmx grompp -f nrg_min_gpu.mdp -c start.gro -p topol.top -o em.tpr -maxwarn 1
$gmx mdrun -v -deffnm em -nt 1 -nb gpu
echo "========================================================"

# 13. perform molecular dynamics simulations
echo "========================================================"
$gmx grompp -f prod_md.mdp -c em.gro -p topol.top -o run.tpr -maxwarn 1
$gmx mdrun -v -deffnm run -nt 1 -nb gpu -nsteps 10000
echo "========================================================"

# 14. re-centre the system in a box of the size of the capsule.
echo "========================================================"
echo 0 0|$gmx trjconv -f run.xtc -s run.tpr -o cen.xtc -box 12.31491 12.31491 44.71485 -center
echo 0 0|$gmx trjconv -f run.gro -s run.tpr -o cen.gro -box 12.31491 12.31491 44.71485 -center
echo "========================================================"

rm *#
