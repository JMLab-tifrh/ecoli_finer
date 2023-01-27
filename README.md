# _E. coli_ chromosome model at 500 base pair resolution
This repository contains scripts and code for generating and simulating _E. coli_ chromosome at 500 bp inside a cell with ribosomes being present as the cytoplasmic particles.<br>

This repository contains two folders named ```G1.0``` and ```G1.6```.<br>
  * ```G1.0``` contains files for simulating a single, unreplicated chromosome in an _E. coli_ cell.<br>
  * ```G1.6``` contains files for simulating a partially replicated chromosome in an _E. coli_ cell.<br>
  
## List and descriptions of files required for generation and simulation of chromosome.
**GROMACS** has been abbreviated as **GMX** from now on.<br><br>

  * **RNA-Seq data** : The RNA-Seq data in present in a file named ```rnaByDNA_scholz.csv```[[1]](https://www.sciencedirect.com/science/article/pii/S2405471219300389).
  * **expt_pij_wt30MM_5kb.mat** : The Hi-C contact probability matrix at 5000 bp resolution for WT cells grown at 30 &deg;C in M9 minimal media[[2]](https://www.sciencedirect.com/science/article/pii/S0092867417315076).
  * **box.dat** : This is an optiional file. If present, ```build.py``` will detect it and use it to define the spherocylinder and simulation box dimensions. Else the script will ask for explicit input from the user. The values are in reduced units. $\sigma = 68.21$ nm, where $\sigma$ is the value of unit length in our simulations. Therefore the cell dimensions, radius $(r) = 0.82$ $\mu m$ and total length $(L) = 3.05$ $\mu m$, reduce to $r = 12.314910 \sigma$ as the $L = 44.714850 \sigma$.<br>
  * **pol.gro** : polymer conformation for the hyper-branched circular polymer genetared by ```build.py```.
  * **labels.txt** : bead labels (0 --> PFR, -1 --> plectoneme bead from which a plectoneme originates, 1 --> plectoneme, 2 -->hyper-branch bead from which the hyper-branch originates, 3 --> hyper-branch) generated by ```build.py```.
  * **bonds.txt** : bond connections required to maintain the determined hyper-branched architecture. generated by ```build.py``.
  * **index.txt** : this file maps each bead to its parent chromosome region. This ensures that no matter in what order the polymer beads are present in ```pol.gro```, they can always be mapped onto the genome.
  * **pol_fork.gro** : generated by ```genRepfork.py``` by reading ```pol.gro```. This script determines which region has been replicated by reading the header of ```pol.gro```. The header contains the information about the replication status of the chromosome. It also updates ```bonds.txt```, ```labels.txt``` and ```index.txt```.
  * **labels_fork.txt** : updated labels for the polymer where the labels of the replicated regions have been appended.
  * **bonds_fork.txt** : updated labels for the polymer where the information regarding the bonds of the replicated regions have been appended.
  * **indices_fork.txt** : updated indices such that the replicated region of the chromosome can now be mapped to the parent region, from where it was replicated.
  
  
  #### files for ribosome parameters and topology
  * **pol.itp** : needed by GMX. contains the architecture (bonds, particle numbers, system definition) for the polymer.
  * **hic.itp** : needed by GMX. contains information on _Hi-C_ bonds.
 
  #### files for ribosome parameters and topology
  * **monomer30S.itp** : needed by GMX. topology for a single 30S ribosomal subunit.
  * **monomer50S.itp** : needed by GMX. topology for a single 50S ribosomal subunit.
  * **polysome.itp** : needed by GMX. topology for a 13-mer of 7S ribosomal subunits.
  * **monomer30S.gro** : needed by GMX. contains initial configuration for a single 30S ribosomal subunit.
  * **monomer50S.gro** : needed by GMX. contains initial configuration for a single 50S ribosomal subunit.
  * **polymer.gro** : needed by GMX. contains initial configuration for a 13-mer of 70S ribosomal subunits.

  #### files defining the complete system topology and all non-bonded interactions
  * **topol.top** : needed by GMX. defines the complete system topology.
  * **param.itp** : needed by GMX. contains information regarding non-bonded interactions among various bead types.

  #### simulation configuration files (.mdp)
  * **nrg_min.mdp** : energy minimization file when PBC has been turned off. Cannot use paralleization of simulation.
  * **nrg_min_gpu.mdp** : energy minimization file when PBC is present. Can use paralleization and GPU acceleration, if present.
  * **prod_md_init.mdp** : MD configuration file when PBC has been turned off. Cannot use paralleization of simulation.
  * **prod_md.mdp** : MD configuration file when PBC is present. Can use paralleization and GPU acceleration, if present.

## Files required for proper visualization of trajectories and conformation
* **str.psf** : PSF for chromosome architecture only. generated by ```gentopol.py```
* **full.psf** : PSF containing information regarding both chromosome and ribosomes. generated by running ```rib_patch_psf.py```
  

## General protocol and usage of each script
The following is the generic protocol for any simulation.<br>
1. use ```rnaSeq_to_labels.py``` to convert RNA-Seq data to PARs and PFRs information. The file generated from this step is ```labels_rnaSeq.txt```. This step is required only once. In case one needs to re-run simulations, one may start from step-2 directly if ```labels_rnaSeq.txt``` is already present. The file ```rnaByDNA_scholz.csv``` is required in this step.<br>
2. use ```build.py``` to generate the initial architecture of the chromosome. Here the files required are ```labels_rnaSeq.txt```, ```expt_pij_wt30MM_5kb.mat``` and optionally ```box.dat```. ```build.py``` will require two arguments that need to be passed during execution. The first is the path to the file ``labels_rnaSeq.txt``` and secondly the _G value_ of the chromosome, i.e. its replication status. Upon successful execution, the following files will be generated: ```pol.gro```Though this script will not generate the achitecture of the replicating chromosome, it will write the information in the first line of the output file ```pol.gro```. Another script ```genRepfork.py``` is then required to be run for generation of a replication fork. ```genRepfork.py``` will read the first line of ```pol.gro``` to determine the replication status and generate another file called 
3. 


** Documentation is being compiled.

## References
[[1]](https://www.sciencedirect.com/science/article/pii/S2405471219300389)  Scholz, S.A., Diao, R., Wolfe, M.B., Fivenson, E.M., Lin, X.N. and Freddolino, P.L., 2019. High-resolution mapping of the Escherichia coli chromosome reveals positions of high and low transcription. Cell systems, 8(3), pp.212-225.<br>
[[2]](https://www.sciencedirect.com/science/article/pii/S0092867417315076) Lioy, V.S., Cournac, A., Marbouty, M., Duigou, S., Mozziconacci, J., Espéli, O., Boccard, F. and Koszul, R., 2018. Multiscale structuring of the E. coli chromosome by nucleoid-associated and condensin proteins. Cell, 172(4), pp.771-783.
