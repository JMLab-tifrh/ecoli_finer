# _E. coli_ chromosome model at 500 base pair resolution
This repository contains scripts and code for generating and simulating _E. coli_ chromosome at 500 bp inside a cell with ribosomes being present as the cytoplasmic particles.<br>

This repository contains two folders named ```G1.0``` and ```G1.6```.<br>
  * ```G1.0``` contains files for simulating a single, unreplicated chromosome in an _E. coli_ cell.<br>
  * ```G1.6``` contains files for simulating a partially replicated chromosome in an _E. coli_ cell.<br>
  
## List and descriptions of files required for generation and simulation of chromosome.
  * **RNA-Seq data** : The RNA-Seq data in present in a file named ```rnaByDNA_scholz.csv```[[1]](https://www.sciencedirect.com/science/article/pii/S2405471219300389).
  * **expt_pij_wt30MM_5kb.mat** : The Hi-C contact probability matrix at 5000 bp resolution for WT cells grown at 30 &deg;C in M9 minimal media[[2]](https://www.sciencedirect.com/science/article/pii/S0092867417315076).
  * **box.dat** : This is an optiional file. If present, ```build.py``` will detect it and use it to define the spherocylinder and simulation box dimensions. Else the script will ask for explicit input from the user. The values are in reduced units. $\sigma = 68.21$ nm, where $\sigma$ is the value of unit length in our simulations. Therefore the cell dimensions, radius $(r) = 0.82$ $\mu m$ and total length $(L) = 3.05$ $\mu m$, reduce to $r = 12.314910 \sigma$ as the $L = 44.714850 \sigma$.<br>

## General protocol
The gene


** Documentation is being compiled.

## References
[[1]](https://www.sciencedirect.com/science/article/pii/S2405471219300389)  Scholz, S.A., Diao, R., Wolfe, M.B., Fivenson, E.M., Lin, X.N. and Freddolino, P.L., 2019. High-resolution mapping of the Escherichia coli chromosome reveals positions of high and low transcription. Cell systems, 8(3), pp.212-225.<br>
[[2]](https://www.sciencedirect.com/science/article/pii/S0092867417315076) Lioy, V.S., Cournac, A., Marbouty, M., Duigou, S., Mozziconacci, J., Espéli, O., Boccard, F. and Koszul, R., 2018. Multiscale structuring of the E. coli chromosome by nucleoid-associated and condensin proteins. Cell, 172(4), pp.771-783.
