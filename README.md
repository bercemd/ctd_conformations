# Simulation and Analysis of RNA Polymerase II C-Terminal Domain Models

This repository includes an example script for replica exchange molecular dynamics (REMD) simulations and scripts for the analysis of the simulations. It also contains the most probably conformations for each system obtained by PCA analysis of the simulation trajectories in PDB format. The list of the sequences and abbreviations are provided under conformations/summary.md file. 

*** Requirements:
```
Python versions 3+

OpenMM

MDAnalysis

MDTraj
```
*** Usage:
```
Python [options] filename.py
```
*** Examples:
```
python REMD.py --ff_path=topologydir --psf_path=psffile --pdb_path=pdbfile --sysinfo=sysinfo.dat --out_dir=outputdir --out_file=outputfile --nstep=500 --niteration=10 --nreplicate=8 --device=GPUdevice
```
*** Citation
```
Amith W. D., Dutagaci, B., Complex Conformational Space of RNA Polymerase II C-Terminal Domain upon Phosphorylation, bioRxiv, doi: https://doi.org/10.1101/2023.04.20.537737, 2023.
```
