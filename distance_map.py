#This python script calculates the distance map of a given PDB file between residues of a protein.
#Modified by Amith, W. D. and Dutagaci, B. 2023, University of California Merced, Merced, CA 95343 USA.
#Based on the examples in the user guide of MDAnalysis package
#Reference:Michaud-Agrawal, N.; Denning, E. J.; Woolf, T. B.; Beckstein, O. MDAnalysis: A toolkit for the analysis of molecular dynamics simulations.
#J.Comp.Chem 2011, 32 (10), 2319-2327. DOI: https://doi.org/10.1002/jcc.21787

import MDAnalysis as mda
from MDAnalysis.tests.datafiles import PDB_small
from MDAnalysis.analysis import distances

import numpy as np
import matplotlib.pyplot as plt


u = mda.Universe("protein.pdb")

seq1_com = u.atoms.center_of_mass(compound='residues')

n_seq1 = len(seq1_com)

print('sequence1 has {} residues'.format(n_seq1))
res_dist = distances.self_distance_array(seq1_com)
res_dist.shape

sq_dist_res = np.zeros((n_seq1, n_seq1))
triu = np.triu_indices_from(sq_dist_res, k=1)
sq_dist_res[triu] = res_dist
sq_dist_res.T[triu] = res_dist

vmin = 0
vmax = 6
 
fig2, ax2 = plt.subplots()
im2 = ax2.pcolormesh(u.residues.resids, u.residues.resids, sq_dist_res, shading='auto', vmin=vmin, vmax=vmax)
 
# plt.pcolor gives a rectangular grid by default
# so we need to make our heatmap square
ax2.set_aspect('equal')
 
# add figure labels and titles
plt.ylabel('Residue IDs')
plt.xlabel('Residue IDs')
 
# colorbar
cbar2 = fig2.colorbar(im2)
cbar2.ax.set_ylabel('Distance (Angstrom)')
plt.show()
