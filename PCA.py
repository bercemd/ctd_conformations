#This python script calculates the first and second principal components based on the cartesian coordinates of backbone atoms
# of a given trajectory.
#Modified by Amith, W. D. and Dutagaci, B. 2023, University of California Merced, Merced, CA 95343 USA.
#Based on the examples in the user guide of MDAnalysis package.
#Reference:Michaud-Agrawal, N.; Denning, E. J.; Woolf, T. B.; Beckstein, O. MDAnalysis: A toolkit for the analysis of molecular dynamics simulations.
#J.Comp.Chem 2011, 32 (10), 2319-2327. DOI: https://doi.org/10.1002/jcc.21787

import MDAnalysis as mda
import MDAnalysis.analysis.pca as pca
import MDAnalysis.analysis.align as align
from MDAnalysis.tests.datafiles import PSF, DCD, PDB
import numpy as np
import pandas as pd

u = mda.Universe("step1_pdbreader.psf", "protein_1.dcd", in_memory=True)

protein = u.select_atoms("protein")


init_frame = mda.Universe("initial-frame.pdb")
prealigner = align.AlignTraj(u, init_frame, select="protein", in_memory=True).run()

#builing average structure
reference_coordinates = u.trajectory.timeseries(asel=protein).mean(axis=1)
reference = mda.Merge(protein).load_new(reference_coordinates[:, None, :], order="afc")


aligner = align.AlignTraj(u, reference, select='backbone',
                          in_memory=True).run()


#for ts in u.trajectory:

pc = pca.PCA(u, select='backbone',
            align=True, mean=None,
            n_components=2).run()

backbone = u.select_atoms('backbone')
transformed = pc.transform(backbone, n_components=2)
df = pd.DataFrame(transformed,
         columns=['PC{}'.format(i+1) for i in range(2)])

pd.set_option("display.max_rows", None, "display.max_columns", None)

print("{0}".format(
      df))

