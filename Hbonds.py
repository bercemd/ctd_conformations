#This python script calculates the number of intrapeptide H-bonds of given trajectory.
#Modified by Amith, W. D. and Dutagaci, B. 2023, University of California Merced, Merced, CA 95343 USA.
#Based on the examples in the user guide of MDAnalysis package
#Reference:Michaud-Agrawal, N.; Denning, E. J.; Woolf, T. B.; Beckstein, O. MDAnalysis: A toolkit for the analysis of molecular dynamics simulations. 
#J.Comp.Chem 2011, 32 (10), 2319-2327. DOI: https://doi.org/10.1002/jcc.21787

import MDAnalysis
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import (
  HydrogenBondAnalysis as HBA)
import matplotlib.pyplot as plt
import pandas as pd
import sys
import numpy as np

np.set_printoptions(threshold=sys.maxsize, linewidth=100)
u = MDAnalysis.Universe("step3_input.psf", "replica_1.dcd")

hbonds = HBA(universe=u)

hbonds.hydrogens_sel = hbonds.guess_hydrogens("protein")
hbonds.acceptors_sel = hbonds.guess_acceptors("protein")

hbonds.run()

df = pd.DataFrame(hbonds.count_by_time())

pd.set_option("display.max_rows", None, "display.max_columns", None)

print("{0}".format(
      df))







