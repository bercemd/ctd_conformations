#This python script is used to perform replica-exchange MD simulations on GPU environment.
#Modified by Amith, W. D. and Dutagaci, B. 2023, University of California Merced, Merced CA 95348, USA.
#Based on OpenMM package and OpenMMTools package for replica-exchange molecular dynamics simulations.
#Reference:Eastman, P.; Friedrichs, M. S.; Chodera, J. D.; Radmer, R. J.; Bruns, C. M.; Ku, J. P.; Beauchamp, K. A.; Lane, T. J.; Wang, L.-P.; Shukla, D.; et al.
#OpenMM 4: A Reusable, Extensible, Hardware Independent Library for High Performance Molecular Simulation. 
#Journal of Chemical Theory and Computation 2013, 9 (1), 461-469. DOI: https://doi.org/10.1021/ct300857j  

import math
from simtk import unit
from openmmtools import testsystems, states, mcmc, multistate
from openmmtools.multistate import ReplicaExchangeSampler
from simtk.openmm import *
from simtk.openmm.app import *

from omm_readinputs import *
from omm_readparams import *
from omm_vfswitch import *
from omm_barostat import *
from omm_restraints import *
from omm_rewrap import *

import os,sys,json
import argparse
import tempfile

parser = argparse.ArgumentParser()

#Define the parameters and the paths for files which required for simulation, example: # of replicates, # of iterations, topology, pdb etc.
parser.add_argument(
    '--ff_path',type=str, default='toppar',
    help='Topology path.')
parser.add_argument(
    '--psf_path',type=str, default='protein.psf',
    help='PSF path.')
parser.add_argument(
    '--pdb_path',type=str, default='protein.pdb',
    help='PDB path.')
parser.add_argument(
    '--sysinfo',type=str, default='sysinfo.dat',
    help='Crystal info')
parser.add_argument(
    '--out_dir',type=str, default='ctdtest',
    help='Output path.')
parser.add_argument(
    '--out_file',type=str, default='storage',
    help='Output file.')
parser.add_argument(
    '--nstep',type=int, default=1000,
    help='number of steps for each iteration')
parser.add_argument(
    '--niteration',type=int, default=10,
    help='number of iterations')
parser.add_argument(
    '--nreplicate',type=int, default=3,
    help='number of replicates')
parser.add_argument(
    '--device',type=str, default='0,1,2',
    help='Device indices')

arg = parser.parse_args()

#Loading files

psf = CharmmPsfFile(arg.psf_path)
pdb = PDBFile(arg.pdb_path)
params = CharmmParameterSet('%s/top_all36mw_prot.rtf'%arg.ff_path,'%s/par_all36mw_prot.prm'%arg.ff_path,'%s/toppar_water_ions.str'%arg.ff_path)
sysinfo = json.load(open(arg.sysinfo, 'r'))
boxlx, boxly, boxlz = map(float, sysinfo['dimensions'][:3])
print(boxlx, boxly, boxlz)
psf.setBox(boxlx*angstroms, boxly*angstroms, boxlz*angstroms)
system = psf.createSystem(params,nonbondedMethod=CutoffPeriodic,nonbondedCutoff=1.2*nanometer,switchDistance=1.0*nanometer,constraints=HBonds)

# Set platform

platform = Platform.getPlatformByName('CUDA')
platform.setPropertyDefaultValue('CudaDeviceIndex', arg.device)

n_replicas = arg.nreplicate  # Number of temperature replicas.
T_min = 300.0 * unit.kelvin  # Minimum temperature.
T_max = 500.0 * unit.kelvin  # Maximum temperature.
temperatures = [T_min + (T_max - T_min) * (math.exp(float(i) / float(n_replicas-1)) - 1.0) / (math.e - 1.0)
                 for i in range(n_replicas)]
print("temperatures", temperatures)

thermodynamic_states = [states.ThermodynamicState(system=system, temperature=T)
                         for T in temperatures]
move = mcmc.LangevinSplittingDynamicsMove(timestep=2.0*unit.femtoseconds, n_steps=arg.nstep)
simulation = ReplicaExchangeSampler(mcmc_moves=move, number_of_iterations=arg.niteration)

storage_path = "%s/%s"%(arg.out_dir,arg.out_file)
reporter = multistate.MultiStateReporter(storage_path, checkpoint_interval=10)
default_box_vectors = thermodynamic_states[0].system.getDefaultPeriodicBoxVectors()
print("default_box_vectors",default_box_vectors)

sampler_states=states.SamplerState(pdb.positions)
print("sampler_states",sampler_states.box_vectors)
sampler_states.box_vectors=default_box_vectors
print("sampler_states",sampler_states.box_vectors)
simulation.create(thermodynamic_states=thermodynamic_states,
                   sampler_states=sampler_states,
                   storage=reporter)
simulation.platform = platform

print("REMD started")

simulation.run()  # This runs the simulation

print("REMD ended")



