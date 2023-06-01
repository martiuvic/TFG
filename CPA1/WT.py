from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
prmtop = AmberPrmtopFile('5om9_ph7_protonated_wat.prmtop')
inpcrd = AmberInpcrdFile('5om9_ph7_protonated_wat.inpcrd')
system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
        constraints=HBonds) #PME=psrtivle mesh ewald, correcció electroestatica del periodic boundary conditions
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds) #equacio newton explicat full
simulation = Simulation(prmtop.topology, system, integrator)
simulation.context.setPositions(inpcrd.positions)
if inpcrd.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
simulation.minimizeEnergy() #
simulation.reporters.append(DCDReporter('output.dcd', 1000))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
        potentialEnergy=True, temperature=True))
simulation.step(10000) 