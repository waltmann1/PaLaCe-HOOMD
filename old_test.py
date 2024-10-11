from __future__ import division
import numpy as np
from ProteinComplex import ProteinComplex
from ProteinBackbone import ProteinBackbone
from Bonds import Bonds
from Angles import Angles
from Dihedrals import Dihedrals
from Dihedrals_coupling import Dihedrals_coupling
from NonBonded import CoarseGrain
from NonBonded import Constraints
from NonBonded import HydrogenBonding
from NonBonded import Solvation
from RigidBodies import Rigid
import hoomd
from hoomd import md


def zigzag(i):
    if i % 2 == 0:
        return [1, 1, 1]
    else:
        return [0 , 0, 0]

def gentle_curve(i):

    return [i, 100 * np.sqrt(i+1), -5 * np.sqrt(i+1)]

def new_curve(i):

    theta = np.linspace(-16 * np.pi, 16 * np.pi, 1000)
    z = np.linspace(-200, 200, 1000)
    # r = z ** 2 + 1
    r = 10
    x = r * np.sin(theta)
    print(np.sin(x[i]))
    y = r * np.cos(theta)
    return [x[i], y[i], -i]

def helix(i):

    num = 1000
    theta = np.linspace(-5 * np.pi * num, 5 * np.pi * num, 20 * num)
    z = np.linspace(0, -10 * num, 20 * num)

    x = 1 * np.sin(theta)
    y = 1 * np.cos(theta)
    print(y)
    return [x[i+1], y[i+1], z[i+1]]




y = ProteinComplex('1sva.pdb')
#print(y.terminals)
#while len(y.chains) > 3:
#    y.remove_chain(y.chains[3])
#    print(len(y.chains))
#y = y.chains[3].get_subsection(133,139)
#print(len(y.chains[0].res_names) - 1)
y = y.chains[0].get_subsection(0, 204)

#quit()
#y = ProteinBackbone(4)
#y.add_all('GLU')
#print(y.res_names)
#quit()
#y.make_rigids(12, 13)
#print(y.chains[2].connections)
#print(y.chains[1].connections)
#print(y.chains[0].connections)
#y.make_rigid_in_sphere(30)
y.make_connections_rigid()
#quit()
#y.make_rigid(0,12, 0)
#y.make_rigid(0,13, 0)
#y.make_rigid(1,12, 1)
#y.make_rigid(1,13, 1)
#y.make_rigids(0, 1)
#y.dump_gsd('trajectory.gsd')
#quit()
#y = ProteinBackbone(4)
#y.add_all('PHE')

system = y.create_system()

nlist = hoomd.md.nlist.cell()
nlist.reset_exclusions(exclusions=['bond', 'angle', 'dihedral', 'constraint', 'body'])




bref = Bonds()
bref.set_all_harmonic_bonds(system)



a = Angles()
aref = a.set_all_harmonic_angles(system)



d = Dihedrals()
dref = d.set_all_palace_dihedrals(system)


dc = Dihedrals_coupling()
dcref = dc.set_all_dihedrals_coupling(system)

num = 0
r = Rigid()
r.set_rigid_bodies(system)
num = r.get_rigids(system)


#dc.reset_all_dihedrals_coupling(system, .000000001)

c = Constraints()
c.add_constraints(system, num)


hb = HydrogenBonding()
hb.set_all_hb(nlist, system)

cg = CoarseGrain()
cg.set_nb_coarse_grain(nlist, system)
cg.write_parameter_table()

sol = Solvation()
sol.add_all_solvation(system)
sol.set_all_solvation(nlist, system)

all = hoomd.group.all()

to_integrate = hoomd.group.union(name='dof', a=hoomd.group.rigid_center(), b=hoomd.group.nonrigid())

"""
fire = hoomd.md.integrate.mode_minimize_fire(dt=0.005, group=all, ftol=1e-2, Etol=1e-7)
count = 0
while not(fire.has_converged() or count > 20):
   hoomd.run(100)
   count += 1
"""
#up = hoomd.md.update.zero_momentum(period=1, phase=-1)


#fire = hoomd.md.integrate.mode_minimize_fire(dt=0.005, group=to_integrate, ftol=1e-2, Etol=1e-7)
#count = 0
#while not(fire.has_converged() or count > 20):
#   hoomd.run(100)
#   count += 1


hoomd.md.integrate.mode_standard(dt=0.005)
hoomd.md.integrate.langevin(group=to_integrate, kT=0.596, seed=42)


hoomd.analyze.log(filename="log-output.log",
                  quantities=['potential_energy', 'pair_cvcel_electrostatics_energy', 'dihedral_coupling_opls_energy'
                              , 'pair_lj_energy', 'bond_harmonic_energy', 'angle_harmonic_energy',
                              'dihedral_table_energy', 'pair_table_energy', 'pair_lj86_energy'], period=1, overwrite=True)

hoomd.dump.gsd("trajectory.gsd", period=100, group=all, overwrite=True)


#up.disable()

# run simulation
hoomd.run(500000)



