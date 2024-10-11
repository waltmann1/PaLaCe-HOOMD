from __future__ import division
from Bonds import Bonds
from Angles import Angles
from Dihedrals import Dihedrals
from Dihedrals_coupling import Dihedrals_coupling
from NonBonded import CoarseGrain
from NonBonded import Constraints
from NonBonded import HydrogenBonding
from NonBonded import Solvation
from RigidBodies import Rigid
import numpy as np
import hoomd


class Simulation(object):

    def __init__(self, system, temperature=300, name="protein_sim", geo=True):
        self.system = system
        self.nlist = hoomd.md.nlist.cell(check_period=1)
        self.nlist.reset_exclusions(exclusions=['bond', 'angle', 'dihedral', 'constraint', 'body'])
        #self.log_list = ['potential_energy', 'temperature', 'kinetic_energy']
        self.log_list = ['potential_energy', 'ndof', 'kinetic_energy']
        self.log_list.append('temperature')
        self.log_period = 1
        self.dump_period = 1000
        self.temperature = 0.596 / 300 * temperature
        self.name = name

        #self.dt = .102
        self.dt = .005
        self.time_unit = 48.9e-15
        self.num = 0
        self.rigid = Rigid()
        self.rigid.set_rigid_bodies(system)
        self.num = self.rigid.get_rigids(system)

        self.constraints = Constraints()
        self.constraints.add_constraints(self.system, self.num, geo=geo)


        self.bonds = Bonds(self.log_list)
        self.bonds.set_all_harmonic_bonds(system)

        self.angles = Angles(self.log_list)
        self.angles.set_all_harmonic_angles(system)

        self.dihedrals = Dihedrals(self.log_list)
        self.dihedrals.set_all_palace_dihedrals(system)

        self.dihedrals_coupling = Dihedrals_coupling(self.log_list)
        self.dihedrals_coupling.set_all_dihedrals_coupling(system)
        #self.dihedrals_coupling.reset_all_dihedrals_coupling(system, 1)


        self.hydrogen_bonding = HydrogenBonding(self.log_list)
        self.hydrogen_bonding.set_all_hb(self.nlist, system)

        self.coarse_grain = CoarseGrain(self.log_list)
        self.coarse_grain.set_nb_coarse_grain(self.nlist, system)
        # cg.write_parameter_table()

        self.solvation = Solvation(log_list=self.log_list)
        self.solvation.set_all_solvation(self.nlist, system)
        if geo:
            self.constraints.add_force_transfer()






        self.all = hoomd.group.all()

        self.to_integrate = hoomd.group.union(name='dof', a=hoomd.group.rigid_center(), b=hoomd.group.nonrigid())

        all = hoomd.group.all()

        h = hoomd.group.type(name='h', type='H')
        o = hoomd.group.type(name='o', type='O')
        h_o = hoomd.group.union(name='ho', a=h, b=o)

        no_h_o = hoomd.group.difference(name='no_h_o', a=self.to_integrate, b=h_o)

        if geo:
            self.to_integrate = no_h_o

        hoomd.md.integrate.mode_standard(dt=self.dt)
        self.nve = hoomd.md.integrate.nve(group=self.to_integrate, limit=.001)
        self.nve.disable()
        self.langevin = hoomd.md.integrate.langevin(group=self.to_integrate, kT=self.temperature, seed=42)
        #for type in system.particles.types:
            #self.langevin.set_gamma(str(type), 10e-3)
        #self.fire = hoomd.md.integrate.mode_minimize_fire(dt=0.005, group=self.to_integrate, ftol=1e-2, Etol=1e-7)



        log_name = self.name + ".log"
        self.logger = hoomd.analyze.log(filename=log_name, quantities=self.log_list, period=self.log_period,
                          overwrite=True)

        dump_name = self.name + ".gsd"
        self.dumper = hoomd.dump.gsd(filename=dump_name, period=self.dump_period, group=self.all, overwrite=True)

    def run(self, time):

        #print(self.system.constraints)
        hoomd.run(time)

    def run_nanoseconds(self, time):

        real_time = int(time * 1e-9 / (self.time_unit * self.dt))
        self.run(real_time)

    def nve_relaxation(self, time):


         self.langevin.disable()
         self.nve.enable()

         hoomd.run(time)
         self.nve.set_params(limit=.01)
         #hoomd.run(time / 2)
         #self.nve.set_params(limit=.001)
         #self.nve.set_params(limit=.01)
         #hoomd.run(time)
         #self.nve.set_params(limit=.1)
         #hoomd.run(time)
         #self.nve.set_params(limit=1)
         self.nve.disable()
         self.langevin.enable()

    def set_dt(self, dt):
        hoomd.md.integrate.mode_standard(dt=dt)

    def run_fire(self, time):

        self.langevin.disable()
        self.nve.enable()
        fire = hoomd.md.integrate.mode_minimize_fire(dt=0.1, group=self.to_integrate, ftol=1e-2, Etol=1e-7)
        hoomd.run(time)
        del fire
        self.langevin.enable()
        self.nve.disable()
        hoomd.md.integrate.mode_standard()

    def temp_interp(self, temp1, temp2, time):

        t1 = 0.596 / 300 * temp1
        t2 = 0.596 / 300 * temp2
        self.langevin.set_params(kT=hoomd.variant.linear_interp(points=[(0, t1), (time, t2)]))
        hoomd.run(time)
        self.langevin.set_params(kT=self.temperature)

    def set_temperature(self, t):
        temp = 0.596 / 300 * t
        self.temperature = temp
        self.langevin.set_params(kT=self.temperature)

    def palace_equil(self):
        #self.run_fire(5500)
        self.temp_interp(0, 300, 100000)
        self.run(200000)

    def basic_temp_equil_no_log(self):

        self.logger.disable()
        self.dumper.disable()
        self.set_temperature(0)
        self.run(10000)
        self.temp_interp(0, 300, 100000)
        self.set_temperature(300)
        self.run(10000)
        self.logger.enable()
        self.dumper.enable()

    def set_log_period(self, period):

        self.logger.disable()
        self.log_period = period
        log_name = self.name + ".log"
        self.logger = hoomd.analyze.log(filename=log_name, quantities=self.log_list, period=self.log_period,
                                        overwrite=True)

    def set_dump_period(self, period):

        self.dumper.disable()
        self.dump_period = period
        dump_name = self.name + ".gsd"
        self.dumper = hoomd.dump.gsd(filename=dump_name, period=self.dump_period, group=self.all, overwrite=True)


    def total_kinetic_energy(self):

        ke = 0
        for part in self.system.particles:
            kin = .5 * part.mass * np.linalg.norm(part.velocity) ** 2
            print(part.type, kin)
            ke += kin


        return ke

    def ndof(self):

        return self.total_kinetic_energy() * 2 / self.temperature


class ProteinGSD(Simulation):

    def __init__(self, name, frame=0):

        hoomd.context.initialize("")
        system = hoomd.init.read_gsd(name, frame=frame)
        super(ProteinGSD, self).__init__(system, name=name[:-4] + '_frame' + str(frame))