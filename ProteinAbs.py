from __future__ import division

import hoomd
import numpy as np
import numpy.linalg as la
from Quaternion import QuaternionBetween
from ProteinCollection import ProteinCollection


class ProteinAbs(object):

    def __init__(self, residues, chain_index=0):
        """

        :param residues: the number of residues
        """

        self.residues = residues
        self.num_particles = 0

        self.name = 'name'

        # particle types
        self.type = []

        # particle positions
        self.position = []

        # particle masses in amus
        self.mass = []

        self.charge = []

        self.body = []

        self.bonds = []
        self.bond_names = []

        self.angles = []
        self.angle_names = []

        self.dihedrals = []
        self.dihedral_names = []

        self.dihedrals_coupling = []
        self.dihedral_coupling_names = []

        self.used_alpha_carbons = []
        self.active_list = []

        self.chain_index = chain_index
        self.rigid_centers = []
        self.complexes = []


    def chain_vector(self):
        """Computes the length of the chain

        :return: vector from first to last position
        :rtype: numpy pos array [x,y,z]
        """

        return np.subtract(self.position[np.sum(self.num_particles) - 1], self.position[0])

    def dump_xyz(self, dump_file=None):
        """

        :param dump_file: name of the file to dump the xyz to
        :return: nothing
        """
        filename = dump_file
        if dump_file is None:
            filename = self.name + '.xyz'
        elif dump_file[-4:] != '.xyz':
            filename = dump_file + '.xyz'
        f_out = open(filename, 'w')
        f_out.write(str(np.sum(self.num_particles)))
        f_out.write('\n\n')
        for i in range(np.sum(self.num_particles)):
            f_out.write('%s %f %f %f\n' % (self.type[i], self.position[i][0], self.position[i][1],
                                           self.position[i][2]))
        f_out.close()

    def create_system(self):
        """

        :return: system object
        """
        hoomd.context.initialize("")
        a_types = self.angle_names
        b_types = self.bond_names
        d_types = self.dihedral_names
        p_types = list(set(self.type)) + ["center" + str(i) for i in range(len(self.rigid_centers))]
        dc_types = self.dihedral_coupling_names
        snap = hoomd.data.make_snapshot(self.num_particles + len(self.rigid_centers), particle_types=p_types,
                                        dihedral_types=d_types, bond_types=b_types, angle_types=a_types,
                                        dihedral_coupling_types=dc_types,
                                        box=hoomd.data.boxdim(L=int(la.norm(self.max_radius())) * 2 + 20))

        cen = self.geometric_center()
        self.shift([c * -1 for c in cen])
        snap.bonds.resize(0)

        for x in range(len(self.rigid_centers)):
            snap.particles.position[x] = self.center_of_mass(x)
            snap.particles.mass[x] = np.sum([self.mass[i] for i in range(self.num_particles) if self.body[i] == x])
            snap.particles.typeid[x] = p_types.index("center" + str(x))
            snap.particles.body[x] = x
            snap.particles.moment_inertia[x] = self.moment_inertia(x)
        tag = len(self.rigid_centers)

        for x in range(len(self.get_bond_types())):
            bodies = [self.body[i] for i in self.get_bonds()[x] if self.body[i] != -1]
            if len(bodies) == len(set(bodies)):
                bond_number = snap.bonds.N + 1
                snap.bonds.resize(bond_number)
                snap.bonds.group[bond_number - 1] = np.add(self.get_bonds()[x], tag)
                snap.bonds.typeid[bond_number - 1] = b_types.index(self.get_bond_types()[x])
        for x in range(len(self.get_angle_types())):
            bodies = [self.body[i] for i in self.get_angles()[x] if self.body[i] != -1]
            if len(bodies) == len(set(bodies)):
                angle_number = snap.angles.N + 1
                snap.angles.resize(angle_number)
                snap.angles.group[angle_number - 1] = np.add(self.get_angles()[x], tag)
                snap.angles.typeid[angle_number - 1] = a_types.index(self.get_angle_types()[x])
        for x in range(len(self.get_dihedral_types())):
            bodies = [self.body[i] for i in self.get_dihedrals()[x] if self.body[i] != -1]
            if len(bodies) == len(set(bodies)):
                dihedral_number = snap.dihedrals.N + 1
                snap.dihedrals.resize(dihedral_number)
                snap.dihedrals.group[dihedral_number - 1] = np.add(self.get_dihedrals()[x], tag)
                snap.dihedrals.typeid[dihedral_number - 1] = d_types.index(self.get_dihedral_types()[x])
        for x in range(len(self.get_dihedral_coupling_types())):
            bodies = [self.body[i] for i in set(self.get_dihedrals_coupling()[x]) if self.body[i] != -1]
            if len(bodies) == len(set(bodies)):
                dihedral_coupling_number = snap.dihedrals_coupling.N + 1
                snap.dihedrals_coupling.resize(dihedral_coupling_number)
                snap.dihedrals_coupling.group[dihedral_coupling_number - 1] = \
                    np.add(self.get_dihedrals_coupling()[x], tag)
                snap.dihedrals_coupling.typeid[dihedral_coupling_number - 1] = \
                    dc_types.index(self.get_dihedral_coupling_types()[x])

        for x in range(np.sum(self.num_particles)):
            snap.particles.position[x + tag] = self.position[x]
            snap.particles.mass[x + tag] = self.mass[x]
            snap.particles.typeid[x + tag] = p_types.index(self.type[x])
            snap.particles.body[x + tag] = self.body[x]
        sys = hoomd.init.read_snapshot(snap)

        return sys

    def dump_gsd(self, dump_file=None):
        """

        :param dump_file: name of the file to dump the xyz to
        :return: nothing
        """

        filename = dump_file
        if dump_file is None:
            filename = self.name + '.gsd'
        elif dump_file[-4:] != '.gsd':
            filename = dump_file + '.gsd'

        sys = self.create_system()

        hoomd.dump.gsd(filename=filename, period=None, group=hoomd.group.all(), overwrite=True)
        return filename

    def dump_getar(self, dump_file=None):
        """

        :param dump_file: name of the file to dump the xyz to
        :return: nothing
        """

        filename = dump_file
        if dump_file is None:
            filename = self.name + '.getar'
        elif dump_file[-4:] != '.getar':
            filename = dump_file + '.getar'

        sys = self.create_system()
        # dump the getar
        print("you don't understand how to dump a getar, remember?")
        return filename

    def get_bonds(self):
        return self.bonds

    def get_angles(self):
        return self.angles

    def get_dihedrals(self):
        return self.dihedrals

    def get_dihedrals_coupling(self):
        return self.dihedrals_coupling

    def get_bond_types(self):
        return self.bond_names

    def get_angle_types(self):
        return self.angle_names

    def get_dihedral_types(self):
        return self.dihedral_names

    def get_dihedral_coupling_types(self):
        return self.dihedral_coupling_names

    def shift(self, vec):
        """

        :param vec: shift vector
        :return:
        """

        self.position = np.add(self.position, vec)

    def align(self, vec):

        q = QuaternionBetween(self.chain_vector(), vec)
        for x in range(len(self.position)):
            self.position[x] = q.orient(self.position[x])

    def geometric_center(self):
        """

        :return: position of the geometric center of the protein
        """

        return np.mean(self.position, axis=0)

    def max_radius(self):

        t = [-1 * c for c in self.geometric_center()]
        self.shift(t)
        rs = [la.norm(pos) for pos in self.position]
        self.shift([-1 * g for g in t])
        return np.max(rs)

    def moment_inertia(self, center_index):
        mass_array = np.array([self.mass[i] for i in range(self.num_particles) if self.body[i] == center_index])
        position_array = np.array([self.position[i] for i in range(self.num_particles) if self.body[i] == center_index])

        cen = self.center_of_mass_arrays(position_array, mass_array)
        position_array = np.subtract(position_array, cen)
        return self.sum_over_xyz(self.pos_x_mass(self.pos_squared(position_array), mass_array))

    def center_of_mass(self, center_index):
        mass_array = np.array([self.mass[i] for i in range(self.num_particles) if self.body[i] == center_index])
        position_array = np.array([self.position[i] for i in range(self.num_particles) if self.body[i] == center_index])

        return np.sum(self.pos_x_mass(position_array, mass_array), axis=0) / np.sum(mass_array)

    def center_of_mass_arrays(self, position_array, mass_array):
        return np.sum(self.pos_x_mass(position_array, mass_array), axis=0) / np.sum(mass_array)

    def pos_x_mass(self, position_array, mass_array):

        y = np.zeros_like(position_array)
        for ind in range(len(mass_array)):
            y[ind][0] = position_array[ind][0] * mass_array[ind]
            y[ind][1] = position_array[ind][1] * mass_array[ind]
            y[ind][2] = position_array[ind][2] * mass_array[ind]
        return y

    def pos_squared(self, position_array):

        y = np.zeros_like(position_array)

        for ind in range(len(position_array)):
            y[ind][0] = position_array[ind][0] * position_array[ind][0]
            y[ind][1] = position_array[ind][1] * position_array[ind][1]
            y[ind][2] = position_array[ind][2] * position_array[ind][2]
        return y

    def sum_over_xyz(self, array):

        final = np.array([0, 0, 0])

        for list in array:
            final[0] += list[0]
            final[1] += list[1]
            final[2] += list[2]
        return final


    def radius_gyration(self):

        cen = self.geometric_center()
        total = 0
        for ind, pos in enumerate(self.position):
            total += self.mass[ind] * la.norm(np.subtract(self.mass[ind], cen)) ** 2

        return total / np.sum(self.mass)

    def get_complex(self):

        return ProteinCollection([self])

    def add_complex(self, com):

        self.complexes.append(com)

    def remove_complex(self, com):

        if com in self.complexes:
            self.complexes.remove(com)

    def update_complexes(self):

        for com in self.complexes:
            com.update()
