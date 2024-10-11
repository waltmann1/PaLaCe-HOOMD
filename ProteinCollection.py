from __future__  import division
import numpy as np
from numpy import linalg as la
import copy as cp
from Quaternion import QuaternionBetween
import hoomd
import random

class ProteinCollection(object):
    """
   seperating the protein list operations from pdb stuff
    """

    def __init__(self, chains):
        """

        :param PDBName: name of the pdb file representing the complex
        """

        self.chains = chains
        self.terminals = []
        self.rigid_count = int(np.sum([chain.rigid_centers for chain in self.chains]))
        if len(chains) > 0:
            self.name = chains[0].name
        else:
            self.name = "from_empty"
        self.dump_context = None
        self.reinit_context = None
        self.num_particles = np.sum([chain.num_particles for chain in self.chains])
        self.ss_bonds = []
        self.ss_bond_names = []
        self.ss_angles = []
        self.ss_angle_names = []
        self.ss_dihedrals = []
        self.ss_dihedral_names = []
        self.gsd_res_map = None

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
        f_out.write(str(self.num_particles))
        f_out.write('\n\n')
        g = 0
        for chain in self.chains:
            for i in range(chain.num_particles):
                # print(chain.type[i], chain.position[i])
                f_out.write('%s %f %f %f\n' % (chain.type[i], chain.position[i][0], chain.position[i][1],
                                               chain.position[i][2]))
        f_out.close()

    def create_res_map(self):

        if self.gsd_res_map is not None:
            return self.gsd_res_map
        self.gsd_res_map = []
        tag = self.rigid_count
        for rmi, chain in enumerate(self.chains):
            if chain.type[0] == 'tN':
                self.gsd_res_map.append([])
                for x in chain.res_indexes:
                    self.gsd_res_map[rmi].append(np.add(x, tag))
                tag += chain.num_particles
            return self.gsd_res_map

    def dump_backbone_pdb(self, dump_file=None):

        filename = dump_file
        if dump_file is None:
            filename = self.name + '_bb.pdb'
        elif dump_file[-4:] != '_bb.pdb':
            filename = dump_file + '_bb.pdb'
        f = open(filename, 'w')
        atom_number = 0
        res_number = 0
        for ind, chain in enumerate(self.chains):
            for ind2, res in enumerate(chain.res_indexes):
                new_res = cp.deepcopy(res)
                new_res = new_res[:5]
                new_res.remove(res[1])
                for atom in new_res:
                    string = "ATOM"
                    for i in range(6 - len(str(atom_number))):
                        string += " "
                    string += str(atom_number)
                    atom_number += 1
                    name = chain.type[atom]
                    if len(name) > 1:
                        name = "CA"
                    else:
                        name += " "
                    for i in range(4 - len(str(name))):
                        string += " "
                    string += name
                    r_name = chain.res_names[ind2]
                    for i in range(6 - len(str(r_name))):
                        string += " "
                    string += r_name
                    string += " "
                    string += chr(ord("A") + ind)
                    for i in range(4 - len(str(ind2))):
                        string += " "
                    string += str(ind2)
                    string += "      "
                    pos = str(str(chain.position[atom][0])[:6])
                    for i in range(6 - len(pos)):
                        pos += " "
                    string += pos + "  "
                    pos = str(str(chain.position[atom][1])[:6])
                    for i in range(6 - len(pos)):
                        pos += " "
                    string += pos + "  "
                    pos = str(str(chain.position[atom][2])[:6])
                    for i in range(6 - len(pos)):
                        pos += " "
                    string += pos + "  "

                    string += "1.00 10.00"
                    for i in range(11):
                        string += " "
                    string += name[:1]
                    string += "  \n"
                    f.write(string)
        f.close()

    def dump_res_map(self, dump_file=None):

        if self.gsd_res_map is None:
            self.create_res_map()
        filename = dump_file
        if dump_file is None:
            filename = self.name + '.resmap'
        elif dump_file[-4:] != '.resmap':
            filename = dump_file + '.resmap'
        f_out = open(filename, 'w')
        for ind, chain in enumerate(self.gsd_res_map):
            f_out.write('chain ' + str(ind))
            f_out.write('\n\n')
            for ind2, res in enumerate(chain):
                string = ""
                for ind3 in res:
                    string += str(ind3)
                    string += " "
                f_out.write(self.chains[ind].res_names[ind2] + " " + string + '\n')
            f_out.write('\n')
        f_out.close()

    def create_system(self):
        """

        :return: system object
        """
        if self.dump_context is None:
            self.dump_context = hoomd.context.initialize("")
        a_types = self.parse_angle_types()
        b_types = self.parse_bond_types() + ['N-H', 'C-O', 'O-H']
        d_types = self.parse_dihedral_types()
        p_types = self.parse_particle_names() + ["center" + str(i) for i in range(self.rigid_count)]
        dc_types = self.parse_dihedral_coupling_names()
        cen = self.geometric_center()
        self.shift([c * -1 for c in cen])
        print(self.rigid_count)
        snap = hoomd.data.make_snapshot(self.num_particles + self.rigid_count, particle_types=p_types,
                                        dihedral_types=d_types, bond_types=b_types, angle_types=a_types,
                                        dihedral_coupling_types=dc_types,
                                        box=hoomd.data.boxdim(L=int(self.max_radius() * 2 + 25)))
        snap.bonds.resize(0)

        for x in range(self.rigid_count):
            snap.particles.position[x] = self.center_of_mass(x)
            snap.particles.mass[x] = np.sum([chain.mass[i] for chain in self.chains for i in range(chain.num_particles)
                                             if chain.body[i] == x])
            snap.particles.typeid[x] = p_types.index("center" + str(x))
            snap.particles.body[x] = x
            snap.particles.moment_inertia[x] = self.moment_inertia(x)
        tag = self.rigid_count

        for chain in self.chains:

            for x in range(len(chain.get_bond_types())):
                bodies = [chain.body[i] for i in chain.get_bonds()[x] if chain.body[i] != -1]
                if len(bodies) == len(set(bodies)):
                    bond_number = snap.bonds.N + 1
                    snap.bonds.resize(bond_number)
                    snap.bonds.group[bond_number - 1] = np.add(chain.get_bonds()[x], tag)
                    snap.bonds.typeid[bond_number - 1] = b_types.index(chain.get_bond_types()[x])
            for x in range(len(chain.get_angle_types())):
                bodies = [chain.body[i] for i in chain.get_angles()[x] if chain.body[i] != -1]
                if len(bodies) == len(set(bodies)):
                    angle_number = snap.angles.N + 1
                    snap.angles.resize(angle_number)
                    snap.angles.group[angle_number - 1] = np.add(chain.get_angles()[x], tag)
                    snap.angles.typeid[angle_number - 1] = a_types.index(chain.get_angle_types()[x])
            for x in range(len(chain.get_dihedral_types())):
                bodies = [chain.body[i] for i in chain.get_dihedrals()[x] if chain.body[i] != -1]
                if len(bodies) == len(set(bodies)):
                    dihedral_number = snap.dihedrals.N + 1
                    snap.dihedrals.resize(dihedral_number)
                    snap.dihedrals.group[dihedral_number - 1] = np.add(chain.get_dihedrals()[x], tag)
                    snap.dihedrals.typeid[dihedral_number - 1] = d_types.index(chain.get_dihedral_types()[x])
            for x in range(len(chain.get_dihedral_coupling_types())):
                bodies = [chain.body[i] for i in set(chain.get_dihedrals_coupling()[x]) if chain.body[i] != -1]
                if len(bodies) == len(set(bodies)):
                    dihedral_coupling_number = snap.dihedrals_coupling.N + 1
                    snap.dihedrals_coupling.resize(dihedral_coupling_number)
                    snap.dihedrals_coupling.group[dihedral_coupling_number - 1] = \
                    np.add(chain.get_dihedrals_coupling()[x], tag)
                    snap.dihedrals_coupling.typeid[dihedral_coupling_number - 1] = \
                    dc_types.index(chain.get_dihedral_coupling_types()[x])
            for x in range(chain.num_particles):
                snap.particles.position[x + tag] = chain.position[x]
                snap.particles.mass[x + tag] = chain.mass[x]
                snap.particles.typeid[x + tag] = p_types.index(chain.type[x])
                snap.particles.body[x + tag] = chain.body[x]
                snap.particles.charge[x + tag] = chain.charge[x]
            tag += chain.num_particles

        for x in range(len(self.ss_bonds)):
                bond_number = snap.bonds.N + 1
                snap.bonds.resize(bond_number)
                snap.bonds.group[bond_number - 1] = np.add(self.ss_bonds[x], self.rigid_count)
                snap.bonds.typeid[bond_number - 1] = b_types.index(self.ss_bond_names[x])

        for x in range(len(self.ss_angles)):
                angle_number = snap.angles.N + 1
                snap.angles.resize(angle_number)
                snap.angles.group[angle_number - 1] = np.add(self.ss_angles[x], self.rigid_count)
                snap.angles.typeid[angle_number - 1] = a_types.index(self.ss_angle_names[x])
        for x in range(len(self.ss_dihedral_names)):
                dihedral_number = snap.dihedrals.N + 1
                snap.dihedrals.resize(dihedral_number)
                snap.dihedrals.group[dihedral_number - 1] = np.add(self.ss_dihedrals[x], self.rigid_count)
                snap.dihedrals.typeid[dihedral_number - 1] = d_types.index(self.ss_dihedral_names[x])

        sys = hoomd.init.read_snapshot(snap)
        self.shift([c for c in cen])

        return sys

    def create_system_shell(self):
        """

        :return: system object
        """
        print('create shell')
        if self.dump_context is None:
            self.dump_context = hoomd.context.initialize("")
        a_types = self.parse_angle_types()
        b_types = self.parse_bond_types() + ['N-H', 'C-O', 'O-H']
        d_types = self.parse_dihedral_types()
        p_types = self.parse_particle_names() + ["center" + str(i) for i in range(self.rigid_count)]
        dc_types = self.parse_dihedral_coupling_names()
        cen = self.geometric_center()
        self.shift([c * -1 for c in cen])
        num = self.rigid_count + np.sum([len(chain.active_list) for chain in self.chains])
        snap = hoomd.data.make_snapshot(num, particle_types=p_types,
                                        dihedral_types=d_types, bond_types=b_types, angle_types=a_types,
                                        dihedral_coupling_types=dc_types,
                                        box=hoomd.data.boxdim(L=int(self.max_radius() * 2 + 25)))
        snap.bonds.resize(0)

        for x in range(self.rigid_count):
            snap.particles.position[x] = self.center_of_mass(x)
            snap.particles.mass[x] = np.sum([chain.mass[i] for chain in self.chains for i in range(chain.num_particles)
                                             if chain.body[i] == x])
            snap.particles.typeid[x] = p_types.index("center" + str(x))
            snap.particles.body[x] = x
            snap.particles.moment_inertia[x] = self.moment_inertia(x)
        tag = self.rigid_count
        for chain in self.chains:
            current = 0
            tags_array = self.active_tags_array(tag, chain.active_list, chain.num_particles)
            for x in range(chain.num_particles):
                if x in chain.active_list:
                    print(current + tag)
                    snap.particles.position[current + tag] = chain.position[x]
                    snap.particles.mass[current + tag] = chain.mass[x]
                    snap.particles.typeid[current + tag] = p_types.index(chain.type[x])
                    snap.particles.body[current + tag] = chain.body[x]
                    snap.particles.charge[current + tag] = chain.charge[x]
                    current += 1
            for x in range(len(chain.get_bond_types())):
                bodies = [chain.body[i] for i in chain.get_bonds()[x] if i in chain.active_list]
                if len(bodies) == 2:
                    bond_number = snap.bonds.N + 1
                    snap.bonds.resize(bond_number)
                    snap.bonds.group[bond_number - 1] = [tags_array[i] for i in chain.get_bonds()[x]]
                    snap.bonds.typeid[bond_number - 1] = b_types.index(chain.get_bond_types()[x])
            for x in range(len(chain.get_angle_types())):
                bodies = [chain.body[i] for i in chain.get_angles()[x] if i in chain.active_list]
                if len(bodies) == 3:
                    angle_number = snap.angles.N + 1
                    snap.angles.resize(angle_number)
                    snap.angles.group[angle_number - 1] = [tags_array[i] for i in chain.get_angles()[x]]
                    snap.angles.typeid[angle_number - 1] = a_types.index(chain.get_angle_types()[x])
            for x in range(len(chain.get_dihedral_types())):
                bodies = [chain.body[i] for i in chain.get_dihedrals()[x] if i in chain.active_list]
                if len(bodies) == 4:
                    dihedral_number = snap.dihedrals.N + 1
                    snap.dihedrals.resize(dihedral_number)
                    snap.dihedrals.group[dihedral_number - 1] = [tags_array[i] for i in chain.get_dihedrals()[x]]
                    snap.dihedrals.typeid[dihedral_number - 1] = d_types.index(chain.get_dihedral_types()[x])
            for x in range(len(chain.get_dihedral_coupling_types())):
                bodies = [chain.body[i] for i in set(chain.get_dihedrals_coupling()[x]) if i in chain.active_list]
                if len(bodies) == 5:
                    dihedral_coupling_number = snap.dihedrals_coupling.N + 1
                    snap.dihedrals_coupling.resize(dihedral_coupling_number)
                    snap.dihedrals_coupling.group[dihedral_coupling_number - 1] = \
                        [tags_array[i] for i in chain.get_dihedrals_coupling()[x]]
                    snap.dihedrals_coupling.typeid[dihedral_coupling_number - 1] = \
                        dc_types.index(chain.get_dihedral_coupling_types()[x])
            tag += len(chain.active_list)

        """
        for x in range(len(self.ss_bonds)):
            bond_number = snap.bonds.N + 1
            snap.bonds.resize(bond_number)
            snap.bonds.group[bond_number - 1] = np.add(self.ss_bonds[x], self.rigid_count)
            snap.bonds.typeid[bond_number - 1] = b_types.index(self.ss_bond_names[x])

        for x in range(len(self.ss_angles)):
            angle_number = snap.angles.N + 1
            snap.angles.resize(angle_number)
            snap.angles.group[angle_number - 1] = np.add(self.ss_angles[x], self.rigid_count)
            snap.angles.typeid[angle_number - 1] = a_types.index(self.ss_angle_names[x])
        for x in range(len(self.ss_dihedral_names)):
            dihedral_number = snap.dihedrals.N + 1
            snap.dihedrals.resize(dihedral_number)
            snap.dihedrals.group[dihedral_number - 1] = np.add(self.ss_dihedrals[x], self.rigid_count)
            snap.dihedrals.typeid[dihedral_number - 1] = d_types.index(self.ss_dihedral_names[x])
        """

        sys = hoomd.init.read_snapshot(snap)
        self.shift([c for c in cen])

        return sys

    def get_active_tag(self, tag, cal, ind):

        return tag + len([active for active in cal if active < ind])

    def active_tags(self, tag, cal, list):

        return [self.get_active_tag(tag, cal, ind) for ind in list]

    def active_tags_array(self, tag, cal, num_particles):
        tags = []
        current = tag
        cal_index =0
        for part in range(num_particles):
            if part > cal[cal_index]:
                cal_index += 1
                current += 1
            tags.append(current)
        return tags

    def parse_angle_types(self):

        a_types = []
        for chain in self.chains:
            for c in chain.angle_names:
                if c not in a_types:
                    a_types.append(c)
        a_types += list(set(self.ss_angle_names))
        return a_types

    def parse_bond_types(self):

        a_types = []
        for chain in self.chains:
            for c in chain.bond_names:
                if c not in a_types:
                    a_types.append(c)
        a_types += list(set(self.ss_bond_names))
        return a_types

    def parse_dihedral_types(self):

        a_types = []
        for chain in self.chains:
            for c in chain.dihedral_names:
                if c not in a_types:
                    a_types.append(c)
        a_types += list(set(self.ss_dihedral_names))
        return a_types

    def parse_dihedral_coupling_names(self):
        a_types = []
        for chain in self.chains:
            for c in chain.dihedral_coupling_names:
                if c not in a_types:
                    a_types.append(c)
        return a_types

    def parse_particle_names(self):
        a_types = []
        for chain in self.chains:
            for c in chain.type:
                if c not in a_types:
                    a_types.append(c)
        return a_types

    def dump_gsd(self, dump_file=None, shell=False):
        """

        :param dump_file: name of the file to dump the xyz to
        :return: nothing
        """

        filename = dump_file
        if dump_file is None:
            filename = self.name + '.gsd'
        elif dump_file[-4:] != '.gsd':
            filename = dump_file + '.gsd'
        #self.center_at_origin()
        if not shell:
            sys = self.create_system()
        else:
            self.create_active_lists()
            sys = self.create_system_shell()

        hoomd.dump.gsd(filename=filename, period=None, group=hoomd.group.all(), overwrite=True)
        return filename

    def is_number(self, s):
        try:
            float(s)
            return True
        except ValueError:
            return False

    def my_int(self, char):

        if self.is_number(char):
            return int(char)
        else:
            if not self.is_number(char[0]):
                return int(ord(char) - 64)
            else:
                ind = 0
                c = char[ind]
                while self.is_number(c):
                    ind += 1
                    c = char[ind]
                return int(char[:ind])

    def update(self):
        self.num_particles = np.sum([chain.num_particles for chain in self.chains])

    def add_chain(self, chain):

        self.chains.append(chain)
        chain.add_complex(self)
        self.update()

    def remove_chain(self, chain):

        self.chains.remove(chain)
        chain.remove_complex(self)
        self.update()

    def remove_chain_by_index(self, index):
        chain = self.chains[index]
        # TODO handle terminals and ss and stuff

    def geometric_center(self):

        pos = [np.mean(c.position, axis=0) for c in self.chains]
        weights = [len(c.position) for c in self.chains]
        total = np.sum(weights)
        weights = [float(weight)/float(total) for weight in weights]
        center = np.array([0, 0, 0])
        for ind, p in enumerate(pos):
            #print(weights[ind])
            center = np.add(center, np.multiply(p, weights[ind]))
        return center

    def shift(self, vec):
        """

        :param vec: vector to shift by
        :return:
        """
        for c in self.chains:
            c.shift(vec)
    def max_radius(self):

        t = [-1 * c for c in self.geometric_center()]
        #print(t)
        self.shift(t)
        rs = [la.norm(pos) for chain in self.chains for pos in chain.position]
        self.shift([-1 * g for g in t])
        #print(np.max(rs))
        #quit()
        return np.max(rs)

    def make_rigid(self, chain, res_index, body_tag=None):

        if body_tag > self.rigid_count:
            raise ValueError("body tags must be consecutive", self.rigid_count, body_tag)
        if body_tag is None or body_tag == self.rigid_count:
            body_tag = self.rigid_count
            self.rigid_count += 1

        self.chains[chain].make_rigid(res_index, body_tag)

    def moment_inertia(self, center_index):
        mass_array = np.array([chain.mass[i] for chain in self.chains for i in range(chain.num_particles)
                               if chain.body[i] == center_index])
        position_array = np.array([chain.position[i] for chain in self.chains for i in range(chain.num_particles)
                                   if chain.body[i] == center_index])

        cen = self.center_of_mass_arrays(position_array, mass_array)
        position_array = np.subtract(position_array, cen)
        return self.sum_over_xyz(self.pos_x_mass(self.pos_squared(position_array), mass_array))

    def center_of_mass(self, center_index):
        mass_array = np.array([chain.mass[i] for chain in self.chains for i in range(chain.num_particles)
                               if chain.body[i] == center_index])
        position_array = np.array([chain.position[i] for chain in self.chains for i in range(chain.num_particles)
                                   if chain.body[i] == center_index])
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

    def make_connections_rigid(self):
        count = self.rigid_count
        for ind, chain in enumerate(self.chains):
            con = chain.connections
            for c in con:
                self.make_rigid(ind, c, body_tag=count)
                self.make_rigid(ind, c - 1, body_tag=count)
                count += 1

    def make_rigid_in_sphere(self, radius):
        tag = self.rigid_count
        center = self.geometric_center()

        for ind, chain in enumerate(self.chains):
            for i in range(chain.residues):
                iac = chain.resnum_to_iac(i)
                if la.norm(np.subtract(chain.position[iac], center)) < radius:
                    self.make_rigid(ind, i, tag)

    def make_rigid_in_cylinder_xy(self, radius, z=100000):
        tag = self.rigid_count
        center = self.geometric_center()

        for ind, chain in enumerate(self.chains):
            for i in range(chain.residues):
                iac = chain.resnum_to_iac(i)
                vec = [np.subtract(chain.position[iac], center)[0],np.subtract(chain.position[iac], center)[1],0]
                if la.norm(vec) < radius and np.abs(np.subtract(chain.position[iac], center)[2]) < z:
                    self.make_rigid(ind, i, tag)

    def make_ends_not_rigid(self,nterm, cterm):
        tag = self.rigid_count
        center = self.geometric_center()

        for ind, chain in enumerate(self.chains):
            for i in range(chain.residues):
                iac = chain.resnum_to_iac(i)
                vec = [np.subtract(chain.position[iac], center)[0],np.subtract(chain.position[iac], center)[1],0]
                if i > nterm and i < chain.residues - cterm:
                    self.make_rigid(ind, i, tag)


    def add_protein_complex(self, complex):

        new_complex = cp.deepcopy(complex)
        for chain in new_complex.chains:
            for ind, body in enumerate(chain.body):
                if body != - 1:
                    #print(chain.body[ind])
                    chain.body[ind] += self.rigid_count
                    #print(chain.body[ind])
        self.rigid_count += new_complex.rigid_count
        for ind, bond in enumerate(new_complex.ss_bonds):
            new_complex.ss_bonds[ind] = list(np.add(self.num_particles, bond))
        for ind, angle in enumerate(new_complex.ss_angles):
            new_complex.ss_angles[ind] = list(np.add(self.num_particles, angle))
        for ind, dihedral in enumerate(new_complex.ss_dihedrals):
            new_complex.ss_dihedrals[ind] = list(np.add(self.num_particles, dihedral))
        for chain in new_complex.chains:
            self.add_chain(chain)

    def orient(self, vec1, vec2):
        temp = self.geometric_center()
        self.shift(np.multiply(-1, temp))
        q = QuaternionBetween(vec1, vec2)
        for chain in self.chains:
            for x in range(chain.num_particles):
                chain.position[x] = q.orient(chain.position[x])
        self.shift(temp)

    def orient_quaternion(self, q):
        temp = self.geometric_center()
        self.shift(np.multiply(-1, temp))
        for chain in self.chains:
            for x in range(chain.num_particles):
                chain.position[x] = q.orient(chain.position[x])
        self.shift(temp)

    def populate_box(self, positions, rotations):

        if not (len(positions) == len(rotations)):
            raise ValueError("Lengths of posiitons and rotations must be equal.", len(positions, len(rotations)))
        num = len(rotations)

        news = []

        self.center_at_origin()
        for ind in range(1, num):
            news.append(cp.deepcopy(self))

        self.orient_quaternion(rotations[0])
        self.shift(positions[0])

        for ind, new in enumerate(news):
            new.orient_quaternion(rotations[ind + 1])
            new.shift(positions[ind + 1])
            print(new.geometric_center())
            self.add_protein_complex(new)

    def center_at_origin(self):
        temp = self.geometric_center()
        self.shift(np.multiply(-1, temp))

    def pair_config(self, random_orientation=False, distance=None):

        self.center_at_origin()
        if distance is None:
            positions = [[-1 * self.max_radius(), 0, 0], [1 * self.max_radius(), 0, 0]]
        else:
            positions = [[-1 * distance/2, 0, 0], [1 * distance/2, 0, 0]]
        rotations = self.get_zero_rotation_array(2)
        if random_orientation:
            rotations = self.get_random_orientation_array(2)
        self.populate_box(positions, rotations)

    def triangle_config(self, random_orientation=False):
        self.center_at_origin()
        positions = [[-1 * self.max_radius(), 0, 0], [1 * self.max_radius(), 0, 0], ]
        rotations = self.get_zero_rotation_array(2)
        if random_orientation:
            rotations = self.get_random_orientation_array(2)
        self.populate_box(positions, rotations)


    def shell_configuration(self, r, n=None):

        self.center_at_origin()
        if n is None:
            single_area = self.max_radius()**2
            sa = 4 * r**2
            n = int(sa / single_area)

        points = self.unit_sphere(n)
        positions = []
        rotations = []

        for point in points:
            print(point)
            positions.append(np.multiply(point, r))
            rotations.append(QuaternionBetween([0, 0, 1], point))

        print("populate")
        print(positions)

        self.populate_box(positions, rotations)

    def unit_sphere(self, n, randomize=False):
        #rnd = 1.
        #if randomize:
        #    rnd = random.random() * n

        points = []
        offset = 2. / n
        increment = np.pi * (3. - np.sqrt(5.));

        for i in range(n):
            y = ((i * offset) - 1) + (offset / 2);
            r = np.sqrt(1 - pow(y, 2))

            #phi = ((i + rnd) % n) * increment
            phi = i * increment

            x = np.cos(phi) * r
            z = np.sin(phi) * r

            points.append([x, y, z])

        return points

    def get_zero_rotation_array(self, len):

        return [QuaternionBetween([1,0,0], [1,0,0]) for _ in range(len)]

    def get_random_orientation_array(self, len):

        array = []
        for i in range(len):
            y = np.random.randint(-100, 100, size=6)
            one = np.divide(y[:3], la.norm(y[:3]))
            two = np.divide(y[3:], la.norm(y[3:]))
            array.append(QuaternionBetween(one, two))
        return array

    def set_gsd_positions(self, gsd_name, frame):

        if self.reinit_context is None:
            self.reinit_context = hoomd.context.initialize("")
        sys = hoomd.init.read_gsd(gsd_name, frame=frame)
        skip = self.rigid_count
        chain = 0
        pos_counter = 0
        for ind in range(skip, skip + self.num_particles):
            part = sys.particles.get(ind)
            if self.chains[chain].num_particles == pos_counter:
                chain += 1
                pos_counter = 0
            self.chains[chain].position[pos_counter] = list(part.position)

    def create_active_lists(self):

        for chain in self.chains:
            for ind in range(chain.num_particles):
                if chain.body[ind] == -1:
                    chain.active_list.append(ind)
                else:
                    for bond in chain.bonds:
                        if ind in bond and (chain.body[bond[0]] == -1 or chain.body[bond[1]] == -1):
                            chain.active_list.append(ind)
            chain.active_list = list(set(chain.active_list))
            chain.active_list.sort()

    def radius_gyration(self):

        cen = self.geometric_center()
        total = 0
        total_mass = 0
        for chain in self.chain:
            for ind, pos in enumerate(self.position):
                total_mass += chain.mass[ind]
                total += chain.mass[ind] * la.norm(np.subtract(chain.mass[ind], cen)) ** 2

        return total / total_mass
