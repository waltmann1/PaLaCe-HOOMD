from __future__ import division
import numpy as np
from ProteinAbs import ProteinAbs
from Quaternion import QuaternionBetween
from ProteinCollection import ProteinCollection


class ChainTooShortError(IndexError):
    def __init__(self,*args,**kwargs):
        IndexError.__init__(self,*args,**kwargs)


class ProteinBackbone(ProteinAbs):

    def __init__(self, residues, chain_index=0, list_470=[]):

        if residues < 1:
            raise ChainTooShortError('must be at least one residue')

        super(ProteinBackbone, self).__init__(residues, chain_index=chain_index)
        current = [0, 0, 0]

        cn_vector = [1.3 * np.cos(np.deg2rad(33)), -1.3 * np.sin(np.deg2rad(33)), 0]
        nh_vector = [0, -1, 0]
        co_vector = [0, -1.24, 0]
        nc_vector = [1.5 * np.cos(np.deg2rad(24)), 1.5 * np.sin(np.deg2rad(24)), 0]
        cc_vector = [1.5 * np.cos(np.deg2rad(31)), -1.5 * np.sin(np.deg2rad(31)), 0]
        count = 0
        self.res_names = ['NOR' for _ in range(self.residues)]
        self.res_indexes = [[] for _ in range(self.residues)]
        self.list_470 = list_470
        for i in range(residues):

            current = np.add(current, cn_vector)

            self.mass.append(14)
            self.position.append(current)
            self.type.append('N')
            self.body.append(-1)
            self.res_indexes[i].append(count)
            count += 1

            self.mass.append(1)
            self.position.append(np.add(current, nh_vector))
            self.type.append('H')
            self.body.append(-1)
            self.res_indexes[i].append(count)
            count += 1

            current = np.add(current, nc_vector)
            self.mass.append(13)
            self.position.append(current)
            self.type.append('alphaC')
            self.body.append(-1)
            self.res_indexes[i].append(count)
            count += 1

            current = np.add(current, cc_vector)
            self.mass.append(12)
            self.position.append(current)
            self.type.append('C')
            self.body.append(-1)
            self.res_indexes[i].append(count)
            count += 1

            self.mass.append(16)
            self.position.append(np.add(current, co_vector))
            self.type.append('O')
            self.body.append(-1)
            self.res_indexes[i].append(count)
            count += 1

            nc_vector[1] = -nc_vector[1]
            nh_vector[1] = -nh_vector[1]
            co_vector[1] = -co_vector[1]
            cn_vector[1] = -cn_vector[1]
            cc_vector[1] = -cc_vector[1]
        self.num_particles = len(self.mass)
        self.name = 'backbone_' + str(self.residues)
        self.add_angles()
        self.add_bonds()
        self.add_dihedrals()
        self.add_dihedrals_coupling()
        y = np.mean(self.position, axis=0)
        self.position = np.subtract(self.position, y)
        self.charge = [0] * self.num_particles
        self.connections = []
        self.add_hb_exclusion_bonds()
        self.type[0] = 'tN'
        self.type[len(self.type) - 2] = 'tC'

        #quit()

    def add_bonds(self):

        for i in range(self.residues):
            self.bonds.append([5 * i, 5*i + 2])
            self.bond_names.append('N-alphaC')
            self.bonds.append([5 * i + 2, 5 * i + 3])
            self.bond_names.append('alphaC-C')
            if i < self.residues - 1:
                self.bonds.append([5 * i + 3, 5 * i + 5])
                self.bond_names.append('C-N')

    def add_hb_exclusion_bonds(self):

        for res_ind, res in enumerate(self.res_indexes[:-1]):
            self.exclude_hb_residues_to_the_right(res_ind)
            self.exclude_hb_resiudes_own_residue(res_ind)
        self.exclude_hb_resiudes_own_residue(len(self.res_indexes) - 1)


    def exclude_hb_resiudes_own_residue(self, res_ind):

        types = ['C', 'N', 'O', 'H', 'tC', 'tN']

        for i, ind in enumerate(self.res_indexes[res_ind]):
            for ind2 in self.res_indexes[res_ind][i + 1:]:
                if self.type[ind] in types and self.type[ind2] in types:
                    self.bonds.append([ind, ind2])
                    #print(self.type[ind], self.type[ind2])
                    self.bond_names.append('exclusion')
                    #quit()




    def exclude_hb_residues_to_the_right(self, res_ind):

        types = ['C', 'N', 'O', 'H', 'tC', 'tN']
        for ind in self.res_indexes[res_ind]:
            for ind2 in self.res_indexes[res_ind + 1]:
                if self.type[ind] in types and self.type[ind2] in types:
                    self.bonds.append([ind, ind2])
                    print([ind,ind2])
                    self.bond_names.append('exclusion')



    def add_angles(self):

        for i in range(self.residues):
            self.angles.append([5 * i, 5*i + 2, 5*i + 3])
            self.angle_names.append('N-alphaC-C')
            if i < self.residues - 1:
                self.angles.append([5 * i + 2, 5 * i + 3, 5*i + 5])
                self.angle_names.append('alphaC-C-N')
                self.angles.append([5 * i + 3, 5 * i + 5, 5 * i + 7])
                self.angle_names.append('C-N-alphaC')

    def add_dihedrals(self):

        self.dihedrals = [[5 * i + 2, 5 * i + 3, 5 * i + 5, 5 * i + 7] for i in range(self.residues - 1)]
        self.dihedral_names = ['alphaC-C-N-alphaC'] * (self.residues - 1)

    def add_dihedrals_coupling(self):

        self.dihedrals_coupling = [[5 * i + 3, 5 * i + 5, 5 * i + 7, 5 * i + 8, 5 * i + 5, 5 * i + 7, 5 * i + 8,
                                    5 * i + 10] for i in range(self.residues - 2)]
        self.dihedral_coupling_names = ['NOR'] * (self.residues - 2)

    def add_residue(self, name, index_alpha_carbon, data=None):

        # if name == 'HOH':
            # return

        if index_alpha_carbon in self.used_alpha_carbons:
                raise ValueError("This alpha carbon has already been assigned an amino acid")

        if self.type[index_alpha_carbon] != 'alphaC':
            raise ValueError("index " + str(index_alpha_carbon) + " is not the index of an alpha carbon. "
                             + str(self.type[index_alpha_carbon]))

        self.used_alpha_carbons.append(index_alpha_carbon)
        self.res_names[self.iac_to_resnum(index_alpha_carbon)] = name

        imp = __import__(name, fromlist=[''])
        error_470 = self.iac_to_resnum(index_alpha_carbon) in self.list_470
        res = getattr(imp, name)(data=data, error_470=error_470)
        if not error_470:
            res.add_to_backbone(index_alpha_carbon, self)
        else:
            res.add_backbone_positions(index_alpha_carbon, self)
        # print(len(self.mass), len(self.position), self.num_particles, res.code)

    def dump_xyz2(self):
        f = open('backbone_' + str(self.residues) + '.xyz', 'w')
        beads = len(self.mass)
        f.write(str(beads) + '\n' + '\n')
        for i in range(beads):
            f.write(self.type[i] + ' ' + str(self.position[i][0]) + ' ' + str(self.position[i][1]) + ' ' +
                    str(self.position[i][2]) + '\n')
        f.close()

    def iac_to_resnum(self, iac):

        if self.type[iac][:6] != 'alphaC':
            raise ValueError(str(iac) + "is not the index of an alpha carbon")

        res_count = 0
        for i in range(0, iac):
            if self.type[i][:6] == 'alphaC':
                res_count += 1

        return res_count

    def resnum_to_iac(self, resnum):

        if resnum >= self.residues:
            raise ValueError("There are less than " + str(resnum+1) + " residues.")
        res_count = 0
        i = 0
        while res_count < resnum + 1:
            i += 1
            if self.type[i][:6] == 'alphaC':
                res_count += 1
        return i

    def populate(self, resname):
        """
        makes all residues that given residue
        :param resname: name of resiude
        :return:
        """

        for i in range(self.residues):
            self.add_residue(resname, self.resnum_to_iac(i))

    def use_chain(self, chain, denature=False):
        """

        :param chain: biopython chain object
        :return:
        """

        gen = chain.get_residues()
        residues = [r for r in gen]
        for i in range(len(self.residues)):
            data = [a for a in residues[i]]
            if denature:
                data = None
            if i not in self.used_alpha_carbons:
                self.add_residue(residues[i].resname, self.resnum_to_iac(i), data=data)

    def add_chain_to_C_terminal(self, bb2):
        """

        :return: void chain changed internally
        """
        self.connections.append(self.residues)
        terminal_index = 0
        hit = False
        for ind in range(self.num_particles - 1, 0, -1):
            if self.type[ind] == 'tC' and not hit:
                terminal_index = ind
                hit = True

        if self.type[terminal_index] != 'tC':
            raise ValueError('Terminal index must be a terminal carbon, tC not a ' + str(self.type[terminal_index]), terminal_index)


        self.type[terminal_index] = 'C'
        bb2.type[0] = 'N'

        # calculate the bond vector between the two chains
        terminal_position = self.position[terminal_index]
        final_bond = np.subtract(self.position[terminal_index - 2], terminal_position)
        new_bond = np.multiply(final_bond, -1.3/1.5)

        # reset the second chain pos to start at 0 and then apply the new bond vector
        bb2.position = np.subtract(bb2.position, bb2.position[1])
        new_position = np.add(terminal_position, new_bond)
        bb2.position = np.add(bb2.position, new_position)

        self.residues += bb2.residues
        self.res_names += bb2.res_names
        temp = []
        for ind, res in enumerate(bb2.res_indexes):
            temp.append(list(np.add(res, self.num_particles)))
        self.res_indexes += temp
        self.used_alpha_carbons += list(np.add(bb2.used_alpha_carbons, self.num_particles))
        self.dihedrals_coupling += list(np.add(bb2.dihedrals_coupling, self.num_particles))
        self.dihedral_coupling_names += bb2.dihedral_coupling_names
        self.dihedrals += list(np.add(bb2.dihedrals, self.num_particles))
        self.dihedral_names += bb2.dihedral_names
        self.angles += list(np.add(bb2.angles, self.num_particles))
        self.angle_names += bb2.angle_names
        self.bonds += list(np.add(bb2.bonds, self.num_particles))
        self.bond_names += bb2.bond_names
        self.type += bb2.type
        self.position = list(self.position) + list(bb2.position)
        self.mass += list(bb2.mass)
        self.charge += bb2.charge
        for bod in bb2.body:
            if bod != -1 and bod != 666:
                bod += self.num_particles
        self.body += bb2.body

        self.bonds.append([terminal_index, self.num_particles])
        self.bond_names.append('C-N')
        self.angles.append([terminal_index - 1, terminal_index, self.num_particles])
        self.angle_names.append('alphaC-C-N')
        self.angles.append([terminal_index, self.num_particles, self.num_particles + 2])
        self.angle_names.append('C-N-alphaC')
        self.dihedrals.append([terminal_index - 1, terminal_index, self.num_particles, self.num_particles + 2])
        if bb2.type[3] == 'alphaCp':
            self.dihedral_names.append('alphaC-C-N-alphaCproline')
        else:
            self.dihedral_names.append('alphaC-C-N-alphaC')
        ti = terminal_index
        imp = __import__(self.res_names[self.iac_to_resnum(ti - 1)], fromlist=[''])
        rdue = getattr(imp, self.res_names[self.iac_to_resnum(ti - 1)])()
        self.dihedrals_coupling.append([ti - 5, ti - 3, ti - 1, ti, ti - 3, ti - 1, ti, self.num_particles])
        self.dihedral_coupling_names.append(rdue.dc_name)

        imp = __import__(self.res_names[self.residues - bb2.residues], fromlist=[''])
        rdue = getattr(imp, self.res_names[self.residues - bb2.residues])()
        if rdue.dc_name == 'PRO':
            self.dihedral_coupling_names[-1] = 'PRP'
        if bb2.residues > 1:
            self.dihedral_coupling_names.append(rdue.dc_name)
            self.dihedrals_coupling.append([ti, self.num_particles, self.num_particles + 2, self.num_particles + 3,
                                        self.num_particles, self.num_particles + 2, self.num_particles + 3,
                                        self.num_particles + 5])
        self.num_particles += bb2.num_particles
        self.exclude_hb_residues_to_the_right(self.connections[-1] - 1)
        self.update_complexes()

    def add_chain_to_N_terminal(self, bb2):

        """

        :return: void chain changed internally
        """

        for connect in self.connections:
            connect += bb2.residues
        self.connections.append(bb2.residues)

        terminal_index = 0
        hit = False
        for ind in range(bb2.num_particles - 1, 0, -1):
            if bb2.type[ind] == 'tC' and not hit:
                terminal_index = ind
                hit = True

        if bb2.type[terminal_index] != 'tC':
            raise ValueError('Terminal index must be a terminal carbon, tC not a ' + str(bb2.type[terminal_index]),
                             terminal_index)



        self.type[0] = 'N'


        # calculate the bond vector between the two chains
        terminal_position = bb2.position[terminal_index]
        final_bond = np.subtract(self.position[2], self.position[0])
        new_bond = np.multiply(final_bond, -1)

        # reset this chain pos to start at 0 and then apply the new bond vector
        bb2.position = np.subtract(bb2.position, bb2.position[terminal_index])
        new_position = np.add(self.position[0], new_bond)
        bb2.position = np.add(bb2.position, new_position)
        # avoid 180 degree angle
        bb2.position[terminal_index] = np.add(bb2.position[terminal_index], [0,0,.2])

        self.residues += bb2.residues
        self.res_names = bb2.res_names + self.res_names
        for ind, res in enumerate(self.res_indexes):
            self.res_indexes[ind] = list(np.add(res, bb2.num_particles))
        self.res_indexes = bb2.res_indexes + self.res_indexes
        self.used_alpha_carbons = bb2.used_alpha_carbons + list(np.add(self.used_alpha_carbons, bb2.num_particles))
        self.dihedrals_coupling = bb2.dihedrals_coupling + list(np.add(self.dihedrals_coupling, bb2.num_particles))
        self.dihedral_coupling_names = bb2.dihedral_coupling_names + self.dihedral_coupling_names
        self.dihedrals = bb2.dihedrals + list(np.add(self.dihedrals, bb2.num_particles))
        self.dihedral_names = bb2.dihedral_names + self.dihedral_names
        self.angles = bb2.angles + list(np.add(self.angles, bb2.num_particles))
        self.angle_names = bb2.angle_names + self.angle_names
        # print('bonds before ', len(self.bonds))
        self.bonds = bb2.bonds + list(np.add(self.bonds, bb2.num_particles))
        # print('bonds after ', len(self.bonds))
        # print('names before ', len(self.bond_names))
        self.bond_names = bb2.bond_names + self.bond_names
        # print('names after ', len(self.bond_names))
        self.position = list(bb2.position) + list(self.position)
        self.mass = list(bb2.mass) + list(self.mass)
        self.charge = bb2.charge + self.charge
        for bod in self.body:
            if bod != -1 and bod != 666:
                bod += self.num_particles
        self.body = bb2.body + self.body

        self.bonds.append([terminal_index, bb2.num_particles])
        self.bond_names.append('C-N')
        self.angles.append([terminal_index - 1, terminal_index, bb2.num_particles])
        self.angle_names.append('alphaC-C-N')
        self.angles.append([terminal_index, bb2.num_particles, bb2.num_particles + 2])
        self.angle_names.append('C-N-alphaC')
        self.dihedrals.append([terminal_index - 1, terminal_index, bb2.num_particles, bb2.num_particles + 2])
        #print([terminal_index - 1, terminal_index, bb2.num_particles, bb2.num_particles + 2])
        if self.type[2] == 'alphaCp':
            self.dihedral_names.append('alphaC-C-N-alphaCproline')
        else:
            self.dihedral_names.append('alphaC-C-N-alphaC')

        ti = terminal_index
        imp = __import__(bb2.res_names[bb2.iac_to_resnum(ti - 1)], fromlist=[''])
        rdue = getattr(imp, bb2.res_names[bb2.iac_to_resnum(ti - 1)])()
        self.dihedrals_coupling.append([ti - 5, ti - 3, ti - 1, ti, ti - 3, ti - 1, ti, bb2.num_particles])
        self.dihedral_coupling_names.append(rdue.dc_name)

        imp = __import__(self.res_names[bb2.residues], fromlist=[''])
        rdue = getattr(imp, self.res_names[bb2.residues])()
        if rdue.dc_name == 'PRO':
            self.dihedral_coupling_names[-1] = 'PRP'
        if self.residues - bb2.residues > 1:
            self.dihedral_coupling_names.append(rdue.dc_name)
            self.dihedrals_coupling.append([ti, bb2.num_particles, bb2.num_particles + 2, bb2.num_particles + 3,
                                            bb2.num_particles, bb2.num_particles + 2, bb2.num_particles + 3,
                                            bb2.num_particles + 5])
        self.num_particles += bb2.num_particles
        self.type = bb2.type + self.type
        self.exclude_hb_residues_to_the_right(self.connections[-1] - 1)
        #quit()
        self.update_complexes()

        self.type[terminal_index] = 'C'


    def increment_num_particles(self, index_alpha_carbon):

        self.res_indexes[self.iac_to_resnum(index_alpha_carbon)].append(self.num_particles)
        self.num_particles += 1

    def curve(self, function, start=[0, 0, 0]):

        self.position = np.subtract(self.position, self.position[self.resnum_to_iac(0)])
        self.position = list(np.add(self.position, start))
        current_vec = np.subtract(self.position[3], self.position[1])
        for i in range(self.residues):
            for x in range(3):
                pos = np.subtract(function(3 * i + x), function(3 * i + x - 1))
                q = QuaternionBetween(pos, current_vec)
                current_vec = pos
                if x == 0:
                    for ind in self.res_indexes[i]:
                        self.position[ind] == q.orient(self.position[ind])
                if x == 1:
                    for ind in self.res_indexes[i][3:]:
                        self.position[ind] == q.orient(self.position[ind])
                if x == 2:
                    self.position[self.res_indexes[i][4]] == q.orient(self.position[self.res_indexes[i][4]])
                    for ind in self.res_indexes[i][7:]:
                        self.position[ind] == q.orient(self.position[ind])
                for arr in self.res_indexes[i + 1:]:
                    for index in arr:
                        self.position[index] = q.orient(self.position[index])

    def add_all(self, resname):

        for i in range(self.residues):
            self.add_residue(resname, self.resnum_to_iac(i))


    def get_subsection(self, res1, res2):
        new = ProteinBackbone(res2 - res1 + 1)
        for i in range(res1, res2 + 1, 1):
            data = [self.position[index] for index in self.res_indexes[i]]
            new.add_residue(self.res_names[i], new.resnum_to_iac(i - res1), data=data)
            if i in self.connections:
                new.connections.append(i)
        return new

    def make_connections_rigid(self):
        count = len(self.rigid_centers)
        for ind, c in enumerate(self.connections):
            self.make_rigid(c, body_tag=(count + ind))
            self.make_rigid(c - 1, body_tag=(count + ind))

    def make_rigid(self, res_index, body_tag=None):

        if body_tag is None:
            body_tag = self.rigid_centers[-1] + 1
        if body_tag not in self.rigid_centers:
            if not len(self.rigid_centers) == 0:
                if not body_tag > self.rigid_centers[-1]:
                    raise ValueError("body tag must be increasing", self.rigid_centers[-1], body_tag)
            self.rigid_centers.append(body_tag)
        if res_index >= len(self.res_indexes):
            raise ValueError("res_index too large", len(self.res_indexes), res_index)
        if res_index < 0:
            raise ValueError("negative res_index not allowed", res_index)
        for i in self.res_indexes[res_index]:
            self.body[i] = body_tag

    def make_rigids(self, start, end, body_tag=None):
        if body_tag is None:
            body_tag = len(self.rigid_centers)
        for i in range(start, end + 1, 1):
            self.make_rigid(i, body_tag)
        return body_tag

    def get_collection(self):

        return ProteinCollection([self])
