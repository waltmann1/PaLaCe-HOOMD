from __future__ import division
from ProteinChain import ProteinChain
from ProteinBackbone import ChainTooShortError
from ProteinBackbone import ProteinBackbone
from ProteinCollection import ProteinCollection
import hoomd
import Bio.PDB
import numpy as np
import numpy.linalg as la
from Quaternion import QuaternionBetween
import copy as cp
from Utils import all_amino_acids


class ProteinComplex(ProteinCollection):

    """
    Biopython has an idea of the protein being made of models and the models being made of chains,
    However the protein is normally made of one model, so this class will instantiate the first.
    Perhaps a di/tri/N - mer class could handle multiple models
    """
    def __init__(self, PDBName, denature=False):
        """

        :param PDBName: name of the pdb file representing the complex
        """

        self.chains = []
        super(ProteinComplex, self).__init__(self.chains)
        self.terminals = []
        self.rigid_count = 0
        self.name = PDBName[:-4]
        self.dump_context = None
        self.reinit_context = None
        self.PDBName = PDBName

        parser = Bio.PDB.PDBParser()
        self.struct = parser.get_structure(PDBName[:-4], PDBName)
        bbs, terminals, indexes = self.read_remark465(PDBName)
        list_470 = self.read_remark470(PDBName)
        #print(list_470)
        for ind, chain in enumerate(list_470):
                for ind2 in range(len(chain)):
                    chain[ind2] = int(chain[ind2] - np.sum([1 + ter[1] - ter[0] for ind3, ter in enumerate(terminals)
                                           if indexes[ind3] == ind and terminals[ind3][0] <= chain[ind2]]))
        #print(list_470, terminals, len(self.struct[0]))
        if len(self.struct[0]) < len(list_470):
            hit = False
            while len(self.struct[0]) < len(list_470) and not hit:
                if len(list_470[0]) == 0:
                    list_470.remove(list_470[0])
                else:
                    hit = True
        #print(list_470)
        for model in self.struct:
            for ind, chain in enumerate(model):
                try:
                    if ind + 1 <= len(list_470):
                        c = ProteinChain(chain, chain_index=ind, denature=denature, list_470=list_470[ind])
                    else:
                        c = ProteinChain(chain, chain_index=ind, denature=denature)
                    ters = self.main_chain_terminals(PDBName)[ind]
                    if ind in indexes:
                        ters = self.main_chain_terminals(PDBName)[ind]
                        #print(ters)
                        for count, bb in enumerate(bbs):
                            if indexes[count] == ind:
                                if terminals[count][0] - 1 == ters[1]:
                                    bb.align([0, 0, -1])
                                    # print(len(c.res_names), c.res_indexes[-1], c.position[-1])
                                    c.add_chain_to_C_terminal(bb)
                                    # print(len(c.res_names[-1]), c.res_indexes[-1], c.position[-1])
                                    print('added ' + str(bb.residues) + ' resiudes to C terminal of chain ' + str(ind))
                                    ters[1] = terminals[count][1]
                                if terminals[count][1] + 1 == ters[0]:
                                    bb.align([0, 0, 1])
                                    c.add_chain_to_N_terminal(bb)
                                    print('added ' + str(bb.residues) + ' resiudes to N terminal of chain ' + str(ind))
                                    ters[0] = terminals[count][0]
                    self.add_chain(c)
                    self.terminals.append(ters)
                    #c.dump_xyz(str(0) + str(ind) +'.xyz')
                except ChainTooShortError:
                    print("Chain with less than 1 residue being ignored")
        if PDBName == '1sva.pdb':
            self.remove_chain(self.chains[-1])

        if denature:
            for ind, chain in enumerate(self.chains):
                chain.shift([0, ind * 10, 0])

        self.num_particles = np.sum([chain.num_particles for chain in self.chains])

        self.indices, self.chas = self.get_ss_indices(PDBName, len(self.chains))

        self.ss_bonds = [[self.calc_particle_index(self.chas[h][0], ind[0]), self.calc_particle_index(self.chas[h][1], ind[1])]
                         for h, ind in enumerate(self.indices)]
        self.ss_bond_names = ['Sc-Sc' for _ in range(len(self.ss_bonds))]

        self.ss_angles = [[self.calc_particle_index(self.chas[h][0], ind[0] - 1),
                           self.calc_particle_index(self.chas[h][0], ind[0]),
                           self.calc_particle_index(self.chas[h][1], ind[1])] for h, ind in enumerate(self.indices)] + \
                         [[self.calc_particle_index(self.chas[h][1], ind[1] - 1),
                           self.calc_particle_index(self.chas[h][1], ind[1]),
                           self.calc_particle_index(self.chas[h][0], ind[0])] for h, ind in enumerate(self.indices)]

        self.ss_angle_names = ['betaC-Sc-Sc' for _ in range(len(self.ss_angles))]

        self.ss_dihedrals = [[self.calc_particle_index(self.chas[h][0], ind[0] - 1),
                              self.calc_particle_index(self.chas[h][0], ind[0]),
                              self.calc_particle_index(self.chas[h][1], ind[1]),
                              self.calc_particle_index(self.chas[h][1], ind[1] - 1)] for h, ind in enumerate(self.indices)]

        self.ss_dihedral_names = ['betaC-Sc-Sc-betaC' for _ in range(len(self.ss_dihedrals))]

        self.gsd_res_map = None



    def is_number(self, s):
        try:
            float(s)
            return True
        except ValueError:
            return False

    def find_disulphide_bonds(self, PDBName, num_chains):
        """

        :return: list of sulphide bonds as [[chain1, res1, chain2, res2], ...]
        """
        ss_bonds = []
        f = open(PDBName)
        data = f.readlines()
        for line in data:
            l = line.split() + [' ']
            if l[0] == 'SSBOND':
                chain1 = int(self.my_int(l[3]) - 1)
                res1 = int(l[4])
                chain2 = int(self.my_int(l[6]) - 1)
                res2 = int(l[7])

                ss_bonds.append([chain1, res1, chain2, res2])
        if len(ss_bonds) == 0:
            return None

        max_chain = np.max([[bond[0], bond[2]] for bond in ss_bonds])
        diff = max_chain - num_chains + 1
        if diff > 0:
            for bond in ss_bonds:
                bond[0] -= diff
                bond[2] -= diff
        return ss_bonds

    def get_ss_indices(self, PDBName, num_chains):

        indices = []
        chains = []
        ss_bonds = self.find_disulphide_bonds(PDBName, num_chains)
        if ss_bonds is None:
            return indices, chains
        for bond in ss_bonds:
            index_1 = self.chains[bond[0]].res_indexes[bond[1] - self.terminals[bond[0]][0]][-1]
            index_2 = self.chains[bond[2]].res_indexes[bond[3] - self.terminals[bond[2]][0]][-1]
            indices.append([index_1, index_2])
            chains.append([bond[0], bond[2]])
            if not self.chains[bond[0]].type[index_1] == 'Sc' and self.chains[bond[2]].type[index_2] == 'Sc':
                raise ValueError('disuplhihde bond indices not found correctly on chain' + str(bond[0]) + ' and ' +
                                 str(bond[2]) + 'for indices ' + str(index_1) + ' and' + str(index_2) + ' found types ',
                              self.chains[bond[0]].type[index_1], self.chains[bond[2]].type[index_2])
        return indices, chains

    def read_remark470(self, PDBName):

        f = open(PDBName)
        data = f.readlines()
        aa = all_amino_acids()
        chains = []
        chain_indices = []
        chain_info = []
        for line in data:
            s = line.split() + [' '] * 2
            if s[0] == 'REMARK' and s[1] == '470':
                if s[2] in aa:
                    chain_info.append([s[2], int(self.my_int(s[3]) - 1), int(self.my_int(s[4]) - 1)])
        for ind, x in enumerate(chain_info):
            if x[1] not in chain_indices:
                chain_indices.append(x[1])
                chains.append([x])
            else:
                chains[chain_indices.index(x[1])].append(x)
        hit_indexes = list(set([res[1] for chain in chains for res in chain]))
        if len(hit_indexes) != 0:
            v = int(max(hit_indexes) + 1)
        else:
            v = 0
        final = [[] for _ in range(v)]
        for chain in chains:
            for res in chain:
                final[res[1]].append(res[2])
        return final


    def read_remark465(self, PDBName):

        f = open(PDBName)
        data = f.readlines()
        chains = []
        chain_indices = []
        chain_info = []
        for line in data:
            s = line.split() + [' '] * 2
            if s[0] == 'REMARK' and s[1] == '465':
                try:
                    imp = __import__(s[2], fromlist=[''])
                    chain_info.append([s[2], self.my_int(s[3]), self.my_int(s[4])])
                except ImportError:
                    pass
        for ind, x in enumerate(chain_info):
            if x[1] not in chain_indices:
                chain_indices.append(x[1])
                chains.append([x])
            else:
                chains[chain_indices.index(x[1])].append(x)
        bbs = []
        # print(chains)
        terminals = []
        indexes = []
        for chain_ind, part in enumerate(chains):
            length = 0
            for ind, res in enumerate(part):
                if ind == 0:
                    length += 1
                elif self.my_int(self.my_int(res[2]) - 1) != part[ind-1][2]:
                    bb = ProteinBackbone(length)
                    for ind2, res2 in enumerate(part[ind - length:ind]):
                        bb.add_residue(res2[0], bb.resnum_to_iac(ind2))
                    bbs.append(bb)
                    terminals.append([part[ind-length][2], part[ind - 1][2]])
                    indexes.append(self.my_int(chain_indices[chain_ind] - 1))
                    length = 1
                else:
                    length += 1
            bb = ProteinBackbone(length)
            for ind3, info in enumerate(part[-length:]):
                bb.add_residue(info[0], bb.resnum_to_iac(ind3))
            bbs.append(bb)
            terminals.append([part[ind + 1 - length][2], part[ind][2]])
            # print([part[ind + 1- length][2], part[ind - 1][2]])
            indexes.append(self.my_int(chain_indices[chain_ind] - 1))

        return bbs, terminals, indexes

    def main_chain_terminals(self, PDBName):
        """

        :param PDBname:
        :return:
        """
        f = open(PDBName)
        data = f.readlines()
        terminals = [[]]
        chain_index = 0
        res_index = -10000
        aa = all_amino_acids()
        go = False
        for line in data:
            s = line.split() + [' ']
            if s[0] == 'ATOM':
                if len(s[2]) > 3:
                    if s[2][-3:] in aa:
                        s.insert(3, s[2][3:])
                        s[2] = s[2][:3]
                if len(s[4]) > 3:
                    s.insert(5, s[4][2:])
                    s[4] = s[4][0]
                if not go:
                    go = True
                    #terminals.append([self.my_int(s[5])])
                    self.add_to_terminal(self.my_int(s[4]) - 1, self.my_int(s[5]), terminals)
                elif self.my_int(s[4]) != chain_index:
                    #terminals.append([res_index])
                    self.add_to_terminal(chain_index - 1, self.my_int(s[5]), terminals)
                elif res_index != self.my_int(s[5]) and res_index != self.my_int(self.my_int(s[5]) - 1):
                    #terminals[chain_index - 1].append(res_index)
                    self.add_to_terminal(chain_index - 1, self.my_int(s[5]), terminals)
                else:
                    pass
                res_index = self.my_int(s[5])
                chain_index = self.my_int(s[4])
            elif s[0] != 'ATOM' and go:
                terminals[chain_index - 1].append(res_index)
                go = False
            else:
                pass

        if len(self.struct[0]) < len(terminals):
            hit = False
            while len(self.struct[0]) < len(terminals) and not hit:
                if len(terminals[0]) == 0:
                    terminals.remove(terminals[0])
                else:
                    hit = True
        return terminals

    def add_to_terminal(self, int, to_add, terminals):
        #int is chain index starting at 0 or chain_index - 1

        if len(terminals) == int:
            terminals.append([to_add])
        elif len(terminals) < int:
            while len(terminals) < int:
                terminals.append([])
            terminals.append([to_add])
        else:
            terminals[int].append(to_add)

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

    def calc_particle_index(self, chain, index):

        pi = 0
        for i in range(chain):
            pi += self.chains[i].num_particles
        return pi + index


    def remove_remark_465(self):

        mc = self.main_chain_terminals(self.PDBName)
        subs = [[chain[0]-self.terminals[ind][0], chain[1]-self.terminals[ind][1]] for ind, chain in enumerate(mc)]
        []
        for ind, chain in enumerate(self.chains):
            self.remove_chain(chain)
            self.add_chain(chain.get_subsection(subs[ind][0], subs[ind[1]]))

