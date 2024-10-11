import unittest
import shutil
import os
import numpy as np
from ProteinBackbone import ProteinBackbone


class TestResidues(unittest.TestCase):

    def setUp(self):
        self.dir = "dumpfolder/"
        directory = os.path.dirname(self.dir)
        if not os.path.exists(directory):
            os.makedirs(directory)

    def tearDown(self):
        shutil.rmtree(self.dir)

    def test_ALA(self):
        pb = ProteinBackbone(1)
        pb.add_residue('ALA', 2)
        bonds = pb.bonds[2:]
        angles = pb.angles[1:]
        dihedrals = pb.dihedrals
        dc = pb.dihedrals_coupling
        types = pb.type[5:]
        guess_bonds = [[2, 5]]
        guess_angles = [[0, 2, 5], [3, 2, 5]]
        guess_types = ['Sa']
        guess_dihedrals = []
        guess_dc = []
        self.check_group(bonds, guess_bonds)
        self.check_group(angles, guess_angles)
        self.check_group(dc, guess_dc)
        self.check_group(dihedrals, guess_dihedrals)
        self.check_single(types, guess_types)
        self.assertEquals(pb.charge[-1], 0)

    def test_ARG(self):
        pb = ProteinBackbone(1)
        pb.add_residue('ARG', 2)
        bonds = pb.bonds[2:]
        bond_names = pb.bond_names[2:]
        angles = pb.angles[1:]
        angle_names = pb.angle_names[1:]
        dihedrals = pb.dihedrals
        dihedral_names = pb.dihedral_names
        dc = pb.dihedrals_coupling
        dc_names = pb.dihedral_coupling_names
        types = pb.type[5:]
        guess_bonds = [[2, 5], [6, 5], [6, 2], [7, 6], [7, 5], [8, 7], [9, 8], [10, 9]]
        guess_bond_names = ['alphaC-betaC', 'halfway', 'exclusion', 'halfway', 'betaC-gammaCr', 'gammaCr-deltaCr',
                            'deltaCr-epsilonNr', 'epsilonNr-Tr']
        guess_angles = [[0, 2, 5], [3, 2, 5], [2, 5, 7], [5, 7, 8], [7, 8, 9], [8, 9, 10]]
        guess_angle_names = ['N-alphaC-betaC', 'betaC-alphaC-C', 'alphaC-betaC-gammaCr', 'betaC-gammaCr-deltaCr',
                             'gammaCr-deltaCr-epsilonNr', 'deltaCr-epsilonNr-Tr']
        guess_types = ['betaC', 'Sr', 'gammaCr', 'deltaCr', 'epsilonNr', 'Tr']
        guess_dihedrals = []
        guess_dihedral_names = []
        guess_dc = [[0, 2, 5, 7, 2, 5, 7, 8], [5, 7, 8, 9, 7, 8, 9, 10]]
        guess_dc_names = ['ARG', 'xARG']
        self.check_group(bonds, guess_bonds)
        self.check_single(bond_names, guess_bond_names)
        self.check_group(angles, guess_angles)
        self.check_single(angle_names, guess_angle_names)
        self.check_group(dihedrals, guess_dihedrals)
        self.check_group(dihedral_names, guess_dihedral_names)
        self.check_group(dc, guess_dc)
        self.check_group(dc_names, guess_dc_names)
        self.check_single(types, guess_types)
        self.assertEquals(pb.charge[-1], 1)

    def test_ASN(self):

        pb = ProteinBackbone(1)
        pb.add_residue('ASN', 2)
        bonds = pb.bonds[2:]
        bond_names = pb.bond_names[2:]
        angles = pb.angles[1:]
        angle_names = pb.angle_names[1:]
        dihedrals = pb.dihedrals
        dihedral_names = pb.dihedral_names
        dc = pb.dihedrals_coupling
        dc_names = pb.dihedral_coupling_names
        types = pb.type[5:]
        guess_bonds = [[2, 5], [5, 6]]
        guess_bond_names = ['alphaC-betaC', 'betaC-Sn']
        guess_angles = [[0, 2, 5], [3, 2, 5], [2, 5, 6]]
        guess_angle_names = ['N-alphaC-betaC', 'betaC-alphaC-C', 'alphaC-betaC-Sn']
        guess_types = ['betaC', 'Sn']
        guess_dihedrals = [[0, 2, 5, 6]]
        guess_dihedral_names = ['N-alphaC-betaC-Sn']
        guess_dc = []
        guess_dc_names = []
        self.check_group(bonds, guess_bonds)
        self.check_single(bond_names, guess_bond_names)
        self.check_group(angles, guess_angles)
        self.check_single(angle_names, guess_angle_names)
        self.check_group(dihedrals, guess_dihedrals)
        self.check_group(dihedral_names, guess_dihedral_names)
        self.check_group(dc, guess_dc)
        self.check_group(dc_names, guess_dc_names)
        self.check_single(types, guess_types)
        self.assertEquals(pb.charge[-1], 0)

    def test_ASP(self):

        pb = ProteinBackbone(1)
        pb.add_residue('ASP', 2)
        bonds = pb.bonds[2:]
        bond_names = pb.bond_names[2:]
        angles = pb.angles[1:]
        angle_names = pb.angle_names[1:]
        dihedrals = pb.dihedrals
        dihedral_names = pb.dihedral_names
        dc = pb.dihedrals_coupling
        dc_names = pb.dihedral_coupling_names
        types = pb.type[5:]
        guess_bonds = [[2, 5], [5, 6]]
        guess_bond_names = ['alphaC-betaC', 'betaC-Sd']
        guess_angles = [[0, 2, 5], [3, 2, 5], [2, 5, 6]]
        guess_angle_names = ['N-alphaC-betaC', 'betaC-alphaC-C', 'alphaC-betaC-Sd']
        guess_types = ['betaC', 'Sd']
        guess_dihedrals = [[0, 2, 5, 6]]
        guess_dihedral_names = ['N-alphaC-betaC-Sd']
        guess_dc = []
        guess_dc_names = []
        self.check_group(bonds, guess_bonds)
        self.check_single(bond_names, guess_bond_names)
        self.check_group(angles, guess_angles)
        self.check_single(angle_names, guess_angle_names)
        self.check_group(dihedrals, guess_dihedrals)
        self.check_group(dihedral_names, guess_dihedral_names)
        self.check_group(dc, guess_dc)
        self.check_group(dc_names, guess_dc_names)
        self.check_single(types, guess_types)
        self.assertEquals(pb.charge[-1], -1)

    def test_CYS(self):

        pb = ProteinBackbone(1)
        pb.add_residue('CYS', 2)
        bonds = pb.bonds[2:]
        bond_names = pb.bond_names[2:]
        angles = pb.angles[1:]
        angle_names = pb.angle_names[1:]
        dihedrals = pb.dihedrals
        dihedral_names = pb.dihedral_names
        dc = pb.dihedrals_coupling
        dc_names = pb.dihedral_coupling_names
        types = pb.type[5:]
        guess_bonds = [[2, 5], [5, 6]]
        guess_bond_names = ['alphaC-betaC', 'betaC-Sc']
        guess_angles = [[0, 2, 5], [3, 2, 5], [2, 5, 6]]
        guess_angle_names = ['N-alphaC-betaC', 'betaC-alphaC-C', 'alphaC-betaC-Sc']
        guess_types = ['betaC', 'Sc']
        guess_dihedrals = [[0, 2, 5, 6]]
        guess_dihedral_names = ['N-alphaC-betaC-Sc']
        guess_dc = []
        guess_dc_names = []
        self.check_group(bonds, guess_bonds)
        self.check_single(bond_names, guess_bond_names)
        self.check_group(angles, guess_angles)
        self.check_single(angle_names, guess_angle_names)
        self.check_group(dihedrals, guess_dihedrals)
        self.check_group(dihedral_names, guess_dihedral_names)
        self.check_group(dc, guess_dc)
        self.check_group(dc_names, guess_dc_names)
        self.check_single(types, guess_types)
        self.assertEquals(pb.charge[-1], 0)

    def test_GLN(self):

        pb = ProteinBackbone(1)
        pb.add_residue('GLN', 2)
        bonds = pb.bonds[2:]
        bond_names = pb.bond_names[2:]
        angles = pb.angles[1:]
        angle_names = pb.angle_names[1:]
        dihedrals = pb.dihedrals
        dihedral_names = pb.dihedral_names
        dc = pb.dihedrals_coupling
        dc_names = pb.dihedral_coupling_names
        types = pb.type[5:]
        guess_bonds = [[2, 5], [5, 6], [6, 7], [7, 2]]
        guess_bond_names = ['alphaC-betaC', 'betaC-Sq', 'Sq-Tq', 'exclusion']
        guess_angles = [[0, 2, 5], [3, 2, 5], [2, 5, 6]]
        guess_angle_names = ['N-alphaC-betaC', 'betaC-alphaC-C', 'alphaC-betaC-Sq']
        guess_types = ['betaC', 'Sq', 'Tq']
        guess_dihedrals = []
        guess_dihedral_names = []
        guess_dc = [[0, 2, 5, 6, 2, 5, 6, 7]]
        guess_dc_names = ['GLN']
        self.check_group(bonds, guess_bonds)
        self.check_single(bond_names, guess_bond_names)
        self.check_group(angles, guess_angles)
        self.check_single(angle_names, guess_angle_names)
        self.check_group(dihedrals, guess_dihedrals)
        self.check_group(dihedral_names, guess_dihedral_names)
        self.check_group(dc, guess_dc)
        self.check_group(dc_names, guess_dc_names)
        self.check_single(types, guess_types)
        self.assertEquals(pb.charge[-1], 0)

    def test_GLU(self):

        pb = ProteinBackbone(1)
        pb.add_residue('GLU', 2)
        bonds = pb.bonds[2:]
        bond_names = pb.bond_names[2:]
        angles = pb.angles[1:]
        angle_names = pb.angle_names[1:]
        dihedrals = pb.dihedrals
        dihedral_names = pb.dihedral_names
        dc = pb.dihedrals_coupling
        dc_names = pb.dihedral_coupling_names
        types = pb.type[5:]
        guess_bonds = [[2, 5], [5, 6], [6, 7], [7, 2]]
        guess_bond_names = ['alphaC-betaC', 'betaC-Se', 'Se-Te', 'exclusion']
        guess_angles = [[0, 2, 5], [3, 2, 5], [2, 5, 6]]
        guess_angle_names = ['N-alphaC-betaC', 'betaC-alphaC-C', 'alphaC-betaC-Se']
        guess_types = ['betaC', 'Se', 'Te']
        guess_dihedrals = []
        guess_dihedral_names = []
        guess_dc = [[0, 2, 5, 6, 2, 5, 6, 7]]
        guess_dc_names = ['GLU']
        self.check_group(bonds, guess_bonds)
        self.check_single(bond_names, guess_bond_names)
        self.check_group(angles, guess_angles)
        self.check_single(angle_names, guess_angle_names)
        self.check_group(dihedrals, guess_dihedrals)
        self.check_group(dihedral_names, guess_dihedral_names)
        self.check_group(dc, guess_dc)
        self.check_group(dc_names, guess_dc_names)
        self.check_single(types, guess_types)
        self.assertEquals(pb.charge[-1], -1)

    def test_GLY(self):

        pb = ProteinBackbone(1)
        pb.add_residue('GLY', 2)
        bonds = pb.bonds[2:]
        bond_names = pb.bond_names[2:]
        angles = pb.angles[1:]
        angle_names = pb.angle_names[1:]
        dihedrals = pb.dihedrals
        dihedral_names = pb.dihedral_names
        dc = pb.dihedrals_coupling
        dc_names = pb.dihedral_coupling_names
        types = pb.type[5:]
        guess_bonds = []
        guess_bond_names = []
        guess_angles = []
        guess_angle_names = []
        guess_types = []
        guess_dihedrals = []
        guess_dihedral_names = []
        guess_dc = []
        guess_dc_names = []
        self.check_group(bonds, guess_bonds)
        self.check_single(bond_names, guess_bond_names)
        self.check_group(angles, guess_angles)
        self.check_single(angle_names, guess_angle_names)
        self.check_group(dihedrals, guess_dihedrals)
        self.check_group(dihedral_names, guess_dihedral_names)
        self.check_group(dc, guess_dc)
        self.check_group(dc_names, guess_dc_names)
        self.check_single(types, guess_types)
        self.assertEquals('alphaCg', pb.type[2])

    def test_HIS(self):

        pb = ProteinBackbone(1)
        pb.add_residue('HIS', 2)
        bonds = pb.bonds[2:]
        bond_names = pb.bond_names[2:]
        angles = pb.angles[1:]
        angle_names = pb.angle_names[1:]
        dihedrals = pb.dihedrals
        dihedral_names = pb.dihedral_names
        dc = pb.dihedrals_coupling
        dc_names = pb.dihedral_coupling_names
        types = pb.type[5:]
        guess_bonds = [[2, 5], [5, 6]]
        guess_bond_names = ['alphaC-betaC', 'Sh-Th']
        guess_angles = [[0, 2, 5], [3, 2, 5], [2, 5, 6]]
        guess_angle_names = ['N-alphaC-betaC', 'betaC-alphaC-C', 'alphaC-betaC-Th']
        guess_types = ['Sh', 'Th']
        guess_dihedrals = [[0, 2, 5, 6]]
        guess_dihedral_names = ['N-alphaC-betaC-Th']
        guess_dc = []
        guess_dc_names = []
        self.check_group(bonds, guess_bonds)
        self.check_single(bond_names, guess_bond_names)
        self.check_group(angles, guess_angles)
        self.check_single(angle_names, guess_angle_names)
        self.check_group(dihedrals, guess_dihedrals)
        self.check_group(dihedral_names, guess_dihedral_names)
        self.check_group(dc, guess_dc)
        self.check_group(dc_names, guess_dc_names)
        self.check_single(types, guess_types)
        self.assertEquals(pb.charge[-1], 0)

    def test_ILE(self):

        pb = ProteinBackbone(1)
        pb.add_residue('ILE', 2)
        bonds = pb.bonds[2:]
        bond_names = pb.bond_names[2:]
        angles = pb.angles[1:]
        angle_names = pb.angle_names[1:]
        dihedrals = pb.dihedrals
        dihedral_names = pb.dihedral_names
        dc = pb.dihedrals_coupling
        dc_names = pb.dihedral_coupling_names
        types = pb.type[5:]
        guess_bonds = [[2, 5], [5, 6], [5, 7], [6, 8], [8, 9], [9, 2], [5, 9]]
        guess_bond_names = ['alphaC-betaC', 'betaC-gammaC1i', 'betaC-gammaC2i', 'gammaC1i-deltaCi', 'Si-1', 'exclusion','Si-1']
        guess_angles = [[0, 2, 5], [3, 2, 5], [2, 5, 6], [2, 5, 7], [6, 5, 7], [5, 6, 8]]
        guess_angle_names = ['N-alphaC-betaC', 'betaC-alphaC-C', 'alphaC-betaC-gammaC1i', 'alphaC-betaC-gammaC2i',
                             'gammaC1i-betaC-gammaC2i', 'betaC-gammaC1i-deltaCi']
        guess_types = ['betaC', 'gammaC1i', 'gammaC2i', 'deltaCi', 'Si']
        guess_dihedrals = []
        guess_dihedral_names = []
        guess_dc = [[0, 2, 5, 6, 2, 5, 6, 8]]
        guess_dc_names = ['ILE']
        self.check_group(bonds, guess_bonds)
        self.check_single(bond_names, guess_bond_names)
        self.check_group(angles, guess_angles)
        self.check_single(angle_names, guess_angle_names)
        self.check_group(dihedrals, guess_dihedrals)
        self.check_group(dihedral_names, guess_dihedral_names)
        self.check_group(dc, guess_dc)
        self.check_group(dc_names, guess_dc_names)
        self.check_single(types, guess_types)
        self.assertEquals(pb.charge[-1], 0)

    def test_LEU(self):

        pb = ProteinBackbone(1)
        pb.add_residue('LEU', 2)
        bonds = pb.bonds[2:]
        bond_names = pb.bond_names[2:]
        angles = pb.angles[1:]
        angle_names = pb.angle_names[1:]
        dihedrals = pb.dihedrals
        dihedral_names = pb.dihedral_names
        dc = pb.dihedrals_coupling
        dc_names = pb.dihedral_coupling_names
        types = pb.type[5:]
        guess_bonds = [[2, 5], [5, 6]]
        guess_bond_names = ['alphaC-betaC', 'betaC-Sl']
        guess_angles = [[0, 2, 5], [3, 2, 5], [2, 5, 6]]
        guess_angle_names = ['N-alphaC-betaC', 'betaC-alphaC-C', 'alphaC-betaC-Sl']
        guess_types = ['betaC', 'Sl']
        guess_dihedrals = [[0, 2, 5, 6]]
        guess_dihedral_names = ['N-alphaC-betaC-Sl']
        guess_dc = []
        guess_dc_names = []
        self.check_group(bonds, guess_bonds)
        self.check_single(bond_names, guess_bond_names)
        self.check_group(angles, guess_angles)
        self.check_single(angle_names, guess_angle_names)
        self.check_group(dihedrals, guess_dihedrals)
        self.check_group(dihedral_names, guess_dihedral_names)
        self.check_group(dc, guess_dc)
        self.check_group(dc_names, guess_dc_names)
        self.check_single(types, guess_types)
        self.assertEquals(pb.charge[-1], 0)

    def test_LYS(self):
        pb = ProteinBackbone(1)
        pb.add_residue('LYS', 2)
        bonds = pb.bonds[2:]
        bond_names = pb.bond_names[2:]
        angles = pb.angles[1:]
        angle_names = pb.angle_names[1:]
        dihedrals = pb.dihedrals
        dihedral_names = pb.dihedral_names
        dc = pb.dihedrals_coupling
        dc_names = pb.dihedral_coupling_names
        types = pb.type[5:]
        guess_bonds = [[2, 5], [6, 5], [6, 2], [7, 5], [7, 6], [8, 7], [9, 8], [10, 9]]
        guess_bond_names = ['alphaC-betaC', 'halfway', 'exclusion', 'betaC-gammaCk', 'halfway', 'gammaCk-deltaCk',
                            'deltaCk-epsilonCk', 'epsilonCk-Tk']
        guess_angles = [[0, 2, 5], [3, 2, 5], [2, 5, 7], [5, 7, 8], [7, 8, 9], [8, 9, 10]]
        guess_angle_names = ['N-alphaC-betaC', 'betaC-alphaC-C', 'alphaC-betaC-gammaCk', 'betaC-gammaCk-deltaCk',
                             'gammaCk-deltaCk-epsilonCk', 'deltaCk-epsilonCk-Tk']
        guess_types = ['betaC', 'Sk', 'gammaCk', 'deltaCk', 'epsilonCk', 'Tk']
        guess_dihedrals = []
        guess_dihedral_names = []
        guess_dc = [[0, 2, 5, 7, 2, 5, 7, 8], [5, 7, 8, 9, 7, 8, 9, 10]]
        guess_dc_names = ['LYS', 'xLYS']
        self.check_group(bonds, guess_bonds)
        self.check_single(bond_names, guess_bond_names)
        self.check_group(angles, guess_angles)
        self.check_single(angle_names, guess_angle_names)
        self.check_group(dihedrals, guess_dihedrals)
        self.check_group(dihedral_names, guess_dihedral_names)
        self.check_group(dc, guess_dc)
        self.check_group(dc_names, guess_dc_names)
        self.check_single(types, guess_types)
        self.assertEquals(pb.charge[-1], 1)

    def test_MET(self):

        pb = ProteinBackbone(1)
        pb.add_residue('MET', 2)
        bonds = pb.bonds[2:]
        bond_names = pb.bond_names[2:]
        angles = pb.angles[1:]
        angle_names = pb.angle_names[1:]
        dihedrals = pb.dihedrals
        dihedral_names = pb.dihedral_names
        dc = pb.dihedrals_coupling
        dc_names = pb.dihedral_coupling_names
        types = pb.type[5:]
        guess_bonds = [[2, 5], [6, 5], [6, 2], [7, 6], [7, 5], [8, 7], [9, 8], [9, 6]]
        guess_bond_names = ['alphaC-betaC', 'halfway', 'exclusion', 'halfway', 'betaC-gammaCm', 'gammaCm-deltaSm',
                            'deltaSm-Tm', 'exclusion']
        guess_angles = [[0, 2, 5], [3, 2, 5], [2, 5, 7], [5, 7, 8], [7, 8, 9]]
        guess_angle_names = ['N-alphaC-betaC', 'betaC-alphaC-C', 'alphaC-betaC-gammaCm', 'betaC-gammaCm-deltaSm',
                             'gammaCm-deltaSm-Tm']
        guess_types = ['betaC', 'Sm', 'gammaCm', 'deltaSm', 'Tm']
        guess_dihedrals = [[5, 7, 8, 9]]
        guess_dihedral_names = ['betaC-gammaCm-deltaSm-Tm']
        guess_dc = [[0, 2, 5, 7, 2, 5, 7, 8]]
        guess_dc_names = ['MET']
        self.check_group(bonds, guess_bonds)
        self.check_single(bond_names, guess_bond_names)
        self.check_group(angles, guess_angles)
        self.check_single(angle_names, guess_angle_names)
        self.check_group(dihedrals, guess_dihedrals)
        self.check_single(dihedral_names, guess_dihedral_names)
        self.check_group(dc, guess_dc)
        self.check_single(dc_names, guess_dc_names)
        self.check_single(types, guess_types)
        self.assertEquals(pb.charge[-1], 0)

    def test_PHE(self):

        pb = ProteinBackbone(1)
        pb.add_residue('PHE', 2)
        bonds = pb.bonds[2:]
        bond_names = pb.bond_names[2:]
        angles = pb.angles[1:]
        angle_names = pb.angle_names[1:]
        dihedrals = pb.dihedrals
        dihedral_names = pb.dihedral_names
        dc = pb.dihedrals_coupling
        dc_names = pb.dihedral_coupling_names
        types = pb.type[5:]
        guess_bonds = [[2, 5], [5, 6]]
        guess_bond_names = ['alphaC-betaC', 'Sf-Tf']
        guess_angles = [[0, 2, 5], [3, 2, 5], [2, 5, 6]]
        guess_angle_names = ['N-alphaC-betaC', 'betaC-alphaC-C', 'alphaC-betaC-Tf']
        guess_types = ['Sf', 'Tf']
        guess_dihedrals = [[0, 2, 5, 6]]
        guess_dihedral_names = ['N-alphaC-betaC-Tf']
        guess_dc = []
        guess_dc_names = []
        self.check_group(bonds, guess_bonds)
        self.check_single(bond_names, guess_bond_names)
        self.check_group(angles, guess_angles)
        self.check_single(angle_names, guess_angle_names)
        self.check_group(dihedrals, guess_dihedrals)
        self.check_group(dihedral_names, guess_dihedral_names)
        self.check_group(dc, guess_dc)
        self.check_group(dc_names, guess_dc_names)
        self.check_single(types, guess_types)
        self.assertEquals(pb.charge[-1], 0)

    def test_PRO(self):

        pb = ProteinBackbone(1)
        pb.add_residue('PRO', 2)
        bonds = pb.bonds[2:]
        bond_names = pb.bond_names[2:]
        angles = pb.angles[1:]
        angle_names = pb.angle_names[1:]
        dihedrals = pb.dihedrals
        dihedral_names = pb.dihedral_names
        dc = pb.dihedrals_coupling
        dc_names = pb.dihedral_coupling_names
        types = pb.type[5:]
        guess_bonds = [[2, 5], [6, 5], [6, 2]]
        guess_bond_names = ['alphaC-betaC', 'betaC-Sp', 'exclusion']
        guess_angles = [[0, 2, 5], [3, 2, 5]]
        guess_angle_names = ['N-alphaC-betaC', 'betaC-alphaC-C']
        guess_types = ['betaC', 'Sp']
        guess_dihedrals = [[0, 2, 5, 6]]
        guess_dihedral_names = ['N-alphaC-betaC-Sp']
        guess_dc = []
        guess_dc_names = []
        self.check_group(bonds, guess_bonds)
        self.check_single(bond_names, guess_bond_names)
        self.check_group(angles, guess_angles)
        self.check_single(angle_names, guess_angle_names)
        self.check_group(dihedrals, guess_dihedrals)
        self.check_group(dihedral_names, guess_dihedral_names)
        self.check_group(dc, guess_dc)
        self.check_group(dc_names, guess_dc_names)
        self.check_single(types, guess_types)
        self.assertEquals(pb.charge[-1], 0)

    def test_SER(self):

        pb = ProteinBackbone(1)
        pb.add_residue('SER', 2)
        bonds = pb.bonds[2:]
        bond_names = pb.bond_names[2:]
        angles = pb.angles[1:]
        angle_names = pb.angle_names[1:]
        dihedrals = pb.dihedrals
        dihedral_names = pb.dihedral_names
        dc = pb.dihedrals_coupling
        dc_names = pb.dihedral_coupling_names
        types = pb.type[5:]
        guess_bonds = [[2, 5], [5, 6]]
        guess_bond_names = ['alphaC-betaC', 'betaC-Ss']
        guess_angles = [[0, 2, 5], [3, 2, 5], [2, 5, 6]]
        guess_angle_names = ['N-alphaC-betaC', 'betaC-alphaC-C', 'alphaC-betaC-Ss']
        guess_types = ['betaC', 'Ss']
        guess_dihedrals = [[0, 2, 5, 6]]
        guess_dihedral_names = ['N-alphaC-betaC-Ss']
        guess_dc = []
        guess_dc_names = []
        self.check_group(bonds, guess_bonds)
        self.check_single(bond_names, guess_bond_names)
        self.check_group(angles, guess_angles)
        self.check_single(angle_names, guess_angle_names)
        self.check_group(dihedrals, guess_dihedrals)
        self.check_group(dihedral_names, guess_dihedral_names)
        self.check_group(dc, guess_dc)
        self.check_group(dc_names, guess_dc_names)
        self.check_single(types, guess_types)
        self.assertEquals(pb.charge[-1], 0)

    def test_THR(self):

        pb = ProteinBackbone(1)
        pb.add_residue('THR', 2)
        bonds = pb.bonds[2:]
        bond_names = pb.bond_names[2:]
        angles = pb.angles[1:]
        angle_names = pb.angle_names[1:]
        dihedrals = pb.dihedrals
        dihedral_names = pb.dihedral_names
        dc = pb.dihedrals_coupling
        dc_names = pb.dihedral_coupling_names
        types = pb.type[5:]
        guess_bonds = [[2, 5], [5, 6]]
        guess_bond_names = ['alphaC-betaC', 'betaC-St']
        guess_angles = [[0, 2, 5], [3, 2, 5], [2, 5, 6]]
        guess_angle_names = ['N-alphaC-betaC', 'betaC-alphaC-C', 'alphaC-betaC-St']
        guess_types = ['betaC', 'St']
        guess_dihedrals = [[0, 2, 5, 6]]
        guess_dihedral_names = ['N-alphaC-betaC-St']
        guess_dc = []
        guess_dc_names = []
        self.check_group(bonds, guess_bonds)
        self.check_single(bond_names, guess_bond_names)
        self.check_group(angles, guess_angles)
        self.check_single(angle_names, guess_angle_names)
        self.check_group(dihedrals, guess_dihedrals)
        self.check_group(dihedral_names, guess_dihedral_names)
        self.check_group(dc, guess_dc)
        self.check_group(dc_names, guess_dc_names)
        self.check_single(types, guess_types)
        self.assertEquals(pb.charge[-1], 0)

    def test_TRP(self):

        pb = ProteinBackbone(1)
        pb.add_residue('TRP', 2)
        bonds = pb.bonds[2:]
        bond_names = pb.bond_names[2:]
        angles = pb.angles[1:]
        angle_names = pb.angle_names[1:]
        dihedrals = pb.dihedrals
        dihedral_names = pb.dihedral_names
        dc = pb.dihedrals_coupling
        dc_names = pb.dihedral_coupling_names
        types = pb.type[5:]
        guess_bonds = [[2, 5], [5, 6], [7, 2], [6, 7]]
        guess_bond_names = ['alphaC-betaC', 'betaC-Sw', 'exclusion', 'Sw-Tw']
        guess_angles = [[0, 2, 5], [3, 2, 5], [2, 5, 6]]
        guess_angle_names = ['N-alphaC-betaC', 'betaC-alphaC-C', 'alphaC-betaC-Sw']
        guess_types = ['betaC', 'Sw', 'Tw']
        guess_dihedrals = []
        guess_dihedral_names = []
        guess_dc = [[0, 2, 5, 6, 2, 5, 6, 7]]
        guess_dc_names = ['TRP']
        self.check_group(bonds, guess_bonds)
        self.check_single(bond_names, guess_bond_names)
        self.check_group(angles, guess_angles)
        self.check_single(angle_names, guess_angle_names)
        self.check_group(dihedrals, guess_dihedrals)
        self.check_group(dihedral_names, guess_dihedral_names)
        self.check_group(dc, guess_dc)
        self.check_group(dc_names, guess_dc_names)
        self.check_single(types, guess_types)
        self.assertEquals(pb.charge[-1], 0)

    def test_TYR(self):

        pb = ProteinBackbone(1)
        pb.add_residue('TYR', 2)
        bonds = pb.bonds[2:]
        bond_names = pb.bond_names[2:]
        angles = pb.angles[1:]
        angle_names = pb.angle_names[1:]
        dihedrals = pb.dihedrals
        dihedral_names = pb.dihedral_names
        dc = pb.dihedrals_coupling
        dc_names = pb.dihedral_coupling_names
        types = pb.type[5:]
        guess_bonds = [[2, 5], [5, 6]]
        guess_bond_names = ['alphaC-betaC', 'Sy-Ty']
        guess_angles = [[0, 2, 5], [3, 2, 5], [2, 5, 6]]
        guess_angle_names = ['N-alphaC-betaC', 'betaC-alphaC-C', 'alphaC-betaC-Ty']
        guess_types = ['Sy', 'Ty']
        guess_dihedrals = [[0, 2, 5, 6]]
        guess_dihedral_names = ['N-alphaC-betaC-Ty']
        guess_dc = []
        guess_dc_names = []
        self.check_group(bonds, guess_bonds)
        self.check_single(bond_names, guess_bond_names)
        self.check_group(angles, guess_angles)
        self.check_single(angle_names, guess_angle_names)
        self.check_group(dihedrals, guess_dihedrals)
        self.check_group(dihedral_names, guess_dihedral_names)
        self.check_group(dc, guess_dc)
        self.check_group(dc_names, guess_dc_names)
        self.check_single(types, guess_types)
        self.assertEquals(pb.charge[-1], 0)

    def test_VAL(self):

        pb = ProteinBackbone(1)
        pb.add_residue('VAL', 2)
        bonds = pb.bonds[2:]
        bond_names = pb.bond_names[2:]
        angles = pb.angles[1:]
        angle_names = pb.angle_names[1:]
        dihedrals = pb.dihedrals
        dihedral_names = pb.dihedral_names
        dc = pb.dihedrals_coupling
        dc_names = pb.dihedral_coupling_names
        types = pb.type[5:]
        guess_bonds = [[2, 5], [5, 6]]
        guess_bond_names = ['alphaC-betaC', 'betaC-Sv']
        guess_angles = [[0, 2, 5], [3, 2, 5], [2, 5, 6]]
        guess_angle_names = ['N-alphaC-betaC', 'betaC-alphaC-C', 'alphaC-betaC-Sv']
        guess_types = ['betaC', 'Sv']
        guess_dihedrals = [[0, 2, 5, 6]]
        guess_dihedral_names = ['N-alphaC-betaC-Sv']
        guess_dc = []
        guess_dc_names = []
        self.check_group(bonds, guess_bonds)
        self.check_single(bond_names, guess_bond_names)
        self.check_group(angles, guess_angles)
        self.check_single(angle_names, guess_angle_names)
        self.check_group(dihedrals, guess_dihedrals)
        self.check_group(dihedral_names, guess_dihedral_names)
        self.check_group(dc, guess_dc)
        self.check_group(dc_names, guess_dc_names)
        self.check_single(types, guess_types)
        self.assertEquals(pb.charge[-1], 0)

    def check_group(self, bbonds, gbonds):
        self.assertTrue(len(gbonds) == len(bbonds))
        for ind, bond in enumerate(bbonds):
            self.assertTrue(bond in gbonds)
        for ind, bond in enumerate(gbonds):
            self.assertTrue(bond in bbonds)

    def check_single(self, types, gtypes):
        for ind, type in enumerate(types):
            self.assertEquals(type, gtypes[ind])

if __name__ == '__main__':
    unittest.main()
