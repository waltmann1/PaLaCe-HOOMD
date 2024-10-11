import unittest
import shutil
import os
import numpy as np
from ProteinBackbone import ProteinBackbone


class TestProteinBackbone(unittest.TestCase):

    def setUp(self):
        self.dir = "dumpfolder/"
        directory = os.path.dirname(self.dir)
        if not os.path.exists(directory):
            os.makedirs(directory)

    def tearDown(self):
        shutil.rmtree(self.dir)

    def test_basic_backbone(self):

        pb = ProteinBackbone(3)
        bonds = [[0, 2], [2, 3], [3, 5], [5, 7], [7, 8], [8, 10], [10, 12], [12, 13]]
        for ind, bond in enumerate(pb.bonds):
            self.assertTrue(np.allclose(bond, bonds[ind]))
        angles = [[0, 2, 3], [2, 3, 5],  [3, 5, 7],  [5, 7, 8], [7, 8, 10], [8, 10, 12], [10, 12, 13]]
        for ind, angle in enumerate(pb.angles):
            self.assertTrue(np.allclose(angle, angles[ind]))
        dihedrals = [[2, 3, 5, 7], [7, 8, 10, 12]]
        for ind, dih in enumerate(pb.dihedrals):
            self.assertTrue(np.allclose(dih, dihedrals[ind]))
        dihedrals_coupling = [[3, 5, 7, 8, 5, 7, 8, 10]]
        for ind, dc in enumerate(pb.dihedrals_coupling):
            self.assertTrue(np.allclose(dc, dihedrals_coupling[ind]))

    def test_dc_names(self):
        pb = ProteinBackbone(5)
        pb.add_residue('GLY', 2)
        pb.add_residue('GLY', 7)
        pb.add_residue('PRO', 12)
        pb.add_residue('VAL', 17)
        pb.add_residue('GLY', 22)
        names = ['alphaC-C-N-alphaC', 'alphaC-C-N-alphaCproline', 'alphaC-C-N-alphaC', 'alphaC-C-N-alphaC']
        pb.dihedrals.remove(pb.dihedrals[-1])
        for ind, name in enumerate(pb.dihedral_names[:4]):
            self.dihedral_type_check(pb.dihedrals[ind], pb)
            self.assertEquals(names[ind], name)
        self.assertEquals(4, len(pb.dihedrals))
        names = ['PRP', 'PRO', 'NOR']
        for ind, name in enumerate(pb.dihedral_coupling_names):
            self.assertEquals(names[ind], name)
            self.dihedral_coupling_check(pb.dihedrals_coupling[ind], pb)
        self.assertEquals(3, len(pb.dihedrals_coupling))

    def test_dc_names_add_to_C(self):
        pb = ProteinBackbone(3)
        pb.add_residue('GLY', 2)
        pb.add_residue('GLY', 7)
        pb.add_residue('PRO', 12)
        pb2 = ProteinBackbone(2)
        pb2.add_residue('VAL', 2)
        pb2.add_residue('GLY', 7)
        pb.add_chain_to_C_terminal(pb2)
        names = ['alphaC-C-N-alphaC', 'alphaC-C-N-alphaCproline', 'alphaC-C-N-alphaC'
                 , 'alphaC-C-N-alphaC']
        pb.dihedral_names.remove(pb.dihedral_names[3])
        pb.dihedrals = pb.dihedrals[:3] + pb.dihedrals[4:]
        for ind, name in enumerate(pb.dihedral_names):
            self.assertEquals(names[ind], name)
            self.dihedral_type_check(pb.dihedrals[ind], pb)
        self.assertEquals(4, len(pb.dihedrals))
        names = ['PRP', 'PRO', 'NOR']
        for ind, name in enumerate(pb.dihedral_coupling_names):
            self.assertEquals(names[ind], name)
            self.dihedral_coupling_check(pb.dihedrals_coupling[ind], pb)
        self.assertEquals(3, len(pb.dihedrals_coupling))

        pb = ProteinBackbone(4)
        pb.add_residue('GLY', 2)
        pb.add_residue('GLY', 7)
        pb.add_residue('PRO', 12)
        pb.add_residue('VAL', 17)
        pb2 = ProteinBackbone(1)
        pb2.add_residue('GLY', 2)
        pb.add_chain_to_C_terminal(pb2)
        names = ['alphaC-C-N-alphaC', 'alphaC-C-N-alphaCproline', 'alphaC-C-N-alphaC'
                 , 'alphaC-C-N-alphaC']
        pb.dihedral_names.remove(pb.dihedral_names[3])
        pb.dihedrals = pb.dihedrals[:3] + pb.dihedrals[4:]
        for ind, name in enumerate(pb.dihedral_names):
            self.assertEquals(names[ind], name)
            self.dihedral_type_check(pb.dihedrals[ind], pb)
        self.assertEquals(4, len(pb.dihedrals))
        names = ['PRP', 'PRO', 'NOR']
        for ind, name in enumerate(pb.dihedral_coupling_names):
            self.assertEquals(names[ind], name)
            self.dihedral_coupling_check(pb.dihedrals_coupling[ind], pb)
        self.assertEquals(3, len(pb.dihedrals_coupling))

    def test_dc_names_add_to_N(self):
        pb2 = ProteinBackbone(3)
        pb2.add_residue('GLY', 2)
        pb2.add_residue('GLY', 7)
        pb2.add_residue('PRO', 12)
        pb = ProteinBackbone(2)
        pb.add_residue('VAL', 2)
        pb.add_residue('GLY', 7)
        pb.add_chain_to_N_terminal(pb2)
        names = ['alphaC-C-N-alphaC', 'alphaC-C-N-alphaCproline', 'alphaC-C-N-alphaC'
                 , 'alphaC-C-N-alphaC']

        pb.dihedral_names.remove(pb.dihedral_names[3])
        pb.dihedrals = pb.dihedrals[:3] + pb.dihedrals[4:]
        for ind, name in enumerate(pb.dihedral_names):
            self.assertEquals(names[ind], name)
            self.dihedral_type_check(pb.dihedrals[ind], pb)
        self.assertEquals(4, len(pb.dihedrals))
        names = ['PRP', 'PRO', 'NOR']
        for ind, name in enumerate(pb.dihedral_coupling_names):
            self.assertEquals(names[ind], name)
            self.dihedral_coupling_check(pb.dihedrals_coupling[ind], pb)
        self.assertEquals(3, len(pb.dihedrals_coupling))

        pb2 = ProteinBackbone(4)
        pb2.add_residue('GLY', 2)
        pb2.add_residue('GLY', 7)
        pb2.add_residue('PRO', 12)
        pb2.add_residue('VAL', 17)
        pb = ProteinBackbone(1)
        pb.add_residue('GLY', 2)
        pb.add_chain_to_N_terminal(pb2)
        names = ['alphaC-C-N-alphaC', 'alphaC-C-N-alphaCproline', 'alphaC-C-N-alphaC'
            , 'alphaC-C-N-alphaC']

        pb.dihedral_names.remove(pb.dihedral_names[3])
        pb.dihedrals = pb.dihedrals[:3] + pb.dihedrals[4:]
        for ind, name in enumerate(pb.dihedral_names):
            self.assertEquals(names[ind], name)
            self.dihedral_type_check(pb.dihedrals[ind], pb)
        self.assertEquals(4, len(pb.dihedrals))
        names = ['PRP', 'PRO', 'NOR']
        for ind, name in enumerate(pb.dihedral_coupling_names):
            self.assertEquals(names[ind], name)
            self.dihedral_coupling_check(pb.dihedrals_coupling[ind], pb)
        self.assertEquals(3, len(pb.dihedrals_coupling))

    def dihedral_type_check(self, dihedral, backbone):

        self.assertTrue(backbone.type[dihedral[0]][:6] == 'alphaC')
        self.assertTrue(backbone.type[dihedral[1]] == 'C')
        self.assertTrue(backbone.type[dihedral[2]] == 'N')
        self.assertTrue(backbone.type[dihedral[3]][:6] == 'alphaC')

    def dihedral_coupling_check(self, dihedral_coupling, backbone):

        self.assertTrue(backbone.type[dihedral_coupling[0]] == 'C')
        self.assertTrue(backbone.type[dihedral_coupling[1]] == 'N')
        self.assertTrue(backbone.type[dihedral_coupling[2]][:6] == 'alphaC')
        self.assertTrue(backbone.type[dihedral_coupling[3]] == 'C')
        self.assertTrue(backbone.type[dihedral_coupling[4]] == 'N')
        self.assertTrue(backbone.type[dihedral_coupling[5]][:6] == 'alphaC')
        self.assertTrue(backbone.type[dihedral_coupling[6]] == 'C')
        self.assertTrue(backbone.type[dihedral_coupling[7]] == 'N')

    def test_body(self):
        pb = ProteinBackbone(5)
        pb.add_residue('GLY', 2)
        pb.add_residue('GLY', 7)
        pb.add_residue('PRO', 12)
        pb.add_residue('VAL', 17)
        pb.add_residue('GLY', 22)
        pb.make_rigid(2, 0)
        for i in pb.res_indexes[2]:
            self.assertEquals(pb.body[i], 0)


if __name__ == '__main__':
    unittest.main()
