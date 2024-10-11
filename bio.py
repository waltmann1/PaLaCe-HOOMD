import Bio.PDB as p
from ProteinChain import ProteinChain

"""
parser = p.PDBParser()


struct = parser.get_structure('1sva','1sva.pdb')
#struct = parser.get_structure('1qo1','1qo1.pdb')


for model in struct:
        for ind,chain in enumerate(model):
            for ind2,residue in enumerate(chain):
                if ind == 0 and residue.resname == "LEU":
                    print(residue.resname)
                    print(" ")
                    atoms = []
                    for atom in residue:
                        atoms.append(atom)
                    print(atoms)
                    print(" ")
"""

parser = p.PDBParser()


struct = parser.get_structure('1sva','1sva.pdb')
#struct = parser.get_structure('1qo1','1qo1.pdb')


pro = 0
for model in struct:
        for ind,chain in enumerate(model):
            if ind == 0:
                pro = ProteinChain(chain,denature=True)
                pro.dump_xyz('test.xyz')