from __future__ import division
from ProteinBackbone import ProteinBackbone


class ProteinChain(ProteinBackbone):

    def __init__(self, chain, chain_index=0, denature=False, list_470=[]):
        """

        :param chain: biopython chain object
        """
        gen = chain.get_residues()
        residues = [r for r in gen]
        remove = []
        res_remove = []
        for i in range(len(residues)):
            try:
                __import__(residues[i].resname, fromlist=[''])
            except ImportError:
                if residues[i].resname not in res_remove:
                    print('Found illegal residue ' + str(residues[i].resname) + " in chain " + str(chain_index))
                    res_remove.append(residues[i].resname)
                remove.append(i)

        remove.sort()
        for ind in reversed(remove):
            del residues[ind]

        super(ProteinChain, self).__init__(len(residues), chain_index=chain_index, list_470=list_470)

        for i in range(len(residues)):
            data = [a for a in residues[i]]
            if denature:
                data = None
            self.add_residue(residues[i].resname, self.resnum_to_iac(i), data=data)
