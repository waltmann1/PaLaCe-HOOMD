from __future__ import division
from ProteinAbs import ProteinAbs


class RigidRod(ProteinAbs):

    def __init__(self, length, charge_density):
        """

        :param chain: biopython chain object
        """

        super(RigidRod, self).__init__(length)

        self.rigid_count = 1
        self.num_particles = length

        for ind in range(self.residues):
            self.position.append([10 * ind, 0, 0])
            self.body.append(0)
            self.charge.append(charge_density * 10)
            self.type.append('bead')
            self.mass.append(116)
            self.body.append(0)
