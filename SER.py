from __future__ import division
import numpy as np
from ResidueAbs import ResidueAbs


class SER(ResidueAbs):

    def __init__(self, data=None, error_470=False):

        self.name = 'Serine'
        self.letter = 's'
        self.code = 'SER'
        self.bio_atoms = 6
        self.beads = 7
        super(SER, self).__init__(data, error_470=error_470)
        self.name = 'Serine'
        self.letter = 's'
        self.code = 'SER'
        self.bio_atoms = 6
        self.beads = 7

    def add_to_backbone(self, index_alpha_carbon, backbone):
        """

        :param index_alpha_carbon: index of the alpha carbon to add to from the backbone
        :param backbone: PaLaCeProtein backbone
        :return: void
        """
        if self.pos is not None:
            self.add_backbone_positions(index_alpha_carbon, backbone)

        down = self.add_beta_carbon(index_alpha_carbon, backbone)

        di = 1
        if down:
            di = -1

        backbone.type.append('Ss')
        backbone.mass.append(17)
        backbone.charge.append(0.0)
        backbone.body.append(-1)
        backbone.bonds.append([backbone.num_particles - 1, backbone.num_particles])
        backbone.bond_names.append('betaC-Ss')
        backbone.angles.append([index_alpha_carbon, backbone.num_particles - 1, backbone.num_particles])
        backbone.angle_names.append('alphaC-betaC-Ss')
        backbone.dihedrals.append([index_alpha_carbon - 2, index_alpha_carbon, backbone.num_particles - 1,
                                backbone.num_particles])
        backbone.dihedral_names.append('N-alphaC-betaC-Ss')
        if self.pos is not None:
            backbone.position.append(self.pos[5])
        else:
            backbone.position.append(np.add(backbone.position[-1], [0.7 * np.cos(np.deg2rad(20.5)),
                                                                    di * 0.7 * np.sin(np.deg2rad(20.5)), 0.0]))
        backbone.increment_num_particles(index_alpha_carbon)

    def read_data(self, data):
        """

        :param data: PDB data which is converted to a position array reflecting the model for a particular amino acid
        :return: the position array
        """

        code = self.check_data(data)
        if code is None:
            return None
        if code == 1:
            return data
        return self.convert_data(data)
