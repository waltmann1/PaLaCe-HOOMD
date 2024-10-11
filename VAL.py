from __future__ import division
import numpy as np
from ResidueAbs import ResidueAbs


class VAL(ResidueAbs):

    def __init__(self, data=None, error_470=False):

        self.name = 'Valine'
        self.letter = 'v'
        self.code = 'VAL'
        self.bio_atoms = 7
        self.beads = 7
        super(VAL, self).__init__(data, error_470=error_470)
        self.name = 'Valine'
        self.letter = 'v'
        self.code = 'VAL'
        self.bio_atoms = 7
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

        backbone.type.append('Sv')
        backbone.mass.append(30)
        backbone.charge.append(0.0)
        backbone.body.append(-1)
        backbone.bonds.append([backbone.num_particles - 1, backbone.num_particles])
        backbone.bond_names.append('betaC-Sv')
        backbone.angles.append([index_alpha_carbon, backbone.num_particles - 1, backbone.num_particles])
        backbone.angle_names.append('alphaC-betaC-Sv')
        backbone.dihedrals.append([index_alpha_carbon - 2, index_alpha_carbon, backbone.num_particles - 1,
                                backbone.num_particles])
        backbone.dihedral_names.append('N-alphaC-betaC-Sv')
        if self.pos is not None:
            backbone.position.append(self.pos[5])
        else:
            backbone.position.append(np.add(backbone.position[-1], [0.6 * np.cos(np.deg2rad(39)),
                                                                    di * 0.6 * np.sin(np.deg2rad(39)), 0.0]))
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
        d = self.convert_data(data[:5])
        d.append(self.get_average_position(data[5:]))
        return d
