from __future__ import division
import numpy as np
from ResidueAbs import ResidueAbs


class PHE(ResidueAbs):

    def __init__(self, data=None, error_470=False):
        self.name = 'Phenylalanine'
        self.letter = 'f'
        self.code = 'PHE'
        self.bio_atoms = 11
        self.beads = 7
        super(PHE, self).__init__(data, error_470=error_470)
        self.name = 'Phenylalanine'
        self.letter = 'f'
        self.code = 'PHE'
        self.bio_atoms = 11
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
        """
        backbone.type.append('Sf')
        backbone.mass.append(12)
        backbone.charge.append(0.0)
        backbone.body.append(666)
        backbone.position.append(backbone.position[-1])
        backbone.increment_num_particles(index_alpha_carbon)
        """
        backbone.type[-1] = 'Sf' 

        backbone.type.append('Tf')
        backbone.mass.append(64)
        backbone.charge.append(0.0)
        backbone.body.append(-1)

        if self.pos is None:

            st_vec = [2.8 * np.cos(np.deg2rad(24.3)), 2.8 * np.sin(np.deg2rad(24.3)), 0]

            if down:
                st_vec[1] = -1 * st_vec[1]
            backbone.position.append(np.add(backbone.position[-1], st_vec))

        else:
            backbone.position.append(self.pos[5])

        backbone.bonds.append([backbone.num_particles-1, backbone.num_particles])
        backbone.bond_names.append('Sf-Tf')

        #backbone.bonds.append([backbone.num_particles-1, backbone.num_particles])
        #backbone.bond_names.append('Sf-Tf')

        backbone.angles.append([index_alpha_carbon, backbone.num_particles-1, backbone.num_particles])
        backbone.angle_names.append('alphaC-betaC-Tf')

        backbone.dihedrals.append([index_alpha_carbon - 2, index_alpha_carbon, backbone.num_particles-1,
                                   backbone.num_particles])
        backbone.dihedral_names.append('N-alphaC-betaC-Tf')
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
        da = self.convert_data(data[:5])
        # da.append(da[-1])
        da.append(self.get_average_position(data[5:]))
        return da
