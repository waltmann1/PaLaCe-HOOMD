from __future__ import division
import numpy as np
from ResidueAbs import ResidueAbs


class ARG(ResidueAbs):

    def __init__(self, data=None, error_470=False):

        self.name = 'Arganine'
        self.letter = 'r'
        self.code = 'ARG'
        self.dc_name = self.code
        self.bio_atoms = 11
        self.beads = 11
        super(ARG, self).__init__(data, error_470=error_470)
        self.name = 'Arganine'
        self.letter = 'r'
        self.code = 'ARG'
        self.dc_name = self.code
        self.bio_atoms = 11
        self.beads = 11

    def add_to_backbone(self, index_alpha_carbon, backbone):
        """

        :param index_alpha_carbon: index of the alpha carbon to add to from the backbone
        :param backbone: PaLaCeProtein backbone
        :return: void
        """
        if self.pos is not None:
            self.add_backbone_positions(index_alpha_carbon, backbone)

        down = self.add_beta_carbon(index_alpha_carbon, backbone)
        backbone.body[-1] = -1
        di = 1
        if down:
            di = -1

        backbone.type.append('Sr')
        backbone.mass.append(28)
        backbone.charge.append(0.0)
        backbone.body.append(-1)
        if self.pos is not None:
            backbone.position.append(self.pos[5])
        else:
            backbone.position.append(np.add(backbone.position[-1], [.75 * np.cos(np.deg2rad(24.3)),
                                                                    .75 * di * np.sin(np.deg2rad(24.3)), 0.0]))
        backbone.bonds.append([backbone.num_particles, backbone.num_particles - 1])
        backbone.bond_names.append("halfway")
        backbone.bonds.append([backbone.num_particles, index_alpha_carbon])
        backbone.bond_names.append("exclusion")
        backbone.increment_num_particles(index_alpha_carbon)

        backbone.type.append('gammaCr')
        backbone.mass.append(12)
        backbone.charge.append(0.0)
        backbone.body.append(-1)
        if self.pos is not None:
            backbone.position.append(self.pos[6])
        else:
            backbone.position.append(np.add(backbone.position[-2], [1.5 * np.cos(np.deg2rad(24.3)),
                                                                    1.5 * di * np.sin(np.deg2rad(24.3)), 0.0]))

        backbone.bonds.append([backbone.num_particles, backbone.num_particles - 1])
        backbone.bond_names.append("halfway")
        backbone.bonds.append([backbone.num_particles, backbone.num_particles - 2])
        backbone.bond_names.append("betaC-gammaCr")
        backbone.angles.append([index_alpha_carbon, backbone.num_particles - 2, backbone.num_particles])
        backbone.angle_names.append('alphaC-betaC-gammaCr')
        backbone.increment_num_particles(index_alpha_carbon)

        backbone.type.append('deltaCr')
        backbone.mass.append(12)
        backbone.charge.append(0.0)
        backbone.body.append(-1)
        if self.pos is not None:
            backbone.position.append(self.pos[7])
        else:
            backbone.position.append(np.add(backbone.position[-1], [-1.5 * np.cos(np.deg2rad(87.7)),
                                                                    1.5 * di * np.sin(np.deg2rad(87.7)), 0.0]))

        backbone.bonds.append([backbone.num_particles, backbone.num_particles - 1])
        backbone.bond_names.append("gammaCr-deltaCr")
        backbone.angles.append([backbone.num_particles - 3, backbone.num_particles - 1, backbone.num_particles])
        backbone.angle_names.append('betaC-gammaCr-deltaCr')
        backbone.dihedrals_coupling.append([index_alpha_carbon - 2, index_alpha_carbon, backbone.num_particles - 3,
                                            backbone.num_particles - 1, index_alpha_carbon, backbone.num_particles - 3,
                                            backbone.num_particles - 1, backbone.num_particles])
        backbone.dihedral_coupling_names.append('ARG')
        backbone.increment_num_particles(index_alpha_carbon)

        backbone.type.append('epsilonNr')
        backbone.mass.append(15)
        backbone.charge.append(0.0)
        backbone.body.append(-1)
        if self.pos is not None:
            backbone.position.append(self.pos[8])
        else:
            backbone.position.append(np.add(backbone.position[-1], [1.5 * np.cos(np.deg2rad(23.6)),
                                                                    1.5 * di * np.sin(np.deg2rad(23.6)), 0.0]))

        backbone.bonds.append([backbone.num_particles, backbone.num_particles - 1])
        backbone.bond_names.append("deltaCr-epsilonNr")
        backbone.angles.append([backbone.num_particles - 2, backbone.num_particles - 1, backbone.num_particles])
        backbone.angle_names.append('gammaCr-deltaCr-epsilonNr')
        backbone.increment_num_particles(index_alpha_carbon)

        backbone.type.append('Tr')
        backbone.mass.append(44)
        backbone.charge.append(1.0)
        backbone.body.append(-1)
        if self.pos is not None:
            backbone.position.append(self.pos[9])
        else:
            backbone.position.append(np.add(backbone.position[-1], [-0.9 * np.cos(np.deg2rad(87.3)),
                                                                    0.9 * di * np.sin(np.deg2rad(87.3)), 0.0]))

        backbone.bonds.append([backbone.num_particles, backbone.num_particles - 1])
        backbone.bond_names.append("epsilonNr-Tr")
        backbone.angles.append([backbone.num_particles - 2, backbone.num_particles - 1, backbone.num_particles])
        backbone.angle_names.append('deltaCr-epsilonNr-Tr')
        spot = backbone.num_particles
        backbone.dihedrals_coupling.append([spot - 5, spot - 3, spot - 2, spot - 1, spot - 3, spot - 2, spot - 1, spot])
        backbone.dihedral_coupling_names.append('xARG')
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
        d.append(self.get_average_position(data[4:6]))
        d += self.convert_data(data[5:8])
        d.append(self.get_average_position(data[8:]))
        return d
