from __future__ import division
import numpy as np
from ResidueAbs import ResidueAbs


class ILE(ResidueAbs):

    def __init__(self, data=None, error_470=False):

        self.name = 'Isoleucine'
        self.letter = 'i'
        self.code = 'ILE'
        self.dc_name = self.code
        self.beads = 10
        self.bio_atoms = 8
        super(ILE, self).__init__(data, error_470=error_470)
        self.name = 'Isoleucine'
        self.letter = 'i'
        self.code = 'ILE'
        self.dc_name = self.code
        self.beads = 10
        self.bio_atoms = 8

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

        backbone.type.append('gammaC1i')
        backbone.mass.append(12)
        backbone.charge.append(0.0)
        backbone.body.append(-1)
        backbone.bonds.append([backbone.num_particles - 1, backbone.num_particles])
        backbone.bond_names.append('betaC-gammaC1i')
        backbone.angles.append([index_alpha_carbon, backbone.num_particles - 1, backbone.num_particles])
        backbone.angle_names.append('alphaC-betaC-gammaC1i')

        if self.pos is not None:
            backbone.position.append(self.pos[5])
        else:
            backbone.position.append(np.add(backbone.position[-1], [1.5 * np.cos(np.deg2rad(20.5)),
                                                                    di * 1.5 * np.sin(np.deg2rad(20.5)), 0.0]))
        backbone.increment_num_particles(index_alpha_carbon)

        backbone.type.append('gammaC2i')
        backbone.mass.append(12)
        backbone.charge.append(0.0)
        backbone.body.append(-1)
        backbone.bonds.append([backbone.num_particles - 2, backbone.num_particles])
        backbone.bond_names.append('betaC-gammaC2i')
        backbone.angles.append([index_alpha_carbon, backbone.num_particles - 2, backbone.num_particles])
        backbone.angle_names.append('alphaC-betaC-gammaC2i')
        backbone.angles.append([backbone.num_particles - 1, backbone.num_particles - 2, backbone.num_particles])
        backbone.angle_names.append('gammaC1i-betaC-gammaC2i')

        if self.pos is not None:
            backbone.position.append(self.pos[6])
        else:
            backbone.position.append(np.add(backbone.position[-2], [-1.5 * np.cos(np.deg2rad(20.8)),
                                                                    di * 1.5 * np.sin(np.deg2rad(20.8)), 0.0]))
        backbone.increment_num_particles(index_alpha_carbon)

        backbone.type.append('deltaCi')
        backbone.mass.append(12)
        backbone.charge.append(0.0)
        backbone.body.append(-1)
        backbone.bonds.append([backbone.num_particles - 2, backbone.num_particles])
        backbone.bond_names.append('gammaC1i-deltaCi')
        backbone.angles.append([backbone.num_particles - 3, backbone.num_particles - 2, backbone.num_particles])
        backbone.angle_names.append('betaC-gammaC1i-deltaCi')
        backbone.dihedrals_coupling.append([index_alpha_carbon - 2, index_alpha_carbon, backbone.num_particles - 3,
                                            backbone.num_particles-2, index_alpha_carbon, backbone.num_particles - 3,
                                            backbone.num_particles-2, backbone.num_particles])
        backbone.dihedral_coupling_names.append('ILE')

        if self.pos is not None:
            backbone.position.append(self.pos[7])
        else:
            backbone.position.append(np.add(backbone.position[-2], [-1.5 * np.cos(np.deg2rad(93.8)),
                                                                    di * 1.5 * np.sin(np.deg2rad(93.8)), 0.0]))
        backbone.increment_num_particles(index_alpha_carbon)

        backbone.type.append('Si')
        backbone.mass.append(48)
        backbone.charge.append(0.0)
        backbone.body.append(-1)

        if self.pos is not None:
            backbone.position.append(self.pos[8])
        else:
            backbone.position.append(np.average(backbone.position[-4:], axis=0))

        backbone.bonds.append([backbone.num_particles - 1, backbone.num_particles])
        backbone.bond_names.append('Si-1')
        backbone.bonds.append([backbone.num_particles, index_alpha_carbon])
        backbone.bond_names.append("exclusion")
        backbone.bonds.append([backbone.num_particles - 4, backbone.num_particles])
        backbone.bond_names.append('Si-1')
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
        d = self.convert_data(data)
        d.append(self.get_average_position(data[4:]))
        return d
