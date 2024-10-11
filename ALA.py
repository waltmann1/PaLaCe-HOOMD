from __future__ import division
from ResidueAbs import ResidueAbs


class ALA(ResidueAbs):

    def __init__(self, data=None, error_470=False):

        self.name = 'Alanine'
        self.letter = 'a'
        self.code = 'ALA'
        self.beads = 6
        self.bio_atoms = 5
        super(ALA, self).__init__(data, error_470=error_470)
        self.name = 'Alanine'
        self.letter = 'a'
        self.code = 'ALA'
        self.beads = 6
        self.bio_atoms = 5

    def add_to_backbone(self, index_alpha_carbon, backbone):
        """

        :param index_alpha_carbon: index of the alpha carbon to add to from the backbone
        :param backbone: PaLaCeProtein backbone
        :return: void
        """
        if self.pos is not None:
            self.add_backbone_positions(index_alpha_carbon, backbone)

        down = self.add_beta_carbon(index_alpha_carbon, backbone)

        backbone.type[-1] = 'Sa'
        # backbone.mass.append(12)
        # backbone.charge.append(0.0)
        # backbone.body[-1] = 666
        # backbone.body.append(666)

        # if self.pos is not None:
            # backbone.position.append(self.pos[4])
        # else:
            # backbone.position.append(backbone.position[-1])
        # backbone.increment_num_particles(index_alpha_carbon)

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
        return d
