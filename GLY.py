from __future__ import division
from ResidueAbs import ResidueAbs


class GLY(ResidueAbs):

    def __init__(self, data=None, error_470=False):

        self.name = 'Glycine'
        self.letter = 'g'
        self.code = 'GLY'
        self.dc_name = self.code
        self.bio_atoms = 4
        self.beads = 5
        super(GLY, self).__init__(data, error_470=error_470)
        self.name = 'Glycine'
        self.letter = 'g'
        self.code = 'GLY'
        self.dc_name = self.code
        self.bio_atoms = 4
        self.beads = 5

    def add_to_backbone(self, index_alpha_carbon, backbone):
        """

        :param index_alpha_carbon: index of the alpha carbon to add to from the backbone
        :param backbone: PaLaCeProtein backbone
        :return: void
        """
        if self.pos is not None:
                self.add_backbone_positions(index_alpha_carbon, backbone)
        backbone.type[index_alpha_carbon] = 'alphaC' + self.letter
        res_number = backbone.iac_to_resnum(index_alpha_carbon)
        if res_number > 0 and res_number < backbone.residues - 1:
            backbone.dihedral_coupling_names[res_number - 1] = 'GLY'

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
