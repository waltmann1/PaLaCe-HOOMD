from ResidueAbs import ResidueAbs
import numpy as np


class PRO(ResidueAbs):

    def __init__(self, data=None, error_470=False):

        self.name = 'Proline'
        self.letter = 'p'
        self.code = 'PRO'
        self.dc_name = self.code
        self.bio_atoms = 7
        self.beads = 7
        super(PRO, self).__init__(data, error_470=error_470)
        self.name = 'Proline'
        self.letter = 'p'
        self.code = 'PRO'
        self.dc_name = self.code
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
        backbone.type[index_alpha_carbon] = 'alphaC' + self.letter
        res_number = backbone.iac_to_resnum(index_alpha_carbon)
        if res_number > 0 and res_number < backbone.residues - 1:
            backbone.dihedral_coupling_names[res_number - 1] = 'PRO'
        res_number -= 1
        if res_number > 0 and res_number < backbone.residues - 1:
            backbone.dihedral_coupling_names[res_number - 1] = 'PRP'
            backbone.dihedral_names[res_number] = 'alphaC-C-N-alphaCproline'

        down = self.add_beta_carbon(index_alpha_carbon, backbone)
        di = 1
        if down:
            di = -1

        backbone.type.append('Sp')
        backbone.mass.append(28)
        backbone.charge.append(0.0)
        backbone.body.append(-1)
        if self.pos is not None:
            backbone.position.append(self.pos[5])
        else:
            backbone.position.append(np.add(backbone.position[-1], [1.2 * np.sin(1), 1.2 * di * np.cos(1), 0.0]))
        backbone.bonds.append([backbone.num_particles, backbone.num_particles - 1])
        backbone.bond_names.append('betaC-Sp')
        backbone.bonds.append([backbone.num_particles, index_alpha_carbon])
        backbone.bond_names.append("exclusion")
        backbone.dihedrals.append([index_alpha_carbon - 2, index_alpha_carbon, backbone.num_particles - 1,
                                   backbone.num_particles])
        backbone.dihedral_names.append('N-alphaC-betaC-Sp')
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
