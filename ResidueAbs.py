import numpy as np


class ResidueAbs(object):

    def __init__(self, data=None, error_470=False):

        self.error_470 = error_470
        self.pos = self.read_data(self.fix_general_data(data))
        self.name = 'residue'
        self.letter = '.'
        self.dc_name = 'NOR'
        # self.bio_atoms = 0
        # self.beads = 0
        self.h_included = False
        code = self.check_data(self.fix_general_data(data))
        if code == 1:
            self.h_included = True

    def add_to_backbone(self, index_alpha_carbon, backbone):

        print('adding ' + self.name + 'to ' + backbone.name + ' at position ' + str(index_alpha_carbon))

    def add_backbone_positions(self, index_alpha_carbon, backbone):

        if self.h_included:
            for i in range(5):
                backbone.position[index_alpha_carbon - 2 + i] = self.pos[i]
            del self.pos[4]
        else:
            backbone.position[index_alpha_carbon - 2] = self.pos[0]
            backbone.position[index_alpha_carbon] = self.pos[1]
            backbone.position[index_alpha_carbon + 1] = self.pos[2]
            backbone.position[index_alpha_carbon + 2] = self.pos[3]
            h_pos = np.add(np.multiply(np.subtract(self.pos[3], self.pos[2]), 1.0/1.24), self.pos[0])
            backbone.position[index_alpha_carbon - 1] = h_pos

        # backbone.position[index_alpha_carbon + 1] = np.add(np.multiply(.425871, self.pos[2]),
                            #                               np.multiply(.571429, self.pos[3]))
        # backbone.position[index_alpha_carbon - 3] = np.add(np.multiply(.93333, self.pos[0]),
                             #                              np.multiply(.06666, h_pos))


    def add_beta_carbon(self, index_alpha_carbon, backbone):

        # particle types
        backbone.type.append('betaC')
        # particle positions
        if self.pos is None:
            pos = [0, 1.5, 0]
            if backbone.position[index_alpha_carbon - 1][1] - backbone.position[index_alpha_carbon - 2][1] > 0:
                    pos = [0, -1.5, 0]
            backbone.position = list(backbone.position)
            backbone.position.append(np.add(backbone.position[index_alpha_carbon], pos))
        else:
            backbone.position = list(backbone.position)
            backbone.position.append(self.pos[4])

        backbone.mass.append(14)
        backbone.charge.append(0)
        backbone.body.append(-1)
        backbone.bonds.append([index_alpha_carbon, backbone.num_particles])
        backbone.bond_names.append('alphaC-betaC')
        backbone.angles.append([index_alpha_carbon - 2, index_alpha_carbon, backbone.num_particles])
        backbone.angle_names.append('N-alphaC-betaC')
        backbone.angles.append([index_alpha_carbon + 1, index_alpha_carbon, backbone.num_particles])
        backbone.angle_names.append('betaC-alphaC-C')

        backbone.increment_num_particles(index_alpha_carbon)

        return backbone.position[index_alpha_carbon - 1][1] - backbone.position[index_alpha_carbon - 2][1] > 0

    def read_data(self, data):
        """

        :param data: PDB data which is converted to a position array reflecting the model for a particular amino acid
        :return: the position array
        """

        print('reading the relevant data')
        # TODO
        return self.convert_data(data)

    def convert_data(self, data):
        """

        :param data: list of biopython atoms
        :return:
        """
        if type(data[0]) is not list and not isinstance(data[0], np.ndarray):
            return [[at.coord[0], at.coord[1], at.coord[2]] for at in data]
        else:
            return data

    def get_average_position(self, data):
        """

        :param data: array of biopython atoms or just positions
        :return:
        """

        if type(data[0]) is not list and not isinstance(data[0], np.ndarray):
            return [np.mean([at.coord[0] for at in data]), np.mean([at.coord[1] for at in data]),
                    np.mean([at.coord[2] for at in data])]
        else:
            return [np.mean([at[0] for at in data]), np.mean([at[1] for at in data]),
                    np.mean([at[2] for at in data])]


    def fix_general_data(self, data):
        """

        :param data: biopython atom data
        :return: the good data
        """
        if data is None:
            return data
        if type(data[0]) is list or isinstance(data[0], np.ndarray):
            return data
        if data[-1].name == 'OXT':
            del data[-1]
        to_remove = []
        for ind, at in enumerate(data):
            if 'H' == at.name[0] or at.name == 'OXT':
                to_remove.append(ind)
        to_remove.sort()
        for nd in reversed(to_remove):
            del data[nd]
        if self.error_470:
            for i in range(10):
                data.append(data[-1])
        return data

    def check_data(self, data):

        if data is None:
            return None
        if self.error_470 and len(data) >= 4:
            return 2
        if self.error_470 and len(data) < 4:
            raise ValueError("even with remark 470 msut have 4 atoms to fill in backbone", len(data), data)

        if type(data[0]) is not list and not isinstance(data[0], np.ndarray):
                if len(data) != self.bio_atoms:
                    raise ValueError(self.code + " residue contains " + str(self.bio_atoms) + " Biopython atoms", len(data), data)
                return 2
        else:
                if len(data) != self.beads:
                    raise ValueError(self.code + " residue contains " + str(self.beads) + " Palace atoms", len(data), data)
                return 1