from __future__ import division
import numpy as np
import hoomd
from hoomd import md
from Loggable import Loggable


def palace_torsion(theta, w_i, phi_i):
    V = 0
    F = 0
    for i in range(0, 3):
        V += w_i[i] * np.cos(-1 * (i * theta + phi_i[i]))
        F += w_i[i] * i * np.sin(-1 * (i*theta + phi_i[i]))
    return V, F


class Dihedrals(Loggable):

    def __init__(self, log_list=None):
        super(Dihedrals, self).__init__(log_list)
        self.log_values = ['dihedral_table_energy']
        self.names = []
        self.w_i = []
        self.phi_i = []
        self.dref = None

        self.names.append('alphaC-C-N-alphaC')
        self.w_i.append([2.2, 9.0, 0])
        self.phi_i.append(np.deg2rad([2.0, 181.5, 0]))

        self.names.append('alphaC-C-N-alphaCproline')
        self.w_i.append([0.9, 10.5, 0.0])
        self.phi_i.append(np.deg2rad([28.3, 178.8, 0]))

        self.names.append('N-alphaC-betaC-Sc')
        self.w_i.append([0.7, 0.7, 2.6])
        self.phi_i.append(np.deg2rad([-40.7, 6.0, -4.1]))

        self.names.append('N-alphaC-betaC-Sd')
        self.w_i.append([1.0, 0.9, 2.1])
        self.phi_i.append(np.deg2rad([-62.9, 56.3, -0.4]))

        self.names.append('N-alphaC-betaC-Sl')
        self.w_i.append([2.1, 1.2, 2.0])
        self.phi_i.append(np.deg2rad([-48.2, 34.1, 0.0]))

        self.names.append('N-alphaC-betaC-Sn')
        self.w_i.append([-1.1, 0.8, 2.1])
        self.phi_i.append(np.deg2rad([116.6, 45.8, -4.4]))

        self.names.append('N-alphaC-betaC-Sp')
        self.w_i.append([-858.5, 222.5, 0.0])
        self.phi_i.append(np.deg2rad([-2.4, -4.8, 0.0]))

        self.names.append('N-alphaC-betaC-Ss')
        self.w_i.append([0.3, 0.4, 1.7])
        self.phi_i.append(np.deg2rad([30.5, 2.6, -1.7]))

        self.names.append('N-alphaC-betaC-St')
        self.w_i.append([-0.2, 0.4, 1.7])
        self.phi_i.append(np.deg2rad([30.5, 2.6, -1.7]))

        self.names.append('N-alphaC-betaC-Sv')
        self.w_i.append([0.3, 0.5, 1.7])
        self.phi_i.append(np.deg2rad([-4.6, 104.7, 178.9]))

        self.names.append('N-alphaC-betaC-Tf')
        self.w_i.append([1.1, 0.7, 2.6])
        self.phi_i.append(np.deg2rad([-48.0, 29.9, 3.3]))

        self.names.append('N-alphaC-betaC-Th')
        self.w_i.append([0.8, 0.5, 2.3])
        self.phi_i.append(np.deg2rad([-50.0, 19.4, -1.7]))

        self.names.append('N-alphaC-betaC-Ty')
        self.w_i.append([0.9, 0.6, 2.3])
        self.phi_i.append(np.deg2rad([-38.9, 10.3, 1.7]))

        self.names.append('betaC-gammaCm-deltaSm-Tm')
        self.w_i.append([0.3, 0.6, 0.9])
        self.phi_i.append(np.deg2rad([2.1, -6.6, -2.9]))

        self.names.append('betaC-Sc-Sc-betaC')
        self.w_i.append([4.3, 4.7, 0])
        self.phi_i.append(np.deg2rad([-3.8, -3.4, 0]))

    def set_all_palace_dihedrals(self, system):

        self.dref = hoomd.md.dihedral.table(width=1000)
        self.add_to_logger()
        snap = system.take_snapshot(all=True)
        for d in snap.dihedrals.types:
        #for d in system.dihedrals:
            name = str(d)
            self.dref.dihedral_coeff.set(name, func=palace_torsion,
                                                    coeff=dict(w_i=self.w_i[self.names.index(name)],
                                                               phi_i=self.phi_i[self.names.index(name)]))
        return self.dref