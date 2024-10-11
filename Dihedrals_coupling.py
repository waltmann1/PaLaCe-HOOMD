from __future__ import division
import numpy as np
import hoomd
from hoomd import md
from Loggable import Loggable

class Dihedrals_coupling(Loggable):

    def __init__(self, log_list=None):
        super(Dihedrals_coupling, self).__init__(log_list)
        self.log_values = ['dihedral_coupling_opls_energy']
        self.names = []
        self.w_ij = []
        self.phi_one_ij = []
        self.phi_two_ij = []
        self.dc_ref = None

        self.names.append('NOR')
        self.w_ij.append([0.0, 2.4, -1.1, 0.8, 5.5, 3.6, 3.5, 1.1, 2.4, 2.5, 5.1, -1.3, 1.5, 0.7, 2.2, -0.5])
        self.phi_one_ij.append(np.deg2rad([0, 0, 0, 0, -71.7, 68.8, -19.6, -244.2, 45.7, 1.5, -79.6, 22.6, -33.3, 81.9,
                                           48.8, 143.3]))
        self.phi_two_ij.append(np.deg2rad([0, 99.2, 60.0, -96.0, 0.0, -53.4, 70.2, 52.5, 0.0, -8.1, 91.3, 158.6, 0.0,
                                           163.9, -64.3, -9.2]))

        self.names.append('PRP')
        self.w_ij.append([0.0, 4.0, -0.4, 0.9, 2.8, -2.1, 3.4, 1.1, 2.0, 1.9, 2.6, 1.1, 0.3, 1.8, 1.5, -0.6])
        self.phi_one_ij.append(np.deg2rad([0.0, 0.0, 0.0, 0.0, -15.5, 78.2, -102.0, 1.1, 4.3, 131.2, -74.8, -118.1,
                                           -1.8, -35.2, -77.3, -59.8]))
        self.phi_two_ij.append(np.deg2rad([0.0, 50.5, -26.3, -70.0, 0.0, 117.5, 108.3, 251.8, 0.0, -78.1, 82.6, 65.4,
                                           0.0, 136.3, 114.1, 46.8]))

        self.names.append('PRO')
        self.w_ij.append([0.0, 81.1, 101.5, -70.1, -500.7, -118.1, -150.9, -108.0, 181.3, -67.2, 68.6, 61.2, 27.4, 23.9,
                          16.9, 19.7])
        self.phi_one_ij.append(np.deg2rad([0.0, 0.0, 0.0, 0.0, 64.1, 251.8, 216.9, -147.7, 126.6, -39.7, 70.2, 60.6,
                                           6.7, -144.0, -79.1, 106.3]))
        self.phi_two_ij.append(np.deg2rad([0.0, -17.7, 211.2, 97.4, 0.0, 865.2, 390.5, 99.0, 0.0, 286.1, 211.9, -76.0,
                                           0.0, 87.8, 204.1, 111.1]))

        self.names.append('GLY')
        self.w_ij.append([0.0, 1.6, 1.9, -0.5, 2.7, -0.7, 1.1, 0.5, 2.8, 1.7, 2.1, 0.5, 1.0, 0.5, -1.0, 0.5])
        self.phi_one_ij.append(np.deg2rad([0.0, 0.0, 0.0, 0.0, 8.7, 57.9, 97.5, 40.1, -6.8, 15.2, -92.8, -34.1, 13.7,
                                           73.2, 101.5, 220.7]))
        self.phi_two_ij.append(np.deg2rad([0.0, -46.8, 180.7, -15.9, 0.0, 153.7, -96.3, 201.6, 0.0, -62.7, 99.0, 88.4,
                                           0.0, -24.4, 81.8, -3.4]))

        self.names.append('ARG')
        self.w_ij.append([0.0, 5.1, 1.5, 0.9, 4.1, 4.6, 2.5, 0.8, -1.1, -2.0, 1.7, 1.2, 2.0, -0.8, 1.2, -0.6])
        self.phi_one_ij.append(np.deg2rad([0.0, 0.0, 0.0, 0.0, -75.1, 84.2, 101.2, 279.6, 36.8, 49.9, -5.7, 11.4, 12.4,
                                           -108.2, -59.7, 82.6]))
        self.phi_two_ij.append(np.deg2rad([0.0, 25.8, 84.0, 26.1, 0.0, -123.0, -48.7, 187.6, 0.0, 30.1, 42.7, 61.8,
                                           0.0, 175.3, 20.5, 49.8]))

        self.names.append('GLN')
        self.w_ij.append([0.0, 3.0, 1.4, 1.7, 2.5, 2.4, 1.7, 1.0, 2.5, 2.4, 1.7, 1.0, 1.5, 1.0, 0.7, 0.6, 2.2, -0.3,
                          0.7, 0.9])
        self.phi_one_ij.append(np.deg2rad([0.0, 0.0, 0.0, 0.0, -30.6, -23.8, 2.7, 17.2, -20.7, -16.1, 57.0, 61.1, 13.5,
                                           28.9, -163.9, 13.2]))
        self.phi_two_ij.append(np.deg2rad([0.0, -8.5, -30.0, 12.7, 0.0, -52.2, -78.0, -80.6, 0.0, 11.7, -35.3, -111.8,
                                           0.0, -1.1, 100.4, 194.3]))

        self.names.append('GLU')
        self.w_ij.append([0.0, 3.1, 1.3, 2.3, 2.1, 1.0, 1.9, 1.7, 0.6, 1.6, 1.0, 0.8, 2.3, -0.6, -0.5, 0.8])
        self.phi_one_ij.append(np.deg2rad([0.0, 0.0, 0.0, 0.0, -64.0, -41.7, 63.4, 56.6, 32.1, 15.6, 15.5, 17.5, 16.5,
                                           116.4, 24.8, -25.5]))
        self.phi_two_ij.append(np.deg2rad([0.0, 0.2, 32.9, 27.7, 0.0, -44.5, -61.5, -90.0, 0.0, 48.5, -1.7, -112.8, 0.0,
                                           61.9, 85.6, 220.2]))

        self.names.append('ILE')
        self.w_ij.append([0.0, 3.1, -1.7, 3.3, 2.3, 1.4, -1.5, 1.5, 1.2, 1.3, 1.0, 1.3, 3.5, -0.6, 0.8, 0.2])
        self.phi_one_ij.append(np.deg2rad([0.0, 0.0, 0.0, 0.0, -95.7, 39.6, 3.3, 179.4, 44.2, -28.5, -2.5, -10.7, 8.0,
                                           87.1, 32.6, -39.5]))
        self.phi_two_ij.append(np.deg2rad([0.0, 22.1, 172.2, 26.6, 0.0, 188.0, 80.7, 128.8, 0.0, 88.5, 13.7, 9.6, 0.0,
                                           154.5, 55.5, 10.7]))

        self.names.append('LYS')
        self.w_ij.append([0.0, 3.2, -0.7, 1.4, 1.7, 2.6, 2.6, 1.0, 1.1, -0.9, 1.2, 0.4, 2.5, 1.0, 0.8, 0.3])
        self.phi_one_ij.append(np.deg2rad([0.0, 0.0, 0.0, 0.0, -40.3, -4.4, -12.1, 176.7, 11.1, -62.2, -38.3, -191.9,
                                           20.2, -102.3, 41.4, 154.6]))
        self.phi_two_ij.append(np.deg2rad([0.0, -13.8, 127.0, 9.9, 0.0, -112.0, -138.4, 52.0, 0.0, -71.2, 146.5, -7.2,
                                           0.0, 127.5, -26.9, 82.4]))

        self.names.append('MET')
        self.w_ij.append([0.0, 3.3, 1.2, 1.8, 1.6, 2.2, 1.8, 1.5, 0.7, 1.2, 2.7, -1.0, 2.1, 1.2, 0.7, 0.8])
        self.phi_one_ij.append(np.deg2rad([0.0, 0.0, 0.0, 0.0, -90.5, 15.7, 136.6, 143.0, 254.2, 70.9, 19.6, 19.2, 7.9,
                                           97.0, 152.1, 191.3]))
        self.phi_two_ij.append(np.deg2rad([0.0, -8.1, -33.6, -8.4, 0.0, -136.5, 47.2, 66.5, 0.0, 133.3, 70.0, -92.9,
                                           0.0, 133.4, -178.7, 46.8]))

        self.names.append('TRP')
        self.w_ij.append([0.0, 3.0, 2.9, 1.0, 0.9, 5.3, 3.8, 1.5, -0.4, 2.5, 2.0, 1.9, 2.8, -1.7, -0.9, -1.1])
        self.phi_one_ij.append(np.deg2rad([0.0, 0.0, 0.0, 0.0, 11.2, 61.0, 54.1, 176.8, -23.8, -16.0, -51.5, 98.5, 11.6,
                                           74.7, 144.2, 156.9]))
        self.phi_two_ij.append(np.deg2rad([0.0, -34.9, 1.7, -8.8, 0.0, -81.6, -114.4, 96.1, 0.0, 175.8, 32.1, 31.5, 0.0,
                                           -26.0, 44.6, 0.7]))

        self.names.append('xARG')
        self.w_ij.append([0.0, 5.5, 3.3, 1.5, 1.6, 2.8, 2.2, 1.6, 0.4, 0.5, 0.7, 0.6, 2.1, -0.3, -0.5, -0.2])
        self.phi_one_ij.append(np.deg2rad([0.0, 0.0, 0.0, 0.0, 49.9, 61.7, 58.7, 76.7, 174.9, 14.8, -19.8, 12.4, 12.0,
                                           126.0, 44.9, 111.0]))
        self.phi_two_ij.append(np.deg2rad([0.0, 7.0, 9.5, 17.6, 0.0, 27.4, 46.5, 91.2, 0.0, 121.1, 94.6, 67.7, 0.0,
                                           42.7, 45.8, -18.0]))

        self.names.append('xLYS')
        self.w_ij.append([0.0, 2.4, 0.5, 1.3, 3.0, 1.5, 0.6, 1.8, 0.9, 1.4, 1.1, 1.0, 1.6, 1.1, 1.1, 0.5])
        self.phi_one_ij.append(np.deg2rad([0.0, 0.0, 0.0, 0.0, 0.4, -10.6, 45.1, 83.6, 43.6, 114.7, 130.7, 71.1, 29.5,
                                           115.3, 104.0, 118.0]))
        self.phi_two_ij.append(np.deg2rad([0.0, 11.5, 17.7, 13.2, 0.0, 24.6, 160.6, 108.2, 0.0, -53.9, -4.4, 64.2, 0.0,
                                           -42.3, -25.6, 56.8]))

    def set_all_dihedrals_coupling(self, system):


        snap = system.take_snapshot(all=True)
        types = snap.dihedrals_coupling.types
        self.dc_ref = hoomd.md.dihedral_coupling.opls_phase()
        self.add_to_logger()
        for dc in types:
        #for dc in system.dihedrals_coupling:
            name = str(dc)
            n = self.names.index(name)
            self.dc_ref.dihedral_coupling_coeff.set(name, Wij=self.w_ij[n], phi1ij=self.phi_one_ij[n],
                                                   phi2ij=self.phi_two_ij[n])
        del snap
        return self.dc_ref

    def reset_all_dihedrals_coupling(self, system, nrg_factor):

        if self.dc_ref is None:
            print("Warning: must call set_all_dihedrals_coupling first. Ignoring reset_all_dihedrals_coupling.")
            return None

        set = []
        for dc in system.dihedrals_coupling:
            if dc.type not in set:
                n = self.names.index(dc.type)
                self.dc_ref.dihedral_coupling_coeff.set(str(dc.type), Wij=list(np.multiply(self.w_ij[n], nrg_factor)),
                                                        phi1ij=self.phi_one_ij[n], phi2ij=self.phi_two_ij[n])
                set.append(dc.type)



    def write_potential(self, name):

        results = []
        for phi in range(0, 180, 1):
            phi_rad = np.deg2rad(phi)
            for psi in range(0, 180, 1):
                psi_rad = np.deg2rad(psi)
                results.append([phi_rad, psi_rad, self.coupled_form(name, phi, psi)])
        f = open(name + str("_ramachandran.txt"))

        f.write('phi, psi, energy')
        for res in results:
            f.write("%f %f %f", res[0], res[1], res[2])
        f.close()

    def coupled_form(self, name, t1, t2):

        if name not in self.names:
            return 0
        y = self.names.index(name)
        total = 0
        for i in range(4):
            for j in range(4):
                a_idx = 4 * i + j
                total += self.w_ij[y][a_idx] * np.cos((i+1) * t1 + self.phi_one_ij[y][a_idx]) * \
                         np.cos((j+1) * t2 + self.phi_two_ij[y][a_idx])

