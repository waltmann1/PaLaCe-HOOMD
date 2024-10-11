from __future__ import division
import numpy as np
import hoomd
from hoomd import md
from Loggable import Loggable


class Bonds(Loggable):

    def __init__(self, log_list=None):

        super(Bonds, self).__init__(log_list)
        self.log_values = ['bond_harmonic_energy']
        self.names = []
        self.k = 170
        self.r0 = []
        self.bond_ref = None

        self.names.append('N-alphaC')
        self.r0.append(1.5)
        self.names.append('alphaC-C')
        self.r0.append(1.5)
        self.names.append('C-N')
        self.r0.append(1.3)
        self.names.append('alphaC-betaC')
        self.r0.append(1.5)

        self.names.append('betaC-Sc')
        self.r0.append(0.9)
        self.names.append('betaC-Sd')
        self.r0.append(1.4)
        self.names.append('betaC-Se')
        self.r0.append(0.8)
        self.names.append('betaC-Sl')
        self.r0.append(1.5)
        self.names.append('betaC-Sn')
        self.r0.append(1.4)
        self.names.append('betaC-Sp')
        self.r0.append(1.2)
        self.names.append('betaC-Sq')
        self.r0.append(0.8)
        self.names.append('betaC-Ss')
        self.r0.append(0.7)
        self.names.append('betaC-St')
        self.r0.append(0.6)
        self.names.append('betaC-Sv')
        self.r0.append(0.6)
        self.names.append('betaC-Sw')
        self.r0.append(0.8)

        self.names.append('Se-Te')
        self.r0.append(2.3)

        # the original version of the forcefield contains bonds for betaC-Tf (3.2) and and Sf-Tf(2.4)
        # despite the fact betaC and Sf should be in the same spot.
        # Instead only Sf atom is being used and the bond length will be 2.8

        self.names.append('Sf-Tf')
        self.r0.append(2.8)
        self.names.append('Sh-Th')
        self.r0.append(2.6)

        # the original version of the forcefield contains bonds for betaC-Th (3.0) and and Sh-Th(2.2)
        # despite the fact betaC and Sh should be in the same spot.
        # Instead only Sh atom is being used and the bond length will be 2.6
        self.names.append('Sq-Tq')
        self.r0.append(2.3)
        self.names.append('Sw-Tw')
        self.r0.append(2.9)

        # the original version of the forcefield contains bonds for betaC-Ty (3.2) and and Sy-Ty(2.4)
        # despite the fact betaC and Sy should be in the same spot.
        # Instead only Sy atom is being used and the bond length will be 2.8
        self.names.append('Sy-Ty')
        self.r0.append(2.8)

        self.names.append('betaC-gammaC1i')
        self.r0.append(1.5)
        self.names.append('betaC-gammaCk')
        self.r0.append(1.5)
        self.names.append('betaC-gammaCm')
        self.r0.append(1.5)
        self.names.append('betaC-gammaCr')
        self.r0.append(1.5)
        self.names.append('betaC-gammaC2i')
        self.r0.append(1.5)

        self.names.append('betaC-Tf')
        self.r0.append(3.2)
        self.names.append('betaC-Th')
        self.r0.append(3.0)
        self.names.append('betaC-Ty')
        self.r0.append(3.2)

        self.names.append('deltaSm-Tm')
        self.r0.append(0.9)
        self.names.append('epsilonCk-Tk')
        self.r0.append(0.6)
        self.names.append('epsilonNr-Tr')
        self.r0.append(0.9)
        self.names.append('gammaCk-deltaCk')
        self.r0.append(1.5)
        self.names.append('gammaCm-deltaSm')
        self.r0.append(1.8)
        self.names.append('gammaCr-deltaCr')
        self.r0.append(1.5)
        self.names.append('gammaC1i-deltaCi')
        self.r0.append(1.5)
        self.names.append('deltaCk-epsilonCk')
        self.r0.append(1.5)
        self.names.append('deltaCr-epsilonNr')
        self.r0.append(1.5)
        self.names.append('Sc-Sc')
        self.r0.append(2.8)
        self.names.append('halfway')
        self.r0.append(.75)
        self.names.append('Si-1')
        self.r0.append(1.3)
        self.names.append('exclusion')
        self.r0.append(0.0)
        self.names.append('C-O')
        self.r0.append(1.24)
        self.names.append('N-H')
        self.r0.append(1.00)
        self.names.append('O-H')
        self.r0.append(4.00)


    def set_all_harmonic_bonds(self, system):
        """

        :param system: the system that needs the parameters set
        :return: reference to the harmonic bond object
        """

        self.bond_ref = hoomd.md.bond.harmonic()
        self.add_to_logger()
        snap = system.take_snapshot(all=True)
        for b in snap.bonds.types:
        #for b in system.bonds:
            #name = str(b.type)
            name = str(b)
            #print(name)
            if name == 'halfway':
                self.bond_ref.bond_coeff.set(name, k=self.k, r0=self.r0[self.names.index(name)])
            elif name[:2] == 'Si':
                self.bond_ref.bond_coeff.set(name, k=self.k, r0=self.r0[self.names.index(name)])
            elif name == 'exclusion':
                self.bond_ref.bond_coeff.set(name, k=0, r0=0)
            else:
                self.bond_ref.bond_coeff.set(name, k=self.k, r0=self.r0[self.names.index(name)])
        del snap
        return self.bond_ref


