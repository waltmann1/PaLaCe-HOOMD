from __future__ import division
import numpy as np
import hoomd
from hoomd import md
from Loggable import Loggable


class Angles(Loggable):

    def __init__(self, log_list=None):

        super(Angles, self).__init__(log_list)
        self.log_values = ['angle_harmonic_energy']
        self.names = []
        self.k = 90
        self.theta = []
        self.angle_ref = None

        self.names.append('N-alphaC-C')
        self.theta.append(np.deg2rad(110.9))
        self.names.append('alphaC-C-N')
        self.theta.append(np.deg2rad(117.3))
        self.names.append('C-N-alphaC')
        self.theta.append(np.deg2rad(122))
        self.names.append('N-alphaC-betaC')
        self.theta.append(np.deg2rad(110.3))
        self.names.append('betaC-alphaC-C')
        self.theta.append(np.deg2rad(109.8))

        self.names.append('alphaC-betaC-Sc')
        self.theta.append(np.deg2rad(113.8))
        self.names.append('alphaC-betaC-Sd')
        self.theta.append(np.deg2rad(114.0))
        self.names.append('alphaC-betaC-Se')
        self.theta.append(np.deg2rad(114.3))
        self.names.append('alphaC-betaC-Sl')
        self.theta.append(np.deg2rad(124.0))
        self.names.append('alphaC-betaC-Sn')
        self.theta.append(np.deg2rad(113.4))
        self.names.append('alphaC-betaC-Sq')
        self.theta.append(np.deg2rad(114.3))
        self.names.append('alphaC-betaC-Ss')
        self.theta.append(np.deg2rad(110.5))
        self.names.append('alphaC-betaC-St')
        self.theta.append(np.deg2rad(126.8))
        self.names.append('alphaC-betaC-Sv')
        self.theta.append(np.deg2rad(129.0))
        self.names.append('alphaC-betaC-Sw')
        self.theta.append(np.deg2rad(114.3))
        self.names.append('alphaC-betaC-Tf')
        self.theta.append(np.deg2rad(114.3))
        self.names.append('alphaC-betaC-Th')
        self.theta.append(np.deg2rad(114.3))
        self.names.append('alphaC-betaC-Ty')
        self.theta.append(np.deg2rad(114.3))
        self.names.append('alphaC-betaC-gammaC1i')
        self.theta.append(np.deg2rad(110.8))
        self.names.append('alphaC-betaC-gammaCk')
        self.theta.append(np.deg2rad(114.3))
        self.names.append('alphaC-betaC-gammaCm')
        self.theta.append(np.deg2rad(114.3))
        self.names.append('alphaC-betaC-gammaCr')
        self.theta.append(np.deg2rad(114.3))
        self.names.append('alphaC-betaC-gammaC2i')
        self.theta.append(np.deg2rad(110.8))
        self.names.append('betaC-gammaCk-deltaCk')
        self.theta.append(np.deg2rad(112.0))
        self.names.append('betaC-gammaCm-deltaSm')
        self.theta.append(np.deg2rad(112.5))
        self.names.append('betaC-gammaCr-deltaCr')
        self.theta.append(np.deg2rad(112.0))
        self.names.append('betaC-gammaC1i-deltaCi')
        self.theta.append(np.deg2rad(114.3))
        self.names.append('gammaCm-deltaSm-Tm')
        self.theta.append(np.deg2rad(100.9))
        self.names.append('deltaCk-epsilonCk-Tk')
        self.theta.append(np.deg2rad(55.7))
        self.names.append('deltaCr-epsilonNr-Tr')
        self.theta.append(np.deg2rad(110.9))
        self.names.append('gammaCk-deltaCk-epsilonCk')
        self.theta.append(np.deg2rad(112.4))
        self.names.append('gammaCr-deltaCr-epsilonNr')
        self.theta.append(np.deg2rad(111.3))
        self.names.append('gammaC1i-betaC-gammaC2i')
        self.theta.append(np.deg2rad(110.8))
        self.names.append('betaC-Sc-Sc')
        self.theta.append(np.deg2rad(121.4))

    def set_all_harmonic_angles(self, system):

        self.angle_ref = hoomd.md.angle.harmonic()
        self.add_to_logger()
        snap = system.take_snapshot(all=True)
        for a in snap.angles.types:
        #for a in system.angles:
            name = str(a)
            self.angle_ref.angle_coeff.set(name, k=self.k, t0=self.theta[self.names.index(name)])
        return self.angle_ref

