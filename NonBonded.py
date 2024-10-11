from __future__ import division
import numpy as np
import hoomd
from hoomd import md
from numpy import linalg as la
from Loggable import Loggable


class Solvation(Loggable):

    def __init__(self, cv_min=5.0, cv_max=6.0, log_list=None):

        super(Solvation, self).__init__(log_list)
        self.log_values = ['pair_cvcel_electrostatics_energy']
        self.names = []
        self.Ni = []

        self.names = ['alphaC', 'alphaCp', 'alphaCg', 'Sa', 'Sc', 'Sd', 'Se', 'Sf', 'Sh', 'Si', 'Sk', 'Sl', 'Sm', 'Sn',
                      'Sp', 'Sq', 'Sr', 'Ss', 'St', 'Sv', 'Sw', 'Sy', 'Te', 'Tf', 'Th', 'Tk', 'Tm', 'Tq', 'Tr', 'Tw',
                      'Ty']

        self.Ni = [6, 6, 6, 4, 5, 6, 6, 3, 3, 13, 6, 13, 6, 8, 9, 6, 6, 5, 8, 10, 4, 3, 3, 11, 8, 10, 5, 5, 12, 14, 12]

        self.gamma = [.0133, -.4117]
        self.delta = [-.0103, .4493]
        self.cv_min = cv_min
        self.cv_max = cv_max
        self.cvcel = None
        self.charge_unit = 18.3 

    def add_all_solvation(self, system):
        """

        :param system: system object
        :return:
        """
        for i in range(len(system.particles)):
            part = system.particles.get(i)
            name = str(part.type)
            if name in self.names:
                one = self.gamma[0] * self.Ni[self.names.index(name)] + self.gamma[1]
                two = self.delta[0] * self.Ni[self.names.index(name)] + self.delta[1]
                part.sol = (one, two)
                part.charge *= self.charge_unit

    def set_all_solvation(self, neighbor_list, system):

        #.15M lamda dh = 7.9
        #.5 M lamda_dh = 4.3
        self.add_all_solvation(system)
        self.add_to_logger()
        length = min(min(system.box.Lz, system.box.Lx), system.box.Ly)
        r_cut = length / 2 - 1
        r_cut=20
        kd = 1 / 7.9
        # dielectric constant of protein / dielectric constant of water
        e_i = 4 / 78.5
        self.cvcel = hoomd.md.pair.solvation(r_cut=r_cut, nlist=neighbor_list, cv_max=self.cv_max,
                                             cv_min=self.cv_min)
        snap = system.take_snapshot()
        for t1 in snap.particles.types:
            for t2 in snap.particles.types:
                self.cvcel.pair_coeff.set(str(t1), str(t2), kappa_DH=kd, eps_in=e_i)
        return self.cvcel


class CoarseGrain(Loggable):
    def __init__(self, log_list=None):

        super(CoarseGrain, self).__init__(log_list)
        self.log_values = ['pair_lj86_energy', 'pair_table_energy']

        self.names = ['alphaC', 'alphaCp', 'alphaCg', 'Sa', 'Sc', 'Sd', 'Se', 'Sf', 'Sh', 'Si', 'Sk', 'Sl', 'Sm', 'Sn',
                      'Sp', 'Sq', 'Sr', 'Ss', 'St', 'Sv', 'Sw', 'Sy', 'Te', 'Tf', 'Th', 'Tk', 'Tm', 'Tq', 'Tr', 'Tw',
                      'Ty']

        self.epsilon = [.002, .002, .6, 2.0, 1.4, .1, 2.0, 1.0, 1.8, 1.2, 2.1, 2.0, .7, .2, 1.2, 1.1, .9, .5, .6, 1.0,
                        4.1, 1.5, .002, .2, .05, .1, .2, .003, .02, .02, .1]
        self.sigma = [9.0, 8.7, 4.2, 3.3, 4.3, 4.4, 4.0, 3.9, 3.8, 5.4, 4.1, 5.5, 4.8, 4.4, 4.8, 4.3, 4.4, 3.9, 4.7,
                      4.9, 3.3, 3.8, 6.9, 6.0, 6.5, 6.0, 5.1, 7.0, 5.6, 8.7, 5.7]
        self.charge = [0] * 31
        self.charge[self.names.index('Sd')] = -1
        self.charge[self.names.index('Tr')] = 1
        self.charge[self.names.index('Tk')] = 1
        self.charge[self.names.index('Te')] = -1

        self.charge_pair = None
        self.lj86_pair = None

    def set_nb_coarse_grain(self, neighbor_list, system):
        cut = 20
        self.lj86_pair = hoomd.md.pair.lj86(r_cut=cut, nlist=neighbor_list)
        self.charge_pair = hoomd.md.pair.table(width=1000, nlist=neighbor_list)
        self.add_to_logger()
        for t1 in system.particles.types:
            for t2 in system.particles.types:
                t1 = str(t1)
                t2 = str(t2)
                if t1 in self.names and t2 in self.names:
                    self.lj86_pair.pair_coeff.set(str(t1), str(t2), epsilon=np.sqrt(self.epsilon[self.names.index(t1)] *
                                                                  self.epsilon[self.names.index(t2)]),
                                          sigma=(self.sigma[self.names.index(t1)] +
                                                        self.sigma[self.names.index(t2)]) / 2)
                    self.charge_pair.pair_coeff.set(str(t1), str(t2), func=cg_charge_pair, rmin=0.0, rmax=cut, coeff=dict(
                                          q1=self.charge[self.names.index(str(t1))],
                                          q2=self.charge[self.names.index(str(t2))]))
                else:
                    self.lj86_pair.pair_coeff.set(str(t1), str(t2), epsilon=0.0, sigma=0.0)
                    self.charge_pair.pair_coeff.set(str(t1), str(t2), func=cg_charge_pair, rmin=0.0, rmax=cut,
                                          coeff=dict(q1=0, q2=0))
        return self.lj86_pair

    def write_parameter_table(self):

        f = open('Coarse_Grain_parameter_table.txt', 'w')
        f.write("name   sigma    epsilon\n")
        for i in range(len(self.names)):
            f.write(self.names[i] + "   " + str(self.sigma[i]) + "    " + str(self.epsilon[i]) + '\n')


class HydrogenBonding(Loggable):

    def __init__(self, log_list=None):

        super(HydrogenBonding, self).__init__(log_list)
        self.log_values = ['pair_lj_energy', 'pair_table_energy']
        self.names = ['O', 'H',  'C', 'N', 'tC', 'tN']
        self.q = np.multiply([-0.5679 * 1.5/1, .2719 * 1.5/1, .5973 / 2, -.4157/ 2, .5973/ 2, -.4157/ 2], 1)
        self.epsilon = [.21, .0157, .086, .17, .086, .17]
        self.sigma = [1.66, 0.6, 1.91, 1.82, 1.91, 1.82]
        #self.epsilon_o = .21
        #self.sigma_o = 1.6612
        #self.epsilon_h = .0157
        #self.sigma_h = .6
        self.charge_pair = None
        self.lj_pair = None

    def set_all_hb(self, neighbor_list, system):
        cut = 10
        #hoomd.md.pair.Geometric()
        self.lj_pair = hoomd.md.pair.lj(r_cut=cut, nlist=neighbor_list, geo=False)
        self.charge_pair = hoomd.md.pair.table(width=1000, nlist=neighbor_list)
        self.add_to_logger()
        for t1 in system.particles.types:
            for t2 in system.particles.types:
                if t1 in self.names and t2 in self.names:
                    index1 = self.names.index(str(t1))
                    index2 = self.names.index(str(t2))
                    self.lj_pair.pair_coeff.set(str(t1), str(t2), epsilon=np.sqrt(self.epsilon[index1] * self.epsilon[index2]),
                                      sigma=(self.sigma[index1] + self.sigma[index2]) / 1)
                    self.charge_pair.pair_coeff.set(str(t1), str(t2), func=hb_charge_pair, rmin=0.0, rmax=cut,
                                      coeff=dict(q1=self.q[index1], q2=self.q[index2]))
                    #print(str(t1),self.q[index1] )
                    #print(str(t2),self.q[index2] )
                else:
                    self.lj_pair.pair_coeff.set(str(t1), str(t2), epsilon=0.0, sigma=0.0)
                    self.charge_pair.pair_coeff.set(str(t1), str(t2), func=hb_charge_pair, rmin=0.0, rmax=cut,
                                      coeff=dict(q1=0.0, q2=0.0))
        return self.lj_pair


class Constraints(object):

    def __init__(self):
        self.name = 'constraints'
        self.constraint = None
        self.m = 0

    def add_constraints(self, system, rigids=0, geo=True):

        if geo:
            hoomd.md.pair.Geometric()
            return
        self.constraint = hoomd.md.constrain.distance()
        self.constraint.set_params(rel_tol=.01)
        #print(list(range(rigids)))
        if len(system.constraints) != 0:
            return
        for ind, part in enumerate(system.particles):
            if part.type == 'H':
                print(part.body, part.type, part.body in list(range(rigids)))
            if part.type == 'H' and part.body not in list(range(rigids)):
                dist_nh = la.norm(np.subtract(list(part.position), list(system.particles.get(ind - 1).position)))
                system.constraints.add(ind, ind - 1, dist_nh)
                system.bonds.add('exclusion', ind, ind + 2)
                system.bonds.add('exclusion', ind, ind + 3)
                #dist_ach = la.norm(np.subtract(list(part.position), list(system.particles.get(ind + 1).position)))
                #system.constraints.add(ind, ind + 1, dist_ach)
                if ind - 2 > -1 and system.particles.get(ind - 2).body not in list(range(rigids)):
                    if str(system.particles.get(ind - 2).type) == 'O':
                        dist_oh = la.norm(np.subtract(list(part.position), list(system.particles.get(ind - 2).position)))
                        system.constraints.add(ind, ind - 2, dist_oh)
                        system.bonds.add('exclusion', ind, ind-2)
                        system.bonds.add('exclusion', ind, ind-3)
                        system.bonds.add('exclusion', ind-1, ind-3)
                        system.bonds.add('exclusion', ind - 2, ind-6)
            elif part.type == 'O' and part.body not in list(range(rigids)):
                dist_co = la.norm(np.subtract(list(part.position), list(system.particles.get(ind - 1).position)))
                print('adding the damn constraint')
                system.constraints.add(ind, ind - 1, dist_co)
                #system.bonds.add('C-O', ind, ind-1)
                #dist_aco = la.norm(np.subtract(list(part.position), list(system.particles.get(ind - 2).position)))
                #system.constraints.add(ind, ind - 2, dist_aco)
        #quit()

    def add_force_transfer(self):
        hoomd.md.pair.GeometricTransfer()
        return


def hb_charge_pair(r, rmin, rmax, q1, q2):
    if r == 0:
        return 0, 0
    #335 in units of kcal /mol
    V = 335 * q1 * q2 / (epsilon_hb(r) * r) #kcal/mol
    #print(q1, q2, r,V)
    F = 335 * q1 * q2 * ((r * espilon_hb_prime(r) + epsilon_hb(r)) / (r**2 * epsilon_hb(r)**2))
    return V, F


def epsilon_hb(r):
    rj = 5.5
    kq = .43
    if r < rj:
        eps = kq * r**2 + 1
    else:
        eps = 2 * kq * r * rj - kq * rj**2 + 1
    return eps


def espilon_hb_prime(r):
    rj = 5.5
    kq = .43
    if r < rj:
        eps = kq * r*2
    else:
        eps = 2 * kq * rj
    return eps


def cg_charge_pair(r, rmin, rmax, q1, q2):
    if r == 0:
        return 0, 0
    V = 335 * q1 * q2 / (epsilon_cg(r) * r) #kcal/mol
    F = 335 * q1 * q2 * ((r * epsilon_cg_prime(r) + epsilon_cg(r)) / (r**2 * epsilon_cg(r)**2))
    #print(V,F)
    return V, F


def epsilon_cg(r):
    ke = 7.99
    return ke * r


def epsilon_cg_prime(r):
    ke = 7.99
    return ke