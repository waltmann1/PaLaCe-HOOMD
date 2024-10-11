from __future__ import division
import numpy as np
import hoomd
from hoomd import md


class Rigid(object):

    def __init__(self):
       self.name = 'rigid'
       self.rigid = None

    def is_center(self, name):
        #print(name)
        if len(name) < 6:
            return False
        if name[:6] != 'center':
            return False
        return True

    def set_rigid_bodies(self, system):

        self.rigid = hoomd.md.constrain.rigid()
        rigids = 0
        while self.is_center(str(system.particles.get(rigids).type)):
            rigids += 1
        type_list = [[] for _ in range(rigids)]
        pos_list = [[] for _ in range(rigids)]
        for i in range(rigids, len(system.particles)):
            part = system.particles.get(i)
            if part.body > -1 and part.body < rigids:
                pos_list[part.body].append(np.array(system.particles.get(i).position))
                type_list[part.body].append(str(part.type))
        for i in range(rigids):
            #print('parmas')
            #print(pos_list[i])
            #print(type_list[i])
            this_pos = np.array(system.particles.get(i).position)
            self.rigid.set_param(system.particles.get(i).type, positions=np.subtract(pos_list[i],this_pos), types=type_list[i])
        self.rigid.create_bodies(create=False)
        # quit()
        return self.rigid

    def get_rigids(self, system):
        rigids = 0
        while self.is_center(str(system.particles.get(rigids).type)):
            rigids += 1
        return rigids
