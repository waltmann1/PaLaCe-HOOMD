from __future__ import division
import numpy as np
from Simulation import Simulation
from ProteinComplex import ProteinComplex
from ProteinBackbone import ProteinBackbone
from Utils import all_amino_acids
import copy as cp


y = ProteinComplex('1sva.pdb')
y.dump_backbone_pdb()
quit()
#y.dump_resmap()
#quit()
#y = ProteinBackbone(30)
#y.add_all("GLY")
#for ind, chain in enumerate(y.chains):
#	print(ind)
#	for ind2, typ in enumerate(chain.type):
#		if typ == 'tN':
#			print(ind2,'tN')
#		if typ == 'tC':
#			print(ind2, 'tC')
		
#y.dump_xyz("check.xyz")
#print(y.ss_angles)
#quit()
#print(len(y.chains))
#quit(

#y = ProteinBackbone(200)

"""
z = ProteinBackbone(1)
z.add_all('GLY')
#while(len(y.chains) > 1):
#   y.remove_chain(y.chains[-1])

y.add_chain_to_N_terminal(z)
print(y.get_dihedrals_coupling())
print(y.get_dihedrals())
print(y.get_angles())
print(y.get_bonds())
print(y.res_indexes)
print(y.num_particles)

for i in y.get_dihedrals_coupling()[-1]:
   print(y.type[i])

print("")
for i in y.get_dihedrals()[-1]:
   print(y.type[i])
print("")
for i in y.get_angles()[-1]:
      print(y.type[i])
print("")
for i in y.get_bonds()[-1]:
         print(y.type[i])

quit()

y.pair_config(random_orientation=True)
"""
y.make_connections_rigid()
#y.dump_resmap()


#y.set_gsd_positions('protein_sim.gsd', 99)




#quit()



sim = Simulation(y.create_system())

sim.temp_interp(10000)
sim.run(50000)
