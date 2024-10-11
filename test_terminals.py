from ProteinComplex import ProteinComplex
from ProteinBackbone import ProteinBackbone


y = ProteinBackbone(200)
y.add_all('LYS')
y.dump_xyz('check1.xyz')
z = ProteinBackbone(1)
z.add_all('PRO')
y.add_chain_to_C_terminal(z)


y.dump_xyz("check.xyz")
print("dumped check")
