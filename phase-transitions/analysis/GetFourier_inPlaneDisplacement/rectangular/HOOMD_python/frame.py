#!/usr/bin/python
import numpy as np
import hoomd
from hoomd import md, deprecated, data, group, init
import gsd.hoomd
import sys
import random

hoomd.context.initialize()

s = hoomd.init.read_gsd("../Sim_dump_frame/init_strip.gsd")

for p in s.particles:
	if p.type == 'A':
		x, y, z = p.position
		z += random.uniform(-0.10,0.10)
		p.position = (x,y,z)

harmonic = md.bond.harmonic()
dih = md.dihedral.harmonic()

#walls.add_plane(origin=(50.5, 50.4, 0),normal=(-1,-1,0),inside=True)

#dih.dihedral_coeff.set('B', k=5.000, d=1, n=1)
#harmonic.bond_coeff.set('B', k=800.000, r0=1.0)

walls = md.wall.group()
walls.add_plane(origin=(-1,-1,0),normal=(1,1,0),inside=True)
print(walls)

dih.dihedral_coeff.set('A', k=10.000, d=1, n=1)
harmonic.bond_coeff.set('A', k=800.000, r0=1.0)

lj=md.wall.lj(walls, r_cut=75.0)
lj.force_coeff.set('B', sigma=1.0,epsilon=10.0, r_cut=75.0)
lj.force_coeff.set('A', sigma=1.0,epsilon=0.0, r_cut=25.0)
lj.force_coeff.set('C', sigma=1.0,epsilon=0.0, r_cut=25.0)
lj.force_coeff.set('D', sigma=1.0,epsilon=0.0, r_cut=25.0)


hoomd.analyze.log(filename="../Sim_dump_frame/observable4.log", quantities=["temperature", "potential_energy","bond_harmonic_energy","kinetic_energy","dihedral_harmonic_energy"], period=5000, header_prefix="#", overwrite=True)

md.integrate.mode_standard(dt=0.0010)

group1 = hoomd.group.type(name='group1', type='A')
group2 = hoomd.group.type(name='group2', type='B')
group12 = hoomd.group.union(name='group12',a=group1,b=group2)
group3 = hoomd.group.type(name='group1', type='D')
group4 = hoomd.group.union(name='group12',a=group12,b=group3)

md.constrain.oneD(group=group2, constraint_vector=[1,1,0])

hoomd.dump.gsd(filename="../Sim_dump_frame/trajectory4.gsd", group=group.all(), period=5000, overwrite=True,static=[])

#print(s.particles[2899])

#all = group.all()
md.integrate.nvt(group=group12,kT=1.0, tau=0.2)

#print(s.particles[5])

hoomd.run(1e7)

#print(s.particles[5])

#hoomd.run(10)
#print(s.particles[5])
#print(s.bonds[5])
