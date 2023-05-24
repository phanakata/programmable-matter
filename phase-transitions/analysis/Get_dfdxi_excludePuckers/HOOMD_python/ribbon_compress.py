#!/usr/bin/python
import numpy as np
import hoomd
from hoomd import md, deprecated, data, group, init
import gsd.hoomd
import sys
import random

k=5.000 #Kappa
e=1440.000*k/2 #Epsilon
nx=101
ny=21
Run=1
s=96.80437093 #Thermalized avg position of slider
d=2.93372430 #Thermalized compression
c=0.2 #Compress ribbon
i=0 #counter to iterate over lines in the last frame position file

print ("Kappa = ",k,"Epsilon = ",e)
obser_file = '../../Sim_dump_ribbon/obser_L'+str(nx)+'W_'+str(ny)+'_e'+str(e)+'_k'+str(k)+'_c'+str(c)+'_r'+str(Run)+'.log'
print ("Observable data is dumped in file: ",obser_file)
traj_file = '../../Sim_dump_ribbon/traj_L'+str(nx)+'W_'+str(ny)+'_e'+str(e)+'_k'+str(k)+'_c'+str(c)+'_r'+str(Run)+'.gsd'
init_strip = '../../Sim_dump_ribbon/init_strip_L'+str(nx)+'_W'+str(ny)+'.gsd'
print ("Initial init strip GSD file :",init_strip)

hoomd.context.initialize()
s = hoomd.init.read_gsd(init_strip)

#Reading lines from file containing frame data (thermalized ribbon) left end free sliding
thermal_ribbon = '../../Sim_dump_ribbon/'
lines = [line.rstrip('\n') for line in open('/home/sourav/Sim_dump_ribbon/test.dat')]

for p in s.particles:
	#print (lines[i])
	if (p.type != 'B' or p.type != 'D'):
		pos = lines[i].split()	
		x = float(pos[0])
		y = float(pos[1])
		z = float(pos[2])
		p.position = (x,y,z)
	if p.type == 'D':
		x, y, z = p.position
		x = x-(d+c)
		p.position = (x,y,z)
	i=i+1

#for p in s.particles:
#	if (p.type == 'A' or p.type == 'E'):
#		x, y, z = p.position
#		z += random.uniform(-0.10,0.10)
#		p.position = (x,y,z)

#for p in s.particles:
#	if p.type == 'D':
#		x, y, z = p.position
#		x = x-(d+c)
#		p.position = (x,y,z)

harmonic = md.bond.harmonic()
dih = md.dihedral.harmonic()

dih.dihedral_coeff.set('A', k=k, d=1, n=1)
dih.dihedral_coeff.set('D', k=k, d=1, n=1)
dih.dihedral_coeff.set('E', k=k, d=1, n=1)

harmonic.bond_coeff.set('A', k=e, r0=1.0)
harmonic.bond_coeff.set('D', k=e, r0=1.0)
harmonic.bond_coeff.set('E', k=e, r0=1.0)

hoomd.analyze.log(filename=obser_file, quantities=["temperature", "potential_energy","bond_harmonic_energy","kinetic_energy","dihedral_harmonic_energy"], period=5000, header_prefix="#", overwrite=True)

md.integrate.mode_standard(dt=0.0010)

group1 = hoomd.group.type(name='group1', type='A')
group2 = hoomd.group.type(name='group2', type='D')
group3 = hoomd.group.type(name='group3', type='E')
group13 = hoomd.group.union(name='group13',a=group1,b=group3)
group12 = hoomd.group.union(name='group12',a=group1,b=group2)
group123 = hoomd.group.union(name='group123',a=group12,b=group3)

#md.constrain.oneD(group=group2, constraint_vector=[1,0,0])

hoomd.dump.gsd(filename=traj_file, group=group.all(), period=5000, overwrite=True)


md.integrate.nvt(group=group13,kT=1.0, tau=0.2)

hoomd.run(1e6)

