import numpy as np
import hoomd
from hoomd import md, deprecated, data, group, init
#import gsd.hoomd
import sys
import random

#elastic parameters
kappa=2.0 #Kappa
kstretch=100 #harmonic bond 
#system parameters
nx=48 #LX
ny=48 #LY 
nz=48 #LZ
spx=2 #spacing in x 
spy=2 #spacing in y 
LLX=nx 
LLY=ny 
LLZ=nz
#prestrain and wraped 
strain=0
wrap=0

#pucker/stitches spin=0/1; AFM/FM order=-1/1; flat/non-flat height/0  
spin=0   
order=-1 #AFM=-1, FM=1
height=0.45
frac=0.1
re_A=1.0
re_B=1.41421356237
re_C = re_A*(1+frac)
re_D = re_A*np.sqrt(2+2*frac+frac*frac)

#for MD 
Run=1       #differen runs 
noise=0
kin=0.100   #temperature

print ("kappa = ",kappa,"kstretch = ",kstretch)
obser_file = 'obser.log'
obser_file0 = 'relax.log'

print ("Observable data is dumped in file: ",obser_file)
traj_file = 'traj.gsd'
traj_file0 = 'relax.gsd'
#traj_fileLast = 'traj.gsd'

init_strip = '../../../../../lattice-configurations/square-sheets//init_pbc_square'+str(nx)+'_sp2_checker.gsd'

print ("Init_strip.gsd file : ",init_strip)

hoomd.context.initialize()
s = hoomd.init.read_gsd(init_strip)

#zero out all the z position 
for p in s.particles:
        x, y, z = p.position
        p.position = (x,y,0)

#initialize the "spins" in the flat phase  
for i in range(spin, nx, spx):
        for j in range(spin, ny, spy):
                x,y,z = s.particles[(j)*nx+(i)].position
                z = height*np.power(order, (i)/spx+(j)/spy)
                s.particles[(j)*nx+(i)].position =(x, y, z)



#change box size
#s.box = data.boxdim(Lx=LLX*(1+strain), Ly=LLX*(1+strain), Lz=LLZ, xy=0, xz=0, yz=0)
#print(s.box)

#shrink/expand the system by small amount
for p in s.particles:
    x, y, z = p.position
    x=(1+strain)*(x)      
    y=(1+strain)*(y)     
    p.position = (x,y,z)


#change box size
s.box = data.boxdim(Lx=LLX*(1+strain), Ly=LLX*(1+strain), Lz=LLZ, xy=0, xz=0, yz=0)
print(s.box)


#wrapped into a cylinder
if wrap==1:
        #update ny (L along the wrapping) and R 
        ny = ny*(1+strain) 
        R = ny/2.0/np.pi #radius 
 
        for p in s.particles:
                x, y, z = p.position
                y1=(R-z)*np.sin(2*np.pi*(y/ny))
                z1=-R+R*(1-np.cos(2*np.pi*(y/ny)))+z*np.cos(2*np.pi*(y/ny))
                p.position = (x,y1,z1)

#add random noise
if noise==1: 
        for p in s.particles:
                x, y, z = p.position
                x+=random.uniform(-0.02,0.02)
                y+=random.uniform(-0.02,0.02)
                z+=random.uniform(-0.02,0.02)
                p.position = (x,y,z)


harmonic = md.bond.harmonic()
dih = md.dihedral.harmonic()

dih.dihedral_coeff.set('A', k=kappa, d=1, n=1)
dih.dihedral_coeff.set('B', k=kappa, d=1, n=1)
dih.dihedral_coeff.set('C', k=kappa, d=1, n=1)
dih.dihedral_coeff.set('D', k=kappa, d=1, n=1)
dih.dihedral_coeff.set('E', k=kappa, d=1, n=1)

harmonic.bond_coeff.set('A', k=kstretch, r0=re_A)
harmonic.bond_coeff.set('B', k=kstretch, r0=re_B)
harmonic.bond_coeff.set('C', k=kstretch, r0=re_C)
harmonic.bond_coeff.set('D', k=kstretch, r0=re_D)

relax_log=hoomd.analyze.log(filename=obser_file0, quantities=["potential_energy","bond_harmonic_energy","dihedral_harmonic_energy", "temperature", "kinetic_energy", "lx", "ly", "pressure_xx", "pressure_yy"], period=100, header_prefix="#", overwrite=True)

relax_gsd=hoomd.dump.gsd(filename=traj_file0, group=group.all(), period=1000, overwrite=True)


#fire=hoomd.md.integrate.mode_minimize_fire(dt=0.005, ftol=1e-4, Etol=1e-10)
fire=hoomd.md.integrate.mode_minimize_fire(dt=0.005, ftol=1e-6, Etol=1e-10)
#NVE no box optimization 
#nve=hoomd.md.integrate.nve(group=group.all())
nph=hoomd.md.integrate.nph(group=group.all(),P=0.0,tauP=1.0, gamma=.5)
while not(fire.has_converged()):
   hoomd.run(1000, quiet=False);
nph.disable()
relax_log.disable()
relax_gsd.disable()
#fire.disable()

#hoomd.dump.gsd(filename=traj_fileLast, group=group.all(), period=1, overwrite=True)
#hoomd.run(1)
#nph.disable()

#print("Done!")


md.integrate.mode_standard(dt=0.0010)
hoomd.analyze.log(filename=obser_file, quantities=["potential_energy","bond_harmonic_energy","dihedral_harmonic_energy", "temperature", "kinetic_energy", "lx", "ly", "pressure_xx", "pressure_yy"], period=10000, header_prefix="#", overwrite=True)
hoomd.dump.gsd(filename=traj_file, group=group.all(), period=10000, overwrite=True)
npt = md.integrate.npt(group=group.all(),kT=kin, tau=0.2, tauP=1.0, P=0.0)
#npt.randomize_velocities(seed=Run)
hoomd.run(1.000001e7)
print("Done!")
