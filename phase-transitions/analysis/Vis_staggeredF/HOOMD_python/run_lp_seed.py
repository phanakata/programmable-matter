#!/usr/bin/python
import numpy as np
import hoomd
from hoomd import md, deprecated, data, group, init
import gsd.hoomd
import sys

# Initialization function
def set_random(L,Np):
  s = gsd.hoomd.Snapshot()
  s.particles.N = Np
  s.particles.types = ['A']
  s.configuration.dimensions = 2
  s.particles.typeid = []
  s.particles.position = []
  s.particles.diameter = []
  for i in range(Np):
    s.particles.typeid.append(0)
    s.particles.diameter.append(rm)
    y = (np.random.uniform() - 0.5)*L
    x = (np.random.uniform() - 0.5)*L
    z = 0.
    s.particles.position.append((x,y,z)) 
  s.configuration.box = [L, L,10,0,0,0]
  fname = "init_strip.gsd"
  gsd.hoomd.create(fname,snapshot=s)
  print "Strip initialized in file: %s" % fname
  return fname

# Parameters


if len(sys.argv) < 5:
  print("Usage: python %s [phi] [L] [lp] [tf]" % sys.argv[0])
  exit()


# Parameters
r0 = 1 
rm = 2*r0
L = float(sys.argv[2])
tf = float(sys.argv[4]) 
tstep = 1e-5
nsteps = tf / tstep
pt = 1
pstep = pt / tstep
phi= float(sys.argv[1]) 
Np = int( phi * L*L / (np.pi*r0**2) )
v0 = 1e2
Dr = v0/float(sys.argv[3])
Dt = 0
mu = 1
k  = 1e4
outgsd = "L%03i_config.gsd" % int(L/r0)

# User output
print("Parameters:")
print("   r0 = %4.2e" % (r0))
print("   rm = %4.2e" % (rm/r0))
print("    L = %4.2e" % (L/r0))
print("  phi = %4.2e" % (phi))
print("    N = %4.2e" % (Np))
print("   v0 = %4.2e" % (v0))
print("   Dr = %4.2e" % (Dr))
print("   Dt = %4.2e" % (Dt))
print("    k = %4.2e" % (k))
print("   mu = %4.2e" % (mu))

# Intialize context and configuration
hoomd.context.initialize()
#system = hoomd.init.read_gsd(initialize_strip(L,L,Np,rm), restart=outgsd)
#deprecated.init.create_random(N=Np,box=data.boxdim(Lx=L,Ly=L,dimensions=2),min_dist=2) 
system = hoomd.init.read_gsd(set_random(L,Np), restart=outgsd)

# Initialize neighborlist and interaction potential, WCA
nl = md.nlist.cell()
spring = md.pair.dpd_conservative(r_cut=rm, nlist=nl)
spring.pair_coeff.set('A', 'A', A=(rm*k))

# Define integration mode, time step, method
seedthermal=np.random.randint(2**16-1)
bd = md.integrate.brownian( group = group.all(), 
                               kT = Dt,
                             seed = seedthermal,
                           dscale = 1,
                      noiseless_t = False,
                      noiseless_r = True)

# Define output files for saving configurations: dcd, gsd
hoomd.dump.gsd(filename=outgsd, period=pstep, group=group.all(), overwrite=True, truncate=False, phase=0, time_step=None,static=['attribute','topology'])
#hoomd.dump.dcd(filename=outdcd, period=pstep, overwrite=True)
#xml = hoomd.deprecated.dump.xml(group=group.all(), filename=outxml, period=pstep)
#xml.set_params(position=True,image=True,orientation=True,velocity=True)


#thermalize
print("\n\n    Thermalizing...\n\n")
dt_therm=1e-6
md.integrate.mode_standard(dt=dt_therm)
hoomd.run(1./dt_therm)
print("... done.")


# Initialize and define Swim force
e = []
for i in range(Np):
  angle = 2*np.pi*np.random.rand()
  ex = np.cos(angle)
  ey = np.sin(angle)
  ez = 0
  e.append((ex,ey,ez)) 
Fs = [tuple(v0/mu*e_ij for e_ij in e[i]) for i in range(Np)]
seedswim=np.random.randint(2**16-1)
swim_force = hoomd.md.force.active( seed = seedswim,
                                   f_lst = Fs, 
                                   group = group.all(), 
                        orientation_link = False, 
                orientation_reverse_link = True, 
                           rotation_diff = Dr, 
                              constraint = None )

f=open("seed.dat","w")
f.write("seedswim: %d\n" % seedswim)
f.write("seedthermal: %d\n" % seedthermal)
f.close()

print("\n\n    Running...\n\n")
# Run simulation
md.integrate.mode_standard(dt=tstep)
hoomd.run(nsteps)
