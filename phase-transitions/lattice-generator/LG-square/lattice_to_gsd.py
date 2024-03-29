#!/usr/bin/python
import numpy as np
import hoomd
from hoomd import md, deprecated, data, group, init
import gsd.hoomd
#import gsd.fl
import sys

# Initialization function
def lattice(L,input_file,output_file):
  s = gsd.hoomd.Snapshot()
  s.particles.types = ['A','B','C','D','E']
  s.configuration.dimensions = 3
  s.particles.typeid = []
  s.particles.position = []
  s.particles.velocity = []
  s.particles.accelaration = []
  s.bonds.types = ['A', 'B', 'C', 'D']
  s.bonds.typeid = []
  s.bonds.group = []
  s.dihedrals.types = ['A']
  s.dihedrals.typeid = []
  s.dihedrals.group = []
  
  #	Reading configuration from lattice.dat file
  array = []
  #lines = [line.rstrip('\n') for line in open('../Sim_dump/lattice.dat')]
  lines = [line.rstrip('\n') for line in open(input_file)]
  for l in lines:
      array.append(l)

  s.particles.N = int(array[0])
  
  for i in range(s.particles.N):
    c = [float(e) for e in array[i+1].split(",")]
    x,y,z = c
    #s.particles.typeid.append(0)
    s.particles.position.append((x,y,z))
 
  for i in range(s.particles.N):
    s.particles.velocity.append((0,0,0))
    s.particles.accelaration.append((1,1,1))
 
  ptr = 1+s.particles.N
  #print((array[ptr]))
  s.bonds.N = int(array[ptr])
 
  for i in range(s.bonds.N):
    c = [int(e) for e in array[ptr+1+i].split(",")]
    p1,p2,p3= c
    #print (p3)
    s.bonds.group.append((p1,p2))	
    #s.bonds.typeid.append((p3))
  
  ptr2 = ptr
  #for i in range(s.bonds.N):
  #  c = [int(e) for e in array[ptr+1+i].split(",")]
  #  p1,p2,p3= c
    #print (p3)
  #  s.bonds.typeid.append((p3))
  
  #print(s.bonds.group)
  #print(s.bonds.typeid)
  #sprin.bonds.typeid[0] = 2
  ptr = ptr+s.bonds.N+1
  #print((array[ptr]))
  s.dihedrals.N = int(array[ptr])
  for i in range(s.dihedrals.N):
     c = [int(e) for e in array[ptr+1+i].split(",")]
     d1,d2,d3,d4 = c
     s.dihedrals.group.append((d1,d2,d3,d4)) 

  ptr = ptr+s.dihedrals.N+1
  print((array[ptr]))
  for i in range(s.particles.N):
     s.particles.typeid.append(int(array[ptr+i]))
  
  for i in range(s.bonds.N):
    c = [int(e) for e in array[ptr2+1+i].split(",")]
    p1,p2,p3= c
    #print (p3)
    s.bonds.typeid.append((p3))
  #for i in range(s.bonds.N):
  #   s.bonds.typeid.append((0))
  
  for i in range(s.dihedrals.N):
     s.dihedrals.typeid.append((0))

  s.configuration.box = [L[0],L[1],L[2],0,0,0]
  #fname = "../Sim_dump/init_strip.gsd"
  fname = output_file;
  gsd.hoomd.create(fname,snapshot=s)
  print ("Strip initialized in file: %s" % fname)
  return fname

# Parameters


if len(sys.argv) < 6:
  print("Usage: python %s [L] " % sys.argv[0])
  print("Usage: python %s [L] " % sys.argv[1])
  print("Usage: python %s [L] " % sys.argv[2])
  print("Usage: python %s [input_file] " % sys.argv[3])
  print("Usage: python %s [output_file] " % sys.argv[4])
  exit()


# Parameters
L=[]
L.append(float(sys.argv[1]))
L.append(float(sys.argv[2]))
L.append(float(sys.argv[3]))
input_file = sys.argv[4]
output_file = sys.argv[5]


# User output
print("Parameters:")
print("    L = %4.2e" % (L[0]))

lattice(L,input_file,output_file)
