## lattice-generator

To run HOOMD simulation you need provide initial positions of nodes. You can either create within HOOMD or create your owen cofiguration generator. 
Lattice generator for square lattice (LG-square) with defects are stable. The rests are under development. 

## directory 
```lattice-generator/```:
* ```configurations/```: (contain dat files for gsd initial configuration + gsd files)
* ```LG-square/```:lattice generation for square + defects (STABLE)
* ```LG-tringular-pbc/```: lattice generation for triangular lattice (under development)
* ```LG-sphere/```: lattice generation for sphere (not developed yet)

## LG-square
1. Use 'make' to compile cpp codes
2. Execute LATINITIAL, Usage: ./LATINITIAL NX NY pbc defect_type spacing_x spacing_y
3. Output file is generated : ../configurations/lattice-file.dat
4. Once lattice-file.dat is generated use python code (lattice_to_gsd.py) to generate GSD file.
5. python lattice_to_gsd.py INPUT PARAMETERS ../configurations/lattice-file.dat ../configurations/init-file.gsd
