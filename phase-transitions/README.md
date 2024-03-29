## Overview 
This project's directory structure is as follows:
* ```lattice-configurations/```: contains initial lattice configurations 
* ```lattice-generator/```: contains codes to generate lattice in dat and GSD format 
* ```planar-sheet/```: NPT and NPH simulations for planar sheets
* ```cylinder/```: NPT and NPH simulations for cylinders 
* ```analysis/```: contains analysis codes (C++ and python) 




## To download 
git clone https://github.com/phanakata/programmable-matter/phase-transitions

## Technical 
* Molecular Dynamics Simulations 
  * HOOMD-blue 2.8.1 CUDA (10.1) DOUBLE HPMC_MIXED TBB SSE SSE2 
  * HOOMD run either on CPU or Tesla V100-PCIE-32GB  80 SM_7.0 @ 1.38 GHz, 32510 MiB DRAM, MNG


## Authors
Paul Hanakata

## Citation

If you use this package/code/dataset, build on  or find our research is useful for your work please cite 
```
@article{PhysRevLett.128.075902,
  title = {Anomalous Thermal Expansion in Ising-like Puckered Sheets},
  author = {Hanakata, Paul Z. and Plummer, Abigail and Nelson, David R.},
  journal = {Phys. Rev. Lett.},
  volume = {128},
  issue = {7},
  pages = {075902},
  numpages = {7},
  year = {2022},
  month = {Feb},
  publisher = {American Physical Society},
  doi = {10.1103/PhysRevLett.128.075902},
  url = {https://link.aps.org/doi/10.1103/PhysRevLett.128.075902}
}
@article{PhysRevMaterials.6.115203,
  title = {Curvature as an external field in mechanical antiferromagnets},
  author = {Plummer, Abigail and Hanakata, Paul Z. and Nelson, David R.},
  journal = {Phys. Rev. Mater.},
  volume = {6},
  issue = {11},
  pages = {115203},
  numpages = {18},
  year = {2022},
  month = {Nov},
  publisher = {American Physical Society},
  doi = {10.1103/PhysRevMaterials.6.115203},
  url = {https://link.aps.org/doi/10.1103/PhysRevMaterials.6.115203}
}
```


## References
* <a href="https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.121.255304" style="color:#268cd7
"> **Paul Z. Hanakata**, Abigail P. Plummer, David R. Nelson , *Anomalous thermal expansion in Ising-like puckered sheets*, Phys. Rev. Lett, 128, 075902  (2022).</a>
* <a href="https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.6.115203" style="color:#268cd7
"> Abigail Plummer*, **Paul Z. Hanakata***, David R. Nelson , *Curvature as an external field in mechanical antiferromagnets*, Phys. Rev. Materials 6, 115203 (2022).</a>
