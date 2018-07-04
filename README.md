# PyESVEP

## Synopsis

This project contains *Python* code for *End-member-based Soil and Vegetation Energy Partitioning* (ESVEP) model for estimating sensible and latent heat flux (evapotranspiration) based on measurements of radiometric surface temperature. The code was developed for the [Sentinels for Evapotranspiration (SEN-ET) project](http://esa-sen4et.org/).

The project requires [PyTSEB](https://github.com/hectornieto/pyTSEB) and its dependencies to run.

## Code Example
### High-level example

You can run ESVEP with the scripts *ESVEP_local_image_main.py* , 
which will read an input configuration file (defaults is *Config_LocalImage.txt*). 
You can edit this configuration file or make a copy to fit your data and site characteristics and either run 
this script in a Python GUI or in a terminal shell:

`python ESVEP_local_image_main.py <configuration file>`  
where \<configuration file> points to a customized configuration file... leave it blank if you want to use the default 
file *Config_LocalImage.txt*


### Low-level example
You can run ESVEP model or any related process in python by importing the module *ESVEP* from the *pyESVEP* package. 

```python
import pyESVEP.ESVEP as ESVEP 
output=ESVEP.ESVEP(Tr_K, vza, Ta_K, u, ea, p, Sn_C, Sn_S, L_dn, LAI, emis_C, emis_S, z_0M, d_0, z_u, z_T)
```

You can type
`help(ESVEP.ESVEP)`
to understand better the inputs needed and the outputs returned

   
## Basic Contents
### High-level modules
- *.pyESVEP/PyESVEP.py*, class (subclass of PyTSEB class from PyTSEB project) for ESVEP scripting

- *ESVEP_local_image_main.py*, high level scripts for running ESVEP through a configuration file (*Config_LocalImage.txt*)

### Low-level modules
The low-level modules in this project are aimed at providing customisation and more flexibility in running TSEB. 
The following modules are included

- *.pyESVEP/ESVEP.py*

core functions for running ESVEP model. 


## Main Scientific References
- R. Tang and Z. L. Li, “An End-Member-Based Two-Source Approach for Estimating Land Surface Evapotranspiration From Remote Sensing Data,” IEEE Transactions on Geoscience and Remote Sensing, vol. 55, no. 10, pp. 5818–5832, Oct. 2017.



## Tests
The folder *./Input* contains examples for images for running ESVEP (*ExampleImage_\< variable >.tif*). Just run the high-level script with the configuration file provided by default and compare the resulting outputs with the files stored in *./Output/*

## Contributors
- **Radoslaw Guzinski** main developer, tester
- **Hector Nieto** main developer of PyTSEB, tester

## License
PyESVEP: a Python End-member-based Soil and Vegetation Energy Partitioning (ESVEP) model

Copyright 2018 Radoslaw Guzinski and contributors.
    
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
