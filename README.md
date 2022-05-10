# MultiFEBE

What is MultiFEBE?
==================

MultiFEBE is a multi-domain integrated Finite Element and Boundary Element solver.

How to install?
===============

## Option 1
Download binaries from [our website](http://www.mmc.siani.es), and run it.

## Option 2
Download the source code from [GitHub](https://github.com/mmc-siani-es/MultiFEBE), unpack the .zip file, and compile it by running the Makefile. 

Requirements:

  * Operating System: GNU/Linux 64 bits
  * Build automation tool: GNU Make
  * Compiler: GNU Fortran Compiler versión 9.4 o superior
  * Linear Algebra:
    * Default: Automatically Tuned Linear Algebra Software (ATLAS)
    * Alternative (need to change Makefile): OpenBLAS 

Optional but recommended:
  * Pre- and post-processor (external tool): Gmsh (https://gmsh.info/)

Instructions:

  * Download source code from [GitHub](https://github.com/mmc-siani-es/MultiFEBE).
  * Unpack the downloaded .zip file.
  * Review the requirements:
    * Install GNU Make and GNU Fortran:
```
$ sudo apt-get install make gfortran
``` 
   * Install ATLAS:
```
$ sudo apt-get install libatlas3-base libatlas-base-dev

```
   * Install OpenBlas (alternative to ATLAS):
```
$ sudo apt-get install libopenblas-base libopenblas-dev
```
  * Edit the `./Makefile` if you want to use OpenBlas instead of ATLAS.
  * Compile it by executing:
  
```
$ make   
```    

  * Once compiled, the executable `multifebe` is locate at the folder `./bin/` 
    
License:
========

GPLv2, please see the [LICENSE](https://github.com/mmc-siani-es/MultiFEBE/blob/main/LICENSE) file for details.

Financing:
========

This work has been developed with the support of research projects:

  * PID2020-120102RB-I00, funded by the Agencial Estatal de Investigación of Spain, MCIN/AEI/10.13039/501100011033.


  * ProID2020010025, funded by Consejerı́a de Economı́a, Conocimiento y Empleo (Agencia Canaria de la Investigación. Innovación y Sociedad de la Información) of the Gobierno de Canarias and FEDER;


  * BIA2017-88770-R, funded by Subdirección General de Proyectos de Investigación of the Ministerio de Economı́a y Competitividad (MINECO) of Spain and FEDER.
