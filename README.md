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

  * OS: GNU/Linux 64 bits
  * Build automation tool: GNU Make
  * Compiler: GNU Fortran Compiler versión 9.4 o superior
  * Linear Algebra:
    * Default: Automatically Tuned Linear Algebra Software (ATLAS)
    * Other (need to change Makefile): OpenBLAS 
  * Pre- and post-processor: Gmsh (https://gmsh.info/)

Instructions:

  * Download source code from [GitHub](https://github.com/mmc-siani-es/MultiFEBE).
  * Unpack the downloaded .zip file.
  * Edit the `./Makefile` if needed.
  * Compile it by executing:
  
```
$ make   
```    

  * Once compiled, the executable `multifebe` is locate at the folder `./bin/` 
    
License:
========

GPLv2, please see the [LICENSE](https://github.com/mmc-siani-es/MultiFEBE/blob/main/LICENSE) file for details.
