# MultiFEBE

What is MultiFEBE?
==================

MultiFEBE is a multi-domain Finite Element and Boundary Element linear mechanics solver. It implements mixed-dimensional couplings between finite elements and boundary elements. It is available for Windows and GNU/Linux 64 bits.

How to install?
===============

## Option 1 (installation)
Download an installer from the [releases](https://github.com/mmc-siani-es/MultiFEBE/releases) and run it. 
More details about the installation can be found in the [manual](https://github.com/mmc-siani-es/MultiFEBE/blob/main/docs/manual.pdf).

## Option 2 (compilation)
Download the source code from [GitHub](https://github.com/mmc-siani-es/MultiFEBE) or [packages](https://github.com/mmc-siani-es/MultiFEBE/packages), and compile it by running the CMake and GNU Make. 
More details about the compilation can be found in Appendix A of the [manual](https://github.com/mmc-siani-es/MultiFEBE/blob/main/docs/manual.pdf).

Requirements:

  * Operating System: Windows or GNU/Linux 64 bits.
  * Build automation tool: CMake, GNU Make.
  * Compiler: GNU Fortran Compiler version 9.4 or superior.
  * Dependencies: OpenBLAS (Linear Algebra).
  * For compilation in Windows: [MSYS2 environment](https://www.msys2.org).

How to use?
===========

MultiFEBE is a command-line program. It takes an input file which defines the case, and, once the case is solved, it creates output files with the results. 

The input file contains data sections defining the type of analysis, mesh, materials, boundary conditions, modelling details, etc. 
It may contain the whole mesh definition (nodes, elements, parts), but it is more appropriate to use a reference to a separate mesh file.
Mesh files can be formatted in a native format, or in [Gmsh](https://gmsh.info) MSH mesh file format 2.2. Results can also be exported to Gmsh file format, so you can fully use Gmsh as pre- and post-processor for MultiFEBE. 
The native mesh file format is very easy to build manually or from custom-built scripts. In this sense, we have created a .bas template file which can be used in [GiD](https://www.gidsimulation.com/) pre- and post-processor to generate our native mesh files.
More details can be found in the [manual](https://github.com/mmc-siani-es/MultiFEBE/blob/main/docs/manual.pdf).

License:
========

GPLv2, please see the [LICENSE](https://github.com/mmc-siani-es/MultiFEBE/blob/main/LICENSE) file for details.

Financing:
========

This work has been developed with the support of research projects:

  * PID2020-120102RB-I00, funded by the Agencial Estatal de Investigación of Spain, MCIN/AEI/10.13039/501100011033.

  <p align="center">
    <img src="docs/img/miciinn-aei.png">
  </p>  

  * ProID2020010025, funded by Consejerı́a de Economı́a, Conocimiento y Empleo (Agencia Canaria de la Investigación. Innovación y Sociedad de la Información) of the Gobierno de Canarias and FEDER;

  <p align="center">
    <img src="docs/img/gobcan-fse.png">
  </p>

  * BIA2017-88770-R, funded by Subdirección General de Proyectos de Investigación of the Ministerio de Economı́a y Competitividad (MINECO) of Spain and FEDER.
  
  <p align="center">
    <img src="docs/img/miciinn-feder-aei.png">
  </p> 

