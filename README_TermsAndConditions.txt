

TERMS and CONDITIONS
====================

  The files contained in this folders are the public distribution of the
BoccaDorata software for many-body nuclear structure. This is meant for
performing microscopic calculations with Green’s function and coupled cluster
theory for finite nuclei and with modern realistic interactions.
  This is freely distributed software and you are welcome to use it for your
own research and teaching, provided you abide to the following conditions:

 - The software is provided as it is, without any guarantee or commitment
  for support.
 
 - You are welcome to use and redistribute this software, as long as it is
  distributed as a whole and this file is always included and unmodified 
  from its original version.  If the software is redistributed with any
  modification form its original form, this should be clearly stated.
 
 - Any publication resulting from the use of the codes contained in this 
  distribution must acknowledge their use. If a modified version of the software
  has been used, this should be mentioned.

 - When using the software for research purposes, I kindly ask you to consider
  referring to any of the publications which led its development:
   Prog. Part. Nucl. Phys. 52, p. 377 (2004),
   Phys. Rev. A76, 052503 (2007),
   Phys. Rev. C79, 064313 (2009),
   Phys. Rev. C89, 024323 (2014).


INSTALLATION
============

  The source code is found in the folder ‘source’. Enter this directory, 
modify the ‘Makefile’ appropriately and type ‘make’. Put the executable
‘BcDor’ in your search $PATH. A C++ compiler and a BLAS/LAPACK library
are needed.

  Once the code is compiled type:

BcDor -v    (for the version and compilation time)

BcDor -h    (for a quick help)


DOCUMENTATION and EXAMPLES
==========================

  A manual of the BoccaDorata package can be found in the ‘docs’ folder.
This includes basics information, the description of file formats requested
for for input and a quick getting started tutorial.

  Auxiliary files for COM corrections and the Coulomb interaction in a 
harmonic oscillator basis are found the the folder ‘datafiles’.

  Some examples of calculated propagators are in the ‘examples’ folder.


HAPPY COMPUTING!


(c) Carlo Barbieri, Surrey, July 2015.