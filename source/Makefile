###############################
#
# GNU complier (gcc version 5.1 (GCC) from Fink)
#
FC= gfortran
FFLAGS= -m64  -funroll-loops -O3

CXX= g++-8
CXXFLAGS=  -m64  -funroll-loops -O3  -ansi -pedantic -Wno-write-strings

LINK= -m64  -O3 

#
#  The ‘/sw/lib’ is needed only if you are using gcc compilers from the
# fink distribution (Mac OSX), otherwise is to be removed. If the LAPACK
# library is not installed in your system the relevant routines
# are included in the attached file ‘BcDor-LAPACKdiags.f’: in this case
# you will not need t link to -llapack  but will need the fortran compiler
# and the -lgfortran library (if you are using gcc) to link together with
# the C++.  In either case, BLAS is still needed for the CCD routines.
#
#  If you have LAPACK installed in your system, use that one.
#  
LINK2= -L/sw/lib   -lblas -llapack 
##LINK2= -L/sw/lib  -lgfortran 




#BcDor-Main.o 
OBJECTS= BcDor-Command_line.o      BcDor-Initialize.o \
         BcDor-Global_variables.o  BcDor-Run_vars.o \
         BcDor-DysonMatrix.o   BcDor-Koltun_MBPT.o  BcDor-Radii.o  BcDor-Miscellanea.o \
         BcDor-Ang_momenta.o  BcDor-Utilities.o \
	 BcDor-Mod_Space.o            BcDor-Mod_Space_UTILS.o \
	 BcDor-Vpp_potential.o        BcDor-Vpp_potential_c.o \
	 BcDor-Sp_Propagator.o        BcDor-Sp_Propagator_a.o \
	 BcDor-Self_energy_ext_Op_a.o BcDor-Self_energy_ext_Op_b.o \
	 BcDor-Self_energy_ext_Op_c.o BcDor-Self_energy_ext_Op_d.o  \
	 BcDor-Self_energy_ext_Op_e.o \
         BcDor-CCD_a.o      BcDor-CCD_b.o \
         BcDor-gII_a.o      BcDor-gII_b.o \
         BcDor-PolProp_a.o  BcDor-PolProp_b.o \
	 BcDor-Lanczos.o 
###	 BcDor-LAPACKdiags.o  ## You will need to uncomment this only if your system does not have LAPACK installed.


BcDor: BcDor-Main.cpp $(OBJECTS)
	$(CXX) $(LINK)  BcDor-Main.cpp  $(OBJECTS) $(LINK2)  -o BcDor
	./BcDor -h |a2ps -1 -r -l137 -o BcDor_help.ps
	ps2pdf BcDor_help.ps
# If you have ‘a2ps’ and ‘ps2pdf’ installed in your system, the last two lines 
# above will generate a print out of the help manual. Otherwise, it is save to 
# either comment them or just ignore the error message.
# 


clean:
	rm -f $(OBJECTS)

veryclean:
	rm -f ./*.o  ./*.bck  ./*.mod

