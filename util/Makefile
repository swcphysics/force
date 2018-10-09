# Bingwei Long  Dec 21/2017

include ../Makefile.common

UTLOBJS = util_cheb.o \
          util_bscint.o \
          util_gauleg.o \
          util_gadgets.o \
          util_inter_extra.o \
          util_io.o \
          util_nr.o \
          util_mesh.o \
          util_rootfinder.o \
          util_spcfnc.o \
          util_datetime.o

utlmods : $(UTLOBJS)

# ****WARNING!****: Must put the .f90 source as the first prerequisite!

util_gauleg.o : util_gauleg.f90 nneft_type.o
	$(FC) $(FFOPTS)

util_inter_extra.o : util_inter_extra.f90 nneft_type.o
	$(FC) $(FFOPTS)

util_mesh.o : util_mesh.f90 nneft_type.o util_gauleg.o
	$(FC) $(FFOPTS)

util_bscint.o : util_bscint.f90 nneft_type.o
	$(FC) $(FFOPTS)

util_gadgets.o : util_gadgets.f90 nneft_type.o nneft_phyconst.o
	$(FC) $(FFOPTS)

util_io.o : util_io.f90 nneft_type.o
	$(FC) $(FFOPTS)

util_spcfnc.o : util_spcfnc.f90 nneft_type.o
	$(FC) $(FFOPTS)

util_rootfinder.o : util_rootfinder.f90 nneft_type.o
	$(FC) $(FFOPTS)

util_nr.o : util_nr.f90 nneft_type.o
	$(FC) $(FFOPTS)

util_cheb.o	: util_cheb.f90 nneft_type.o util_nr.o
	$(FC) $(FFOPTS)

util_datetime.o : util_datetime.f90
	$(FC) $(FFOPTS)
