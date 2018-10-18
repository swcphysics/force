include ../Makefile.common






LNKOPTS = -o $(DBIN)/$@ $(DBIN)/$@.o $(DOBJ)/*.o

800 : 800.o 700
	$(FC) $(LNKOPTS)

700 : 700.o 600
	$(FC) $(LNKOPTS)

600 : 600.o 500
	$(FC) $(LNKOPTS)

500 : 500.o  
	$(FC) $(LNKOPTS)















800.o : 800.f90 700.o $(DOBJ)/8.o $(DOBJ)/9.o
	$(FC) $(ECOPTS)

700.o : 700.f90 600.o $(DOBJ)/7.o
	$(FC) $(ECOPTS)

600.o : 600.f90 500.o $(DOBJ)/8.o $(DOBJ)/9.o
	$(FC) $(ECOPTS)

500.o : 500.f90 $(DOBJ)/7.o
	$(FC) $(ECOPTS)






