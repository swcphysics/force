include ../Makefile.common





9.o : 9.f90 8.o 7.o
	$(FC) $(FFOPTS)


8.o : 8.f90 7.o
	$(FC) $(FFOPTS)


7.o : 7.f90 6.o
	$(FC) $(FFOPTS)

6.o : 6.f90 5.o $(DOBJ)/2.o $(DOBJ)/3.o $(DOBJ)/4.o 
	$(FC) $(FFOPTS)

5.o : 5.f90 $(DOBJ)/4.o $(DOBJ)/2.o
	$(FC) $(FFOPTS)









