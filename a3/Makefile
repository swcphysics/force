include ../Makefile.common





LNKOPTS = -o $(DBIN)/$@ $(DBIN)/$@.o $(DOBJ)/*.o

400 : 400.o 300
	$(FC) $(LNKOPTS)

300 : 300.o 200
	$(FC) $(LNKOPTS)

200 : 200.o 100  
	$(FC) $(LNKOPTS)

100 : 100.o 
	$(FC) $(LNKOPTS)





400.o : 400.f90 300.o $(DOBJ)/8.o $(DOBJ)/9.o
	$(FC) $(ECOPTS)

300.o : 300.f90 200.o $(DOBJ)/7.o
	$(FC) $(ECOPTS)

200.o : 200.f90 100.o $(DOBJ)/8.o $(DOBJ)/9.o
	$(FC) $(ECOPTS)

100.o : 100.f90 $(DOBJ)/7.o
	$(FC) $(ECOPTS)





# 100,200,300,400的顺序能是任意的
# -o $(DBIN)/$@表示把主程序生成的可执行文件放在$(DBIN)中
# $(DBIN)/$@.o表示生成可执行文件需要$(DBIN)中的.o文件
# $(DOBJ)/*.o表示生成可执行文件需要$(DOBJ)中的.o文件
# 主程序生成的.o文件不能与module生成的.o文件放在同一个文件夹中



