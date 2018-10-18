include ../Makefile.common





4.o : 4.f90 3.o 2.o
	$(FC) $(FFOPTS)

3.o : 3.f90 2.o
	$(FC) $(FFOPTS)

2.o : 2.f90 1.o 
	$(FC) $(FFOPTS)

1.o : 1.f90 
	$(FC) $(FFOPTS)





# .f90文件后的.o文件的顺序能是任意的
#总结：
# N.o : N.f90 N-1.o  添加除N-1.o外得到N.o需要的.o文件
#	$(FC) $(FFOPTS)


# N-1.o : N-1.f90 N-2.o    添加除N-2.o外得到N-1.o需要的.o文件
#	$(FC) $(FFOPTS)

#...

# 4.o : 4.f90 3.o   添加除3.o外得到4.o需要的.o文件
#	$(FC) $(FFOPTS)

# 3.o : 3.f90 2.o   添加除2.o外得到3.o需要的.o文件
#	$(FC) $(FFOPTS)

# 2.o : 2.f90 1.o   添加除1.o外得到2.o需要的.o文件
#	$(FC) $(FFOPTS)

# 1.o : 1.f90    得到1.o不需要任何.o文件
#	$(FC) $(FFOPTS)
# 当在module j中有use i时i.o : i.f90 $(FC) $(FFOPTS)必须放在j.o : j.f90 $(FC) $(FFOPTS)后



