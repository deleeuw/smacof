smacofSSWI: Makefile smacofSSWI.c ../common/smacofPrint.c ../common/smacofUtils.c ../common/smacofCoreW.c ../include/smacof.h
	clang -o smacofSSWI smacofSSWI.c ../common/smacofPrint.c ../common/smacofUtils.c ../common/smacofCoreW.c

rlib: 
	R CMD SHLIB smacofSSWI.c ../common/smacofPrint.c ../common/smacofUtils.c ../common/smacofCoreW.c

clean:
	rm -rf *.o ../common/*.o core a.out

pristine:
	rm -rf smacofSSWI smacofSSWI.so

format:
	clang-format -i *.c

