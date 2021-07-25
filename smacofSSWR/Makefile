smacofSSWR: Makefile smacofSSWR.c ../common/smacofPrint.c ../common/smacofUtils.c ../common/smacofCoreW.c ../include/smacof.h
	clang -o smacofSSWR smacofSSWR.c ../common/smacofPrint.c ../common/smacofUtils.c ../common/smacofCoreW.c

rlib: 
	R CMD SHLIB smacofSSWR.c ../common/smacofPrint.c ../common/smacofUtils.c ../common/smacofCoreW.c

clean:
	rm -rf *.o ../common/*.o core a.out

pristine:
	rm -rf smacofSSWRE smacofSSWR.so

format:
	clang-format -i *.c

