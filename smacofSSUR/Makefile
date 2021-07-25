smacofSSUR: Makefile smacofSSUR.c ../common/smacofPrint.c ../common/smacofUtils.c  ../common/smacofCoreU.c ../include/smacof.h
	clang -o smacofSSUR smacofSSUR.c ../common/smacofPrint.c ../common/smacofUtils.c ../common/smacofCoreU.c

rlib: 
	R CMD SHLIB smacofSSUR.c ../common/smacofPrint.c ../common/smacofUtils.c ../common/smacofCoreU.c

clean:
	rm -rf *.o ../common/*.o core a.out

pristine:
	rm -rf smacofSSUR smacofSSUR.so

format:
	clang-format -i *.c

