smacofSSUI: Makefile smacofSSUI.c ../common/smacofPrint.c ../common/smacofUtils.c ../common/smacofCoreU.c ../include/smacof.h
	clang -o smacofSSUI smacofSSUI.c ../common/smacofPrint.c  ../common/smacofUtils.c ../common/smacofCoreU.c

rlib: 
	R CMD SHLIB smacofSSUI.c ../common/smacofPrint.c ../common/smacofUtils.c ../common/smacofCoreU.c

clean:
	rm -rf *.o ../common/*.o core a.out

pristine:
	rm -rf smacofSSUI smacofSSUI.so

format:
	clang-format -i *.c

