CC = gcc
LD = ld

bootsmf : 
	cc -c -KPIC -G bootsmf.c
	ld -G -o bootsmf.so bootsmf.o

clear :
	$(RM) *.o bootsmf bootsmf.so
