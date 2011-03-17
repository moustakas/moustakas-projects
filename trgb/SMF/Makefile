CC = gcc
LD = ld

mlike : 
	cc -c smf_lf.c /deep1/ioannis/lib/c/nrutil.o
	ld -G -o smf_lf.so smf_lf.o /deep1/ioannis/lib/c/nrutil.o

clean :
	$(RM) *.o smf_lf.o smf_lf.so
