CC=gcc
CFLAGS=-Wall -O3 -DHAVE_INLINE
LDFLAGS=-lgsl -lgslcblas -lm
OBJFILES=submesh.o dos.o numstates.o weights.o sum.o

all: submesh.o dos.o numstates.o weights.o sum.o numstates_test.out

clean:
	rm *.o *.out

submesh.o: submesh.c submesh.h
	$(CC) $(CFLAGS) -c submesh.c

dos.o: dos.c dos.h
	$(CC) $(CFLAGS) -c dos.c

numstates.o: numstates.c numstates.h submesh.o
	$(CC) $(CFLAGS) -c numstates.c

weights.o: weights.c weights.h dos.o
	$(CC) $(CFLAGS) -c weights.c

sum.o: sum.c sum.h weights.o
	$(CC) $(CFLAGS) -c sum.c

numstates_test.out: numstates.o numstates_test.c
	$(CC) $(CFLAGS) -c numstates_test.c -o numstates_test.o
	$(CC) numstates_test.o -o numstates_test.out $(OBJFILES) $(LDFLAGS)
