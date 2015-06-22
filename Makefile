CC=gcc
CFLAGS=-Wall -O3 -DHAVE_INLINE
LDFLAGS=-lgsl -lgslcblas -lm
OBJFILES=dos.o numstates.o weights.o

all: dos.o numstates.o weights.o

clean:
	rm *.o *.out

dos.o: dos.c dos.h
	$(CC) $(CFLAGS) -c dos.c

numstates.o: numstates.c numstates.h
	$(CC) $(CFLAGS) -c numstates.c

weights.o: weights.c weights.h dos.o
	$(CC) $(CFLAGS) -c weights.c

#HTightBinding_test.out: HTightBinding.o bstrlib/bstrlib.o bstrlib/bstraux.o HTightBinding_test.c
#	$(CC) $(CFLAGS) -c HTightBinding_test.c -o HTightBinding_test.o
#	$(CC) HTightBinding_test.o -o HTightBinding_test.out $(OBJFILES) $(LDFLAGS)
