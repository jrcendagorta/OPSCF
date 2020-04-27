CFLAGS=-Wall -g
CC = g++

all: clean OPSCF

clean:
	rm -f OPSCF *.o

OPSCF:	OPSCF.o Integrators.o Initial.o Force.o
	$(CC) $(CFLAGS) -o OPSCF OPSCF.o Integrators.o Initial.o Force.o

OPSCF.o: System.h Initial.h Integrators.h Force.h

Initial.o: Initial.h System.h

Force.o: Force.h System.h Initial.h

Integrators.o: Integrators.h System.h
