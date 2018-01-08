install: diva

clean:
	rm *.o

OBJECTS = diva.o functions.o help.o

HEADERS = diva.h

CC	= cc
CFLAGS  = -c -O3
LFLAGS	= -lm 

functions.o : functions.c $(HEADERS) diva.h
	$(CC) $(CFLAGS) $*.c

diva.o : diva.c $(HEADERS)
	$(CC) $(CFLAGS) $*.c


help.o : help.c $(HEADERS)
	$(CC) $(CFLAGS) $*.c



diva : $(OBJECTS) diva.o functions.o help.o
	$(CC) -o $@ $(OBJECTS) $(LFLAGS)


