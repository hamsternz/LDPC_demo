COPTS=-Wall -pedantic -g 
LIBS=-lcurses -lm

ldpc : ldpc.c
	gcc -o ldpc ldpc.c $(COPTS) $(LIBS)
