all:
	# To compile, type 'make compile'
compile:
	gcc -c kmcinput.c
	gcc -c kmclib.c
	gcc -o kmc kmcinput.o kmclib.o -lm -Wall
