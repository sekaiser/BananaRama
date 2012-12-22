# Minimal Makefile for the Apfelmann project
# (c) 2012 by R. Fruth, S. Kaiser, E. Kuhnt
#
# configuration of important variables
CC=mpicc
CFLAGS=-Wall -Wextra -ansi
OBJ=mandelbrot

compile: source/mandelbrot.c
	$(CC) $(CFLAGS) -o $(OBJ) source/mandelbrot.c

.PHONY: clean

clean:	
	rm -f $(OBJ)

