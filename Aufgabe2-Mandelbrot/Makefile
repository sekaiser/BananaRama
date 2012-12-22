# Minimal Makefile for the Apfelmann project
# by R. Fruth, S. Kaiser, E. Kuhnt
#
# configuration of important variables
CC=mpicc
CFLAGS=-Wall -Wextra -ansi
OBJ=mandelbrot

compile: mandelbrot.c 
	$(CC) $(CFLAGS) -o $(OBJ) mandelbrot.c

.PHONY: clean

clean:	
	rm -f $(OBJ)

