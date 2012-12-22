# Minimal Makefile for the Apfelmann project
# (c) 2012 by R. Fruth, S. Kaiser, E. Kuhnt
#
# configuration of important variables
CC=mpicc
CFLAGS=-Wall -Wextra -ansi
OBJ=mandelbrot

compile: mandelbrot.c \
	 mandelbrot.h
	$(CC) $(CFLAGS) -o $(OBJ) mandelbrot.c

.PHONY: clean

clean:	
	rm -f $(OBJ)

