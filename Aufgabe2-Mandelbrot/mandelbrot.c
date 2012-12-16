#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <mpi.h>

typedef struct {
	int start, end;
} sliceT;

typedef struct {
	int Width, Height, N_SLICES;
	unsigned int MaxIterations;
	double MinRe, MaxRe, MinIm, MaxIm, Re_factor, Im_factor;
	char* FileName;
} ImageConfig;

unsigned char* allocateUCharP(int size) {

	unsigned char* tmp = (unsigned char *)malloc(size);
	if(tmp==0) {
		fprintf(stderr, "Out of memory!\n");
		exit(1);
	}
	memset(tmp,0,sizeof(tmp));
	return tmp;
}

unsigned char* bmp_get_fileheader(int filesize) {
	/* reserving some memory */
	unsigned char* header = allocateUCharP(14);
	/* setting the file header */
	header[ 0] = 'B';
	header[ 1] = 'M';
	header[ 2] = (unsigned char)(filesize    );
	header[ 3] = (unsigned char)(filesize>> 8);
	header[ 4] = (unsigned char)(filesize>>16);
	header[ 5] = (unsigned char)(filesize>>24);
	header[ 6] = 0;
	header[ 7] = 0;
	header[ 8] = 0;
	header[ 9] = 0;
	header[10] = 54;
	header[11] = 0;
	header[12] = 0;
	header[13] = 0;
	return header;
}

unsigned char* bmp_get_infoheader(int ImageWidth, int ImageHeight) {
	unsigned char* header = allocateUCharP(40);
        /* setting the info header */
	header[ 0] = 40;
	header[ 1] = 0;
	header[ 2] = 0;
	header[ 3] = 0;
	header[ 4] = (unsigned char)( ImageWidth    );
	header[ 5] = (unsigned char)( ImageWidth>> 8);
	header[ 6] = (unsigned char)( ImageWidth>>16);
	header[ 7] = (unsigned char)( ImageWidth>>24);
	header[ 8] = (unsigned char)( ImageHeight    );
	header[ 9] = (unsigned char)( ImageHeight>> 8);
	header[10] = (unsigned char)( ImageHeight>>16);
	header[11] = (unsigned char)( ImageHeight>>24);
	header[12] = 1;
	header[13] = 0;
	header[14] = 24;
	header[15] = 0;
	return header;
}

FILE* bmp_get_filehandle(const char* filename) {
	FILE* tmp = fopen(filename,"wb");
	if(tmp==0) {
		fprintf(stderr, "Can not write to file '%s'!\n", filename);
		exit(1);
	}
	return tmp;
}

void bmp_write_output(int ImageWidth, int ImageHeight, const char* FileName, 
		 const unsigned char* img) {

	/* computing the filesize */
 	int filesize = 54 + 3*ImageWidth*ImageHeight;
	/* file header */
	unsigned char* fileheader = bmp_get_fileheader(filesize);
	unsigned char* infoheader = bmp_get_infoheader(ImageWidth, ImageHeight);
	unsigned char bmppad[3] = {0,0,0};

	/* open the file */
	FILE* f = bmp_get_filehandle(FileName);
	/* write the file header */
	fwrite(fileheader,1,14,f);
	/* write the info header */
	fwrite(infoheader,1,40,f);

	/* write the image contents */
	int i;
	for(i=0; i<ImageHeight; i++) {
	    fwrite(img+(ImageWidth*(ImageHeight-i-1)*3),3,ImageWidth,f);
	    fwrite(bmppad,1,(4-(ImageWidth*3)%4)%4,f);
	}
	
	/* close the file */
	fclose(f);
}

void compute_slice(ImageConfig* image, sliceT mySlice, unsigned char* img) {

	unsigned y;
	for(y=mySlice.start; y<mySlice.end; ++y) {
		double c_im = (*image).MaxIm - y*(*image).Im_factor;
		unsigned x;
	        for(x=0; x<(*image).Width; ++x) {
			double c_re = (*image).MinRe + x*(*image).Re_factor;
			double Z_re = c_re, Z_im = c_im;
			int isInside = 1;
			double color;
			unsigned n;
			/* iterate */
			for(n=0; n<(*image).MaxIterations; ++n) {
				double Z_re2 = Z_re*Z_re, Z_im2 = Z_im*Z_im;
				if(Z_re2 + Z_im2 > 4) {
					isInside = 0;
					color = n;
					break;
				}
				Z_im = 2*Z_re*Z_im + c_im;
				Z_re = Z_re2 - Z_im2 + c_re;
			}
			if(isInside==1) { 
				// TODO: check if rgb > 255
				img[(x+y*(*image).Width)*3+2] = (unsigned char)(0); // r
				img[(x+y*(*image).Width)*3+1] = (unsigned char)(0); // g
				img[(x+y*(*image).Width)*3+0] = (unsigned char)(0); // b
			} else {
				img[(x+y*(*image).Width)*3+2] = (unsigned char)(color / (*image).MaxIterations * 255);
				img[(x+y*(*image).Width)*3+1] = (unsigned char)(color / (*image).MaxIterations * 255 / 2);
				img[(x+y*(*image).Width)*3+0] = (unsigned char)(color / (*image).MaxIterations * 255 / 2);
			}
		}
	}
}

void initializeSlices(int ImageHeight, int N_SLICES, sliceT* slice_dimensions) {
	double average = (double)ImageHeight / N_SLICES;
	int rest = ImageHeight - (N_SLICES*(int)average);
	/* computing the slice dimensions */
	int slice_n;
	int slice_start = 0;
	for (slice_n=0; slice_n<N_SLICES; slice_n++) {
		int slice_end = slice_start + average;
		if(slice_n==N_SLICES-1) {
			slice_end += rest;
		}
		sliceT mySlice = {slice_start, slice_end};
		slice_dimensions[slice_n] = mySlice;
		slice_start = slice_end;
	}
}

unsigned char* allocateImageBuffer(int ImageWidth, int ImageHeight) {
	unsigned char* img = (unsigned char *)malloc(3*ImageWidth*ImageHeight);
	if(img==0) {
		fprintf(stderr, "Out of memory!\n");
		exit(1);
	}
	memset(img,0,sizeof(img));
	return img;
}

void slave_recv(int rank, ImageConfig* image, sliceT* slice_dimensions, unsigned char* transportbuffer) {
	int slice_n;
	MPI_Status status;
	for (;;) {
		/* receive a new slice to calculate */
		MPI_Recv(&slice_n, 1, MPI_INT, 0, 325, MPI_COMM_WORLD, &status);
		printf("Receiving slice '%d' in process '%d'\n", slice_n, rank);
		/* suicide signal */
		if (slice_n < 0) {
			MPI_Finalize();
			exit(0);
		}
		/* calculate requested slice */
		printf("slice '%d' start-end: %d-%d\n", slice_n, slice_dimensions[slice_n].start, slice_dimensions[slice_n].end);
		compute_slice(image, slice_dimensions[slice_n], transportbuffer);
		/* send results back to master */
		MPI_Ssend(transportbuffer, 3*((*image).Width)*((*image).Height), MPI_CHAR, 0, 327, MPI_COMM_WORLD);
	}
}

void slave(int rank, ImageConfig* image) {
	/* reserving the transport buffer */
	unsigned char* transportbuffer = allocateImageBuffer((*image).Width, (*image).Height);
	/* reserving a buffer for the slice dimensions */
	sliceT slice_dimensions[(*image).N_SLICES];
	/* initialize the slices */
	initializeSlices((*image).Height, (*image).N_SLICES, slice_dimensions);
	/* processing the slices */
	slave_recv(rank, image, slice_dimensions, transportbuffer);
}

void master(int mpiSize, ImageConfig* image) {
	
	/* reserving a buffer for the slice dimensions */
	sliceT slice_dimensions[(*image).N_SLICES];
	/* initialize the slices */
	initializeSlices((*image).Height, (*image).N_SLICES, slice_dimensions);
	/* reserving the image buffer */
	unsigned char* transportbuffer = allocateImageBuffer((*image).Width, (*image).Height);
	/* reserving the output buffer */
	unsigned char* out = allocateImageBuffer((*image).Width, (*image).Height);

	/* 
	    initialize MPI
	    send a slice to each slave
	 */

	int slice_cache[mpiSize];
	int slice_n = 0;
	int process;
	for (process = 1; process < mpiSize; process++) {
		/* generate Mandelbrot set in each process */
		if(slice_n<(*image).N_SLICES) {
			printf("Sending slice '%d' to process '%d'\n", slice_n, process);
			MPI_Ssend(&slice_n, 1, MPI_INT, process, 325, MPI_COMM_WORLD);
			slice_cache[process] = slice_n;
			slice_n++;
		}
	}

	int r_slice_n;
	for(r_slice_n = 0; r_slice_n<(*image).N_SLICES; r_slice_n++) {
		
		int source;
		MPI_Status status;

		MPI_Recv(transportbuffer, 3*((*image).Width)*((*image).Height), 
			 MPI_CHAR, MPI_ANY_SOURCE, 327, MPI_COMM_WORLD, &status);

		/* source of the slice */
		source = status.MPI_SOURCE;
		/* assemble the image */
		printf("received computation for slice '%d' by process '%d'\nso char '%d' to '%d' are copied to out\n", 
			slice_cache[source], source, (*image).Width*
			slice_dimensions[slice_cache[source]].start*3, 
			(*image).Width*slice_dimensions[slice_cache[source]].end*3);
		int x;
		for(	x=(*image).Width*slice_dimensions[slice_cache[source]].start*3; 
			x < (*image).Width*slice_dimensions[slice_cache[source]].end*3; x++) {
			out[x] = transportbuffer[x];
		}	
		/* compute another slice */
		if (slice_n < (*image).N_SLICES) {
			MPI_Ssend(&slice_n, 1, MPI_INT, source, 325, MPI_COMM_WORLD);
			slice_cache[source] = slice_n;
			slice_n++;
		} else {
			int kill = -1;
			MPI_Ssend(&kill, 1, MPI_INT, source, 325, MPI_COMM_WORLD);
		}
	}

	MPI_Finalize();

	/* writing the Mandelbrot set to a bitmap file */
	bmp_write_output((*image).Width, (*image).Height, (*image).FileName, out);
}

/*
 * the main routine
 */
int main(int argc, char **argv) {

	/* initialize MPI */
	int rank, mpiSize;
	MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	/* move this stuff to a dialog */
	ImageConfig image;
	image.Width = 640;
	image.Height = 480;
	image.MaxIterations = 120;
	image.N_SLICES = 11;
	image.MinRe = -2.0;
	image.MaxRe = 1.0;
	image.MinIm = -1.2;
	image.MaxIm = image.MinIm+(image.MaxRe-image.MinRe)*image.Height/image.Width;
	image.Re_factor = (image.MaxRe-image.MinRe)/(image.Width-1);
	image.Im_factor = (image.MaxIm-image.MinIm)/(image.Height-1);
	image.FileName = "mandelbrot.bmp";

	/* doing the work */
	if(rank==0) {
		master(mpiSize, &image);
	} else {
		slave(rank, &image);	
	}
}
