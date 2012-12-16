#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <mpi.h>

typedef struct {
	int start, end;
} sliceT;

void bmp_write_output(int ImageWidth, int ImageHeight, const char* FileName, 
		 const unsigned char* img) {

	/* computing the filesize */
 	int filesize = 54 + 3*ImageWidth*ImageHeight;

	/* writing the file header */
	unsigned char bmpfileheader[14] = {'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0};
	unsigned char bmpinfoheader[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0};
	unsigned char bmppad[3] = {0,0,0};
	bmpfileheader[ 2] = (unsigned char)(filesize    );
	bmpfileheader[ 3] = (unsigned char)(filesize>> 8);
	bmpfileheader[ 4] = (unsigned char)(filesize>>16);
	bmpfileheader[ 5] = (unsigned char)(filesize>>24);

	/* writing the info header*/
	bmpinfoheader[ 4] = (unsigned char)( ImageWidth    );
	bmpinfoheader[ 5] = (unsigned char)( ImageWidth>> 8);
	bmpinfoheader[ 6] = (unsigned char)( ImageWidth>>16);
	bmpinfoheader[ 7] = (unsigned char)( ImageWidth>>24);
	bmpinfoheader[ 8] = (unsigned char)( ImageHeight    );
	bmpinfoheader[ 9] = (unsigned char)( ImageHeight>> 8);
	bmpinfoheader[10] = (unsigned char)( ImageHeight>>16);
	bmpinfoheader[11] = (unsigned char)( ImageHeight>>24);

	/* open the file */
	FILE* f;
	f = fopen(FileName,"wb");

	/* write the file header */
	fwrite(bmpfileheader,1,14,f);

	/* write the info header */
	fwrite(bmpinfoheader,1,40,f);

	/* write the image contents */
	int i;
	for(i=0; i<ImageHeight; i++) {
	    fwrite(img+(ImageWidth*(ImageHeight-i-1)*3),3,ImageWidth,f);
	    fwrite(bmppad,1,(4-(ImageWidth*3)%4)%4,f);
	}
	
	/* close the file */
	fclose(f);
}

void compute_slice(int MaxIterations, int ImageWidth, double MinRe, double MaxIm,
		   double Im_factor, double Re_factor, sliceT mySlice, unsigned char* img) {

	unsigned y;
	for(y=mySlice.start; y<mySlice.end; ++y) {
		double c_im = MaxIm - y*Im_factor;
		unsigned x;
	        for(x=0; x<ImageWidth; ++x) {
			double c_re = MinRe + x*Re_factor;
			double Z_re = c_re, Z_im = c_im;
			int isInside = 1;
			double color;
			unsigned n;
			/* iterate */
			for(n=0; n<MaxIterations; ++n) {
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
				img[(x+y*ImageWidth)*3+2] = (unsigned char)(0); // r
				img[(x+y*ImageWidth)*3+1] = (unsigned char)(0); // g
				img[(x+y*ImageWidth)*3+0] = (unsigned char)(0); // b
			} else {
				img[(x+y*ImageWidth)*3+2] = (unsigned char)(color / MaxIterations * 255);
				img[(x+y*ImageWidth)*3+1] = (unsigned char)(color / MaxIterations * 255 / 2);
				img[(x+y*ImageWidth)*3+0] = (unsigned char)(color / MaxIterations * 255 / 2);
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

unsigned char* allocateImageBuffer(int ImageWidth, int ImageHeight, unsigned char* img) {
	if(img) {
		free(img);
	}
	img = (unsigned char *)malloc(3*ImageWidth*ImageHeight);
	memset(img,0,sizeof(img));
	return img;
}

int main(int argc, char **argv) {

	/* initialize MPI */
	int rank, mpiSize;
	MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	/* TODO: move this stuff to a dialog */
	int ImageWidth = 640;
	int ImageHeight = 480;
	unsigned int MaxIterations = 120;
	int N_SLICES = 11;
	double MinRe = -2.0;
	double MaxRe = 1.0;
	double MinIm = -1.2;
	char* FileName = "mandelbrot.bmp";
	double MaxIm = MinIm+(MaxRe-MinRe)*ImageHeight/ImageWidth;
	double Re_factor = (MaxRe-MinRe)/(ImageWidth-1);
	double Im_factor = (MaxIm-MinIm)/(ImageHeight-1);

	/* reserving a buffer for the slice dimensions */
	sliceT slice_dimensions[N_SLICES];
	/* initialize the slices */
	initializeSlices(ImageHeight, N_SLICES, slice_dimensions);

	if(rank==0) {
		/* reserving the image buffer */
		unsigned char* transportbuffer = NULL;
		transportbuffer = allocateImageBuffer(ImageWidth, ImageHeight, transportbuffer);
		/* reserving the output buffer */
		unsigned char* out = NULL;
		out = allocateImageBuffer(ImageWidth,ImageHeight, out);
		/* 
		    initialize
		    send a slice to each slave
		 */

		int slice_distribution[mpiSize];
		int slice_n = 0;
		int process;
		for (process = 1; process < mpiSize; process++) {
			/* generate Mandelbrot set in each process */
			if(slice_n<N_SLICES) {
				printf("Sending slice '%d' to process '%d'\n", slice_n, process);
				MPI_Ssend(&slice_n, 1, MPI_INT, process, 325, MPI_COMM_WORLD);
				slice_distribution[process] = slice_n;
				slice_n++;
			}
		}

		int r_slice_n;
		for(r_slice_n = 0; r_slice_n<N_SLICES; r_slice_n++) {
		
			int source;
			MPI_Status status;

			MPI_Recv(transportbuffer, 3*ImageWidth*ImageHeight, MPI_CHAR, MPI_ANY_SOURCE, 327, MPI_COMM_WORLD, &status);

			/* source of the slice */
			source = status.MPI_SOURCE;
			/* assemble the image */
			printf("received computation for slice '%d' by process '%d'\nso char '%d' to '%d' are copied to out\n", 
				slice_distribution[source], source, ImageWidth*slice_dimensions[slice_distribution[source]].start*3, 
				ImageWidth*slice_dimensions[slice_distribution[source]].end*3);
			int x;
			for(	x=ImageWidth*slice_dimensions[slice_distribution[source]].start*3; 
				x < ImageWidth*slice_dimensions[slice_distribution[source]].end*3; x++) {
				out[x] = transportbuffer[x];
			}	
			/* compute another slice */
			if (slice_n < N_SLICES) {
				MPI_Ssend(&slice_n, 1, MPI_INT, source, 325, MPI_COMM_WORLD);
				slice_distribution[source] = slice_n;
				slice_n++;
			} else {
				int kill = -1;
				MPI_Ssend(&kill, 1, MPI_INT, source, 325, MPI_COMM_WORLD);
			}
		}

		MPI_Finalize();

		/* writing the Mandelbrot set to a bitmap file */
		bmp_write_output(ImageWidth, ImageHeight, FileName, out);
	} else {
		int slice_n;
		MPI_Status status;
		/* reserving the transport buffer */
		unsigned char* transportbuffer = NULL;
		transportbuffer = allocateImageBuffer(ImageWidth, ImageHeight, transportbuffer);
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
			/* TODO: remove! reserving the image buffer */
			compute_slice(MaxIterations, ImageWidth, MinRe, MaxIm, Im_factor, Re_factor, 
				      slice_dimensions[slice_n], transportbuffer);
			/* send results back to master */
			MPI_Ssend(transportbuffer, 3*ImageWidth*ImageHeight, MPI_CHAR, 0, 327, MPI_COMM_WORLD);
		}	
	}
}
