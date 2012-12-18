#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <mpi.h>

typedef unsigned char Tuchar;

typedef struct {
	int start, end;
} sliceT;

typedef struct {
	int iWidth, iHeight, iSlices;
	unsigned int uiMaxIterations;
	double dReMin, dReMax, dImMin, dImMax, dReFactor, dImFactor;
	char* FileName;
} TImageConfig;

Tuchar* allocateTucharP(int size) {

	Tuchar* tmp = (Tuchar*)malloc(size);
	if(tmp==0) {
		fprintf(stderr, "Out of memory!\n");
		exit(1);
	}
	memset(tmp,0,sizeof(tmp));
	return tmp;
}

Tuchar* bmp_get_fileheader(int filesize) {
	/* reserving some memory */
	Tuchar* header = allocateTucharP(14);
	/* setting the file header */
	header[ 0] = 'B';
	header[ 1] = 'M';
	header[ 2] = (Tuchar)(filesize    );
	header[ 3] = (Tuchar)(filesize>> 8);
	header[ 4] = (Tuchar)(filesize>>16);
	header[ 5] = (Tuchar)(filesize>>24);
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

Tuchar* bmp_get_infoheader(int iImageWidth, int iImageHeight) {
	Tuchar* header = allocateTucharP(40);
        /* setting the info header */
	header[ 0] = 40;
	header[ 1] = 0;
	header[ 2] = 0;
	header[ 3] = 0;
	header[ 4] = (Tuchar)( iImageWidth    );
	header[ 5] = (Tuchar)( iImageWidth>> 8);
	header[ 6] = (Tuchar)( iImageWidth>>16);
	header[ 7] = (Tuchar)( iImageWidth>>24);
	header[ 8] = (Tuchar)( iImageHeight    );
	header[ 9] = (Tuchar)( iImageHeight>> 8);
	header[10] = (Tuchar)( iImageHeight>>16);
	header[11] = (Tuchar)( iImageHeight>>24);
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

void bmp_write_output(int iImageWidth, int iImageHeight, const char* FileName, 
		 const Tuchar* img) {

	/* computing the filesize */
 	int filesize = 54 + 3*iImageWidth*iImageHeight;
	/* file header */
	Tuchar* fileheader = bmp_get_fileheader(filesize);
	Tuchar* infoheader = bmp_get_infoheader(iImageWidth, iImageHeight);
	Tuchar bmppad[3] = {0,0,0};

	/* open the file */
	FILE* f = bmp_get_filehandle(FileName);
	/* write the file header */
	fwrite(fileheader,1,14,f);
	/* write the info header */
	fwrite(infoheader,1,40,f);

	/* write the image contents */
	int i;
	for(i=0; i<iImageHeight; i++) {
	    fwrite(img+(iImageWidth*(iImageHeight-i-1)*3),3,iImageWidth,f);
	    fwrite(bmppad,1,(4-(iImageWidth*3)%4)%4,f);
	}
	
	/* close the file */
	fclose(f);
}

void point_iterate_and_store(TImageConfig* image, int x, int y, double Z_re, double Z_im, double c_re, double c_im, 
	     Tuchar* img) {
	int isInside = 1;
	double color;
	unsigned n;
	/* iterate */
	for(n=0; n<(*image).uiMaxIterations; ++n) {
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
		img[(x+y*(*image).iWidth)*3+2] = (Tuchar)(0); /* r */
		img[(x+y*(*image).iWidth)*3+1] = (Tuchar)(0); /* g */
		img[(x+y*(*image).iWidth)*3+0] = (Tuchar)(0); /* b */
	} else {
		img[(x+y*(*image).iWidth)*3+2] = (Tuchar)(color / (*image).uiMaxIterations * 255);
		img[(x+y*(*image).iWidth)*3+1] = (Tuchar)(color / (*image).uiMaxIterations * 255 / 2);
		img[(x+y*(*image).iWidth)*3+0] = (Tuchar)(color / (*image).uiMaxIterations * 255 / 2);
	}
}

void compute_slice(TImageConfig* image, sliceT mySlice, Tuchar* img) {
	int y;
	for(y=mySlice.start; y<mySlice.end; ++y) {
		double c_im = (*image).dImMax - y*(*image).dImFactor;
		int x;
	        for(x=0; x<(*image).iWidth; ++x) {
			double c_re = (*image).dReMin + x*(*image).dReFactor;
			double Z_re = c_re, Z_im = c_im;
			/* iterate and store a points value */
			point_iterate_and_store(image, x, y, Z_re, Z_im, c_re, c_im, img);
		}
	}
}

void initializeSlices(int iImageHeight, int iN_SLICES, sliceT* slice_dimensions) {
	double average = (double)iImageHeight / iN_SLICES;
	int rest = iImageHeight - (iN_SLICES*(int)average);
	/* computing the slice dimensions */
	int slice_n;
	int slice_start = 0;
	for (slice_n=0; slice_n<iN_SLICES; slice_n++) {
		int slice_end = slice_start + average;
		if(slice_n==iN_SLICES-1) {
			slice_end += rest;
		}
		sliceT mySlice = {slice_start, slice_end};
		slice_dimensions[slice_n] = mySlice;
		slice_start = slice_end;
	}
}

Tuchar* allocateImageBuffer(int iImageWidth, int iImageHeight) {
	Tuchar* img = (Tuchar*)malloc(3*iImageWidth*iImageHeight);
	if(img==0) {
		fprintf(stderr, "Out of memory!\n");
		exit(1);
	}
	memset(img,0,sizeof(img));
	return img;
}

void slave_recv(int rank, TImageConfig* image, sliceT* slice_dimensions, 
	Tuchar* transportbuffer) {

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
		MPI_Ssend(transportbuffer, 3*((*image).iWidth)*((*image).iHeight), MPI_CHAR, 0, 327, MPI_COMM_WORLD);
	}
}

void slave(int rank, TImageConfig* image) {
	/* reserving the transport buffer */
	Tuchar* transportbuffer = allocateImageBuffer((*image).iWidth, (*image).iHeight);
	/* reserving a buffer for the slice dimensions */
	sliceT slice_dimensions[(*image).iSlices];
	/* initialize the slices */
	initializeSlices((*image).iHeight, (*image).iSlices, slice_dimensions);
	/* processing the slices */
	slave_recv(rank, image, slice_dimensions, transportbuffer);
}

int mpi_slices_init(int mpiSize, int iN_SLICES, int* slice_cache) {
	int process, slice_n=0;
	for (process = 1; process < mpiSize; process++) {
		/* generate Mandelbrot set in each process */
		if(slice_n<iN_SLICES) {
			printf("Sending slice '%d' to process '%d'\n", slice_n, process);
			MPI_Ssend(&slice_n, 1, MPI_INT, process, 325, MPI_COMM_WORLD);
			slice_cache[process] = slice_n;
			slice_n++;
		}
	}
	return slice_n;
}

void mpi_slices_process(TImageConfig* image, int slice_n, sliceT* slice_dimensions, 
			int* slice_cache, Tuchar* transportbuffer, Tuchar* out) {

	int r_slice_n;
	for(r_slice_n = 0; r_slice_n<(*image).iSlices; r_slice_n++) {
		
		int source;
		MPI_Status status;

		MPI_Recv(transportbuffer, 3*((*image).iWidth)*((*image).iHeight), 
			 MPI_CHAR, MPI_ANY_SOURCE, 327, MPI_COMM_WORLD, &status);

		/* source of the slice */
		source = status.MPI_SOURCE;
		/* assemble the image */
		printf("received computation for slice '%d' by process '%d'\nso char '%d' to '%d' are copied to out\n", 
			slice_cache[source], source, (*image).iWidth*
			slice_dimensions[slice_cache[source]].start*3, 
			(*image).iWidth*slice_dimensions[slice_cache[source]].end*3);
		int x;
		for(	x=(*image).iWidth*slice_dimensions[slice_cache[source]].start*3; 
			x < (*image).iWidth*slice_dimensions[slice_cache[source]].end*3; x++) {
			out[x] = transportbuffer[x];
		}	
		/* compute another slice */
		if (slice_n < (*image).iSlices) {
			MPI_Ssend(&slice_n, 1, MPI_INT, source, 325, MPI_COMM_WORLD);
			slice_cache[source] = slice_n;
			slice_n++;
		} else {
			int kill = -1;
			MPI_Ssend(&kill, 1, MPI_INT, source, 325, MPI_COMM_WORLD);
		}
	}
}

void master(int mpiSize, TImageConfig* image) {
	
	/* reserving a buffer for the slice dimensions */
	sliceT slice_dimensions[(*image).iSlices];
	/* initialize the slices */
	initializeSlices((*image).iHeight, (*image).iSlices, slice_dimensions);
	/* reserving the transport buffer */
	Tuchar* transportbuffer = allocateImageBuffer((*image).iWidth, (*image).iHeight);
	/* reserving the output buffer */
	Tuchar* out = allocateImageBuffer((*image).iWidth, (*image).iHeight);
	/* initialize MPI: send a slice to each slave */
	int slice_cache[mpiSize];
	int slice_n = mpi_slices_init(mpiSize, (*image).iSlices, slice_cache);
	/* process remaining slices */
	mpi_slices_process(image, slice_n, slice_dimensions, slice_cache, transportbuffer,
			   out);
	/* finalize MPI */
	MPI_Finalize();
	/* writing the Mandelbrot set to a bitmap file */
	bmp_write_output((*image).iWidth, (*image).iHeight, (*image).FileName, out);
}

/*
 * the main routine
 */
int main(int argc, char **argv) {

	/* initialize MPI */
	int rank, iSize;
	MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &iSize);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	/* move this stuff to a dialog */
	TImageConfig image;
	image.iWidth = 640;
	image.iHeight = 480;
	image.uiMaxIterations = 120;
	image.iSlices = 11;
	image.dReMin = -2.0;
	image.dReMax = 1.0;
	image.dImMin = -1.2;
	image.dImMax = image.dImMin + (image.dReMax - image.dReMin) * image.iHeight / image.iWidth;
	image.dReFactor = (image.dReMax - image.dReMin) / (image.iWidth - 1);
	image.dImFactor = (image.dImMax - image.dImMin) / (image.iHeight - 1);
	image.FileName = "mandelbrot.bmp";

	/* doing the work */
	if(rank==0) {
		master(iSize, &image);
	} else {
		slave(rank, &image);	
	}
	return 0;
}
