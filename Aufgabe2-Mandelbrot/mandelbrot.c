#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <mpi.h>

typedef unsigned char Tuchar;

typedef struct {
	int start, end;
} TSlice;

typedef struct {
	int iWidth, iHeight, iSlices;
	unsigned int uiMaxIterations;
	double dReMin, dReMax, dImMin, dImMax, dReFactor, dImFactor;
	char* FileName;
} TImageConfig;

Tuchar* allocateTucharP(int* size) {

	Tuchar* tmp = (Tuchar*)malloc(*size);

	if(tmp==0) {
		fprintf(stderr, "Out of memory!\n");
		exit(1);
	}

	memset(tmp,0,sizeof(tmp));

	return tmp;
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

Tuchar* getBmpFileHeader(int* filesize) {

	/* reserving some memory */
	int iMemory = 14;
	Tuchar* pHeader = allocateTucharP(&iMemory);

	/* setting the file header */
	pHeader[ 0] = 'B';
	pHeader[ 1] = 'M';
	pHeader[ 2] = (Tuchar)(*filesize    );
	pHeader[ 3] = (Tuchar)(*filesize>> 8);
	pHeader[ 4] = (Tuchar)(*filesize>>16);
	pHeader[ 5] = (Tuchar)(*filesize>>24);
	pHeader[ 6] = 0;
	pHeader[ 7] = 0;
	pHeader[ 8] = 0;
	pHeader[ 9] = 0;
	pHeader[10] = 54;
	pHeader[11] = 0;
	pHeader[12] = 0;
	pHeader[13] = 0;

	return pHeader;
}

Tuchar* getBmpInfoHeader(int* iImageWidth, int* iImageHeight) {

	/* reserving some memory */
	int iMemory = 40;
	Tuchar* pHeader = allocateTucharP(&iMemory);

        /* setting the info header */
	pHeader[ 0] = 40;
	pHeader[ 1] = 0;
	pHeader[ 2] = 0;
	pHeader[ 3] = 0;
	pHeader[ 4] = (Tuchar)( *iImageWidth    );
	pHeader[ 5] = (Tuchar)( *iImageWidth>> 8);
	pHeader[ 6] = (Tuchar)( *iImageWidth>>16);
	pHeader[ 7] = (Tuchar)( *iImageWidth>>24);
	pHeader[ 8] = (Tuchar)( *iImageHeight    );
	pHeader[ 9] = (Tuchar)( *iImageHeight>> 8);
	pHeader[10] = (Tuchar)( *iImageHeight>>16);
	pHeader[11] = (Tuchar)( *iImageHeight>>24);
	pHeader[12] = 1;
	pHeader[13] = 0;
	pHeader[14] = 24;
	pHeader[15] = 0;

	return pHeader;
}

FILE* getBmpFileHandler(const char* filename) {

	FILE* tmp = fopen(filename,"wb");

	if(tmp==0) {
		fprintf(stderr, "Can not write to file '%s'!\n", filename);
		exit(1);
	}

	return tmp;
}

void writeBmp(int* iImageWidth, int* iImageHeight, const char* FileName, const Tuchar* pImage) {

	int iLoop;

	/* computing the filesize */
 	int iFileSize = 54 + 3 * (*iImageWidth) * (*iImageHeight);

	/* file header */
	Tuchar* pFileHeader = getBmpFileHeader(&iFileSize);
	Tuchar* pInfoHeader = getBmpInfoHeader(iImageWidth, iImageHeight);
	Tuchar bmppad[3] = {0,0,0};

	/* open the file */
	FILE* aFile = getBmpFileHandler(FileName);
	/* write the file header */
	fwrite(pFileHeader, 1, 14, aFile);
	/* write the info header */
	fwrite(pInfoHeader, 1, 40, aFile);

	/* write the image contents */
	for(iLoop = 0; iLoop < *iImageHeight; iLoop++) {
	    fwrite(pImage + (*iImageWidth * (*iImageHeight - iLoop - 1) * 3), 3, *iImageWidth, aFile);
	    fwrite(bmppad, 1, (4 - (*iImageWidth * 3) % 4) % 4, aFile);
	}
	
	/* close the file */
	fclose(aFile);
}

void iterateAndStoreAPoint(TImageConfig* image, int* x, int* y, double* dZre, double* dZim, double* dCre, 
			   double* dCim, Tuchar* img) {

	int bIsInside = 1;
	double dColor;
	unsigned uN;

	/* iterate */
	for(uN = 0; uN < image->uiMaxIterations; ++uN) {
		double dZre2 = (*dZre) * (*dZre), dZim2 = (*dZim) * (*dZim);
		if(dZre2 + dZim2 > 4) {
			bIsInside = 0;
			dColor = uN;
			break;
		}
		*dZim = 2 * (*dZre) * (*dZim) + (*dCim);
		*dZre = dZre2 - dZim2 + (*dCre);
	}

	if(bIsInside==1) { 
		img[(*x + *y * image->iWidth) * 3 + 2] = (Tuchar)(0); /* r */
		img[(*x + *y * image->iWidth) * 3 + 1] = (Tuchar)(0); /* g */
		img[(*x + *y * image->iWidth) * 3 + 0] = (Tuchar)(0); /* b */
	} else {
		img[(*x + *y * image->iWidth) * 3 + 2] = (Tuchar)(dColor / image->uiMaxIterations * 255);
		img[(*x + *y * image->iWidth) * 3 + 1] = (Tuchar)(dColor / image->uiMaxIterations * 255 / 2);
		img[(*x + *y * image->iWidth) * 3 + 0] = (Tuchar)(dColor / image->uiMaxIterations * 255 / 2);
	}
}

void computeSlice(TImageConfig* image, TSlice* mySlice, Tuchar* img) {
	int y;
	for(y = mySlice->start; y < mySlice->end; ++y) {
		double dCim = image->dImMax - y * image->dImFactor;
		int x;
	        for(x=0; x < image->iWidth; ++x) {
			double dCre = image->dReMin + x * image->dReFactor;
			double dZre = dCre, dZim = dCim;
			/* iterate and store a points value */
			iterateAndStoreAPoint(image, &x, &y, &dZre, &dZim, &dCre, &dCim, img);
		}
	}
}

void initializeSlices(int* iHeight, int* iSlices, TSlice* sliceDimensions) {

	double dAverage = (double) *iHeight / *iSlices;
	int iRest = *iHeight - (*iSlices * (int)dAverage);

	/* computing the slice dimensions */
	int iSlice, iSliceStart = 0;
	for (iSlice = 0; iSlice < *iSlices; iSlice++) {
		int iSliceEnd = iSliceStart + dAverage;
		if(iSlice == *iSlices - 1) {
			iSliceEnd += iRest;
		}
		TSlice mySlice = {iSliceStart, iSliceEnd};
		sliceDimensions[iSlice] = mySlice;
		iSliceStart = iSliceEnd;
	}
}

void slaveReceive(int* rank, TImageConfig* image, TSlice* sliceDimensions, Tuchar* buffer) {

	int iSlice;
	MPI_Status status;

	for (;;) {
		/* receive a new slice to calculate */
		MPI_Recv(&iSlice, 1, MPI_INT, 0, 325, MPI_COMM_WORLD, &status);
		printf("Receiving slice '%d' in process '%d'\n", iSlice, *rank);

		/* suicide signal */
		if (iSlice < 0) {
			MPI_Finalize();
			exit(0);
		}

		/* calculate requested slice */
		printf("slice '%d' start-end: %d-%d\n", iSlice, sliceDimensions[iSlice].start, sliceDimensions[iSlice].end);
		computeSlice(image, &(sliceDimensions[iSlice]), buffer);

		/* send results back to master */
		MPI_Ssend(buffer, 3 * image->iWidth * image->iHeight, MPI_CHAR, 0, 327, MPI_COMM_WORLD);
	}
}

void slave(int* rank, TImageConfig* image) {
	/* reserving the transport buffer */
	Tuchar* buffer = allocateImageBuffer(image->iWidth, image->iHeight);

	/* reserving a buffer for the slice dimensions */
	TSlice sliceDimensions[image->iSlices];

	/* initialize the slices */
	initializeSlices(&(image->iHeight), &(image->iSlices), sliceDimensions);

	/* processing the slices */
	slaveReceive(rank, image, sliceDimensions, buffer);
}

int initializeProcesses(int* iSize, int* iSlices, int* iSliceCache) {
	int iProcess, iSlice=0;
	for (iProcess = 1; iProcess < *iSize; iProcess++) {
		/* generate Mandelbrot set in each process */
		if(iSlice<*iSlices) {
			printf("Sending slice '%d' to process '%d'\n", iSlice, iProcess);
			MPI_Ssend(&iSlice, 1, MPI_INT, iProcess, 325, MPI_COMM_WORLD);
			iSliceCache[iProcess] = iSlice;
			iSlice++;
		}
	}
	return iSlice;
}

void processSlices(TImageConfig* image, int* iSlice, TSlice* sliceDimensions, 
			int* iSliceCache, Tuchar* buffer, Tuchar* out) {

	int r_slice_n;
	for(r_slice_n = 0; r_slice_n<(*image).iSlices; r_slice_n++) {
		
		int source;
		MPI_Status status;

		MPI_Recv(buffer, 3*((*image).iWidth)*((*image).iHeight), 
			 MPI_CHAR, MPI_ANY_SOURCE, 327, MPI_COMM_WORLD, &status);

		/* source of the slice */
		source = status.MPI_SOURCE;
		/* assemble the image */
		printf("received computation for slice '%d' by process '%d'\nso char '%d' to '%d' are copied to out\n", 
			iSliceCache[source], source,
			image->iWidth * sliceDimensions[iSliceCache[source]].start * 3, 
			image->iWidth * sliceDimensions[iSliceCache[source]].end * 3);

		int x;
		for(x = image->iWidth * sliceDimensions[iSliceCache[source]].start * 3; 
			x < image->iWidth * sliceDimensions[iSliceCache[source]].end * 3; x++) {
			out[x] = buffer[x];
		}	

		/* compute another slice */
		if (*iSlice < image->iSlices) {
			MPI_Ssend(iSlice, 1, MPI_INT, source, 325, MPI_COMM_WORLD);
			iSliceCache[source] = *iSlice;
			(*iSlice)++;
		} else {
			int kill = -1;
			MPI_Ssend(&kill, 1, MPI_INT, source, 325, MPI_COMM_WORLD);
		}
	}
}

void master(int* iSize, TImageConfig* image) {
	
	/* reserving a buffer for the slice dimensions */
	TSlice sliceDimensions[image->iSlices];

	/* initialize the slices */
	initializeSlices(&(image->iHeight), &(image->iSlices), sliceDimensions);

	/* reserving the transport buffer */
	Tuchar* buffer = allocateImageBuffer(image->iWidth, image->iHeight);

	/* reserving the output buffer */
	Tuchar* out = allocateImageBuffer(image->iWidth, image->iHeight);

	/* initialize MPI: send a slice to each slave */
	int iSliceCache[*iSize];
	int iSlice = initializeProcesses(iSize, &(image->iSlices), iSliceCache);

	/* process remaining slices */
	processSlices(image, &iSlice, sliceDimensions, iSliceCache, buffer, out);

	/* finalize MPI */
	MPI_Finalize();

	/* writing the Mandelbrot set to a bitmap file */
	writeBmp(&(image->iWidth), &(image->iHeight), image->FileName, out);
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
		master(&iSize, &image);
	} else {
		slave(&rank, &image);	
	}
	return 0;
}
