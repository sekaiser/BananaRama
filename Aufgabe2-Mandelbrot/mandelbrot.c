#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <mpi.h>
#include <getopt.h>

/* #define __mandelbrot__debug__ */

/* defines an alias for unsigned char */
typedef unsigned char Tuchar;

/* defines a slice */
typedef struct {
	int start, end;
} TSlice;

/* defines the image configuration */
typedef struct {
	int iWidth, iHeight, iSlices;
	unsigned int uiMaxIterations;
	double dReMin, dReMax, dImMin, dImMax, dReFactor, dImFactor;
	char* FileName;
} TImageConfig;

Tuchar* allocateTucharP(int* size) {

	/* allocate the buffer */
	Tuchar* tmp = (Tuchar*)malloc(*size);

	if(tmp==0) {
		fprintf(stderr, "Out of memory!\n");
		exit(1);
	}

	memset(tmp,0,sizeof(tmp));

	return tmp;
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

	/* create the file handler */
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

	int bIsInside = 1, iPosition;
	double dColor;
	unsigned uN;

	/* iterate the point */
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
	/* set the color config */
	iPosition = (*x + *y * image->iWidth) * 3;	
	if(bIsInside==1) { 
		img[iPosition] 	 = (Tuchar)(0); /* b */
		img[++iPosition] = (Tuchar)(0); /* g */
		img[++iPosition] = (Tuchar)(255); /* r */
	} else {
		/*dColor = dColor / image->uiMaxIterations * 255;*/
		img[iPosition]   = (Tuchar)(0);
		img[++iPosition] = (Tuchar)(0);
		img[++iPosition] = (Tuchar) dColor;
	}
}

void computeSlice(TImageConfig* image, TSlice* mySlice, Tuchar* img) {
	int y;
	/* loop over the pixels */
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

#ifdef __mandelbrot__debug__
void slaveReceiveSlices(int* rank, TImageConfig* image, TSlice* sliceDimensions, Tuchar* buffer) {
#else
void slaveReceiveSlices(TImageConfig* image, TSlice* sliceDimensions, Tuchar* buffer) {
#endif

	int iSlice;
	MPI_Status status;

	for (;;) {
		/* receive a new slice to calculate */
		MPI_Recv(&iSlice, 1, MPI_INT, 0, 325, MPI_COMM_WORLD, &status);
		#ifdef __mandelbrot__debug__
			printf("Receiving slice '%d' in process '%d'\n", iSlice, *rank);
		#endif
		/* suicide signal */
		if (iSlice < 0) {
			MPI_Finalize();
			exit(0);
		}

		/* calculate requested slice */
		#ifdef __mandelbrot__debug__
		printf("slice '%d' start-end: %d-%d\n", iSlice, sliceDimensions[iSlice].start, sliceDimensions[iSlice].end);
		#endif
		computeSlice(image, &(sliceDimensions[iSlice]), buffer);

		/* send results back to master */
		MPI_Ssend(buffer, 3 * image->iWidth * image->iHeight, MPI_CHAR, 0, 327, MPI_COMM_WORLD);
	}
}

void slaveReceiveConfig(TImageConfig* config) {
	
	/* creating the transport buffer */
	int iBuffer;
	unsigned int uiBuffer;
	double dBuffer;
	/* get iWidth */
	MPI_Recv(&iBuffer, 	1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	config->iWidth = iBuffer;
	/* get iHeight */
	MPI_Recv(&iBuffer, 	1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	config->iHeight = iBuffer;
	/* get uiMaxIterations */
	MPI_Recv(&uiBuffer, 	1, MPI_UNSIGNED, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	config->uiMaxIterations = uiBuffer;
	/* get iSlices */
	MPI_Recv(&iBuffer, 	1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	config->iSlices = iBuffer;	
	/* get dReMin */
	MPI_Recv(&dBuffer, 	1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	config->dReMin = dBuffer;
	/* get dReMax */	
	MPI_Recv(&dBuffer, 	1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	config->dReMax = dBuffer;
	/* get dImMin */
	MPI_Recv(&dBuffer, 	1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	config->dImMin = dBuffer;
	/* get dImMax */
	 MPI_Recv(&dBuffer, 	1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	config->dImMax = dBuffer;
	/* get dReFactor */
	MPI_Recv(&dBuffer, 	1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	config->dReFactor = dBuffer;
	/* get dImFactor */
	MPI_Recv(&dBuffer, 	1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	config->dImFactor = dBuffer;

}

#ifdef __mandelbrot__debug__
void slave(int* rank) {
#else
void slave() {
#endif
	/* get the config */
	TImageConfig image;
	slaveReceiveConfig(&image);
	
	/* reserving a buffer for the slice dimensions */
	TSlice sliceDimensions[image.iSlices];

	/* reserving the transport buffer */
	int iImgSize = 3 * image.iWidth * image.iHeight;	
	Tuchar* buffer = allocateTucharP(&iImgSize);

	/* initialize the slices */
	initializeSlices(&(image.iHeight), &(image.iSlices), sliceDimensions);

	/* processing the slices */
	#ifdef __mandelbrot__debug__
	slaveReceiveSlices(rank, &image, sliceDimensions, buffer);
	#else
	slaveReceiveSlices(&image, sliceDimensions, buffer);
	#endif
}

void configureProcesses(int* iSize, TImageConfig* image) {
	int iProcess;
	/* spread the config */
	for(iProcess = 1; iProcess < *iSize; iProcess++) {
		MPI_Ssend(&(image->iWidth), 	1, 	MPI_INT, 	iProcess, 	1, 	MPI_COMM_WORLD);
		/* get iHeight */
		MPI_Ssend(&(image->iHeight), 	1, 	MPI_INT, 	iProcess, 	1, 	MPI_COMM_WORLD);
		/* get uiMaxIterations */
		MPI_Ssend(&(image->uiMaxIterations), 1, MPI_UNSIGNED, 	iProcess, 	1, 	MPI_COMM_WORLD);
		/* get iSlices */
		MPI_Ssend(&(image->iSlices), 	1, 	MPI_INT, 	iProcess, 	1, 	MPI_COMM_WORLD);
		/* get dReMin */
		MPI_Ssend(&(image->dReMin), 	1, 	MPI_DOUBLE, 	iProcess, 	1, 	MPI_COMM_WORLD);
		/* get dReMax */
		MPI_Ssend(&(image->dReMax), 	1, 	MPI_DOUBLE, 	iProcess, 	1, 	MPI_COMM_WORLD);
		/* get dImMin */
		MPI_Ssend(&(image->dImMin), 	1, 	MPI_DOUBLE, 	iProcess, 	1, 	MPI_COMM_WORLD);
		/* get dImMax */
		MPI_Ssend(&(image->dImMax), 	1, 	MPI_DOUBLE, 	iProcess, 	1, 	MPI_COMM_WORLD);
		/* get dReFactor */
		MPI_Ssend(&(image->dReFactor), 	1, 	MPI_DOUBLE, 	iProcess, 	1, 	MPI_COMM_WORLD);
		/* get dImFactor */
		MPI_Ssend(&(image->dImFactor), 	1, 	MPI_DOUBLE, 	iProcess, 	1, 	MPI_COMM_WORLD);
	}
}

int initializeProcesses(int* iSize, int* iSlices, int* iSliceCache) {
	int iProcess, iSlice=0;
	for (iProcess = 1; iProcess < *iSize; iProcess++) {
		/* send first slice to each process */
		if(iSlice<*iSlices) {
			#ifdef __mandelbrot__debug__
			printf("Sending slice '%d' to process '%d'\n", iSlice, iProcess);
			#endif
			MPI_Ssend(&iSlice, 1, MPI_INT, iProcess, 325, MPI_COMM_WORLD);
			iSliceCache[iProcess] = iSlice;
			iSlice++;
		}
	}
	return iSlice;
}

void processSlices(TImageConfig* image, int* iSlice, TSlice* sliceDimensions, 
			int* iSliceCache, Tuchar* buffer, Tuchar* out) {

	int iRecvSlice;
	/* receive slices */
	for(iRecvSlice = 0; iRecvSlice < image->iSlices; iRecvSlice++) {
		int iSource;
		MPI_Status status;
		/* receive the slice */
		MPI_Recv(buffer, 3 * image->iWidth * image->iHeight, 
			 MPI_CHAR, MPI_ANY_SOURCE, 327, MPI_COMM_WORLD, &status);

		/* source of the slice */
		iSource = status.MPI_SOURCE;
		/* assemble the image */
		#ifdef __mandelbrot__debug__
		printf("received computation for slice '%d' by process '%d'\nso char '%d' to '%d' are copied to out\n", 
			iSliceCache[iSource], iSource,
			image->iWidth * sliceDimensions[iSliceCache[iSource]].start * 3, 
			image->iWidth * sliceDimensions[iSliceCache[iSource]].end * 3);
		#endif
		int x;
		/* cache the slice values */
		for(x = image->iWidth * sliceDimensions[iSliceCache[iSource]].start * 3; 
			x < image->iWidth * sliceDimensions[iSliceCache[iSource]].end * 3; x++) {
			out[x] = buffer[x];
		}	

		/* send another slice */
		if (*iSlice < image->iSlices) {
			MPI_Ssend(iSlice, 1, MPI_INT, iSource, 325, MPI_COMM_WORLD);
			iSliceCache[iSource] = *iSlice;
			(*iSlice)++;
		} else {
			/* send the kill signal */
			int kill = -1;
			MPI_Ssend(&kill, 1, MPI_INT, iSource, 325, MPI_COMM_WORLD);
		}
	}
}

void master(int* iSize, TImageConfig* image) {
	
	/* configure the processes */
	configureProcesses(iSize, image);	
	
	/* reserving a buffer for the slice dimensions */
	TSlice sliceDimensions[image->iSlices];

	/* initialize the slices */
	initializeSlices(&(image->iHeight), &(image->iSlices), sliceDimensions);

	int iImgSize = 3 * image->iWidth * image->iHeight;
	/* reserving the transport buffer */
	Tuchar* buffer = allocateTucharP(&iImgSize);

	/* reserving the output buffer */
	Tuchar* out = allocateTucharP(&iImgSize);

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

void initializeConfig(TImageConfig* image) {

	image->iWidth = 640;
	image->iHeight = 480;
	image->uiMaxIterations = 120;
	image->iSlices = 11;
	image->dReMin = -2.0;
	image->dReMax = 1.0;
	image->dImMin = -1.2;
	image->FileName = "mandelbrot.bmp";

}

void finalizeConfig(TImageConfig* image) {

	image->dImMax = image->dImMin + (image->dReMax - image->dReMin) * image->iHeight / image->iWidth;
	image->dReFactor = (image->dReMax - image->dReMin) / (image->iWidth - 1);
	image->dImFactor = (image->dImMax - image->dImMin) / (image->iHeight - 1);

}

void bestpictureConfig(TImageConfig* config) {

	config->iWidth = 800;
	config->iHeight= 600;
	config->uiMaxIterations = 255;
	config->dReMin = -0.3871904296875;
	config->dReMax = -0.3858232421874879;
	config->dImMin = 0.6238156738281248;
	config->FileName = "mandelbrot.bmp";

}

void dialog(TImageConfig* config) {
	printf("here would the dialog start\n");	
}

void parseCommandLineParameters(int argc, char** argv, TImageConfig* config) {

	int iIndex, bDialogFlag = 0, bBestPicFlag = 0;

	/* we will write directly to the config in order to prevent more memory allocations */ 
	while (1) {
		/* getopt_long stores the option index here. */
		int iOptionIndex = 0;
		/* options have to be declared here -- otherwise they will not work as expected */
		struct option oLongOptions[] = {
			/* These options set a flag. */
			{"dialog", 	no_argument,       	0, 	'd'},
			{"bestpic",	no_argument,	   	0, 	'p'},
			/* These options don't set a flag.
			   We distinguish them by their indices. */
			{"width",     	required_argument,      0, 	'w'},
			{"height",  	required_argument,      0, 	'h'},
			{"iterations",  required_argument, 	0, 	'n'},
			{"immin",  	required_argument, 	0, 	'a'},
			/* skipping one indice in order be able to improve later */
			{"remin",	required_argument,	0, 	'i'},
			{"remax",	required_argument,	0, 	'j'},
			{"filename",    required_argument, 	0,	'f'},
		       	{0, 		0, 			0,	0}
		};
		/* get the next option */
		iIndex = getopt_long (argc, argv, "dpw:h:n:a:i:j:f:", oLongOptions, &iOptionIndex);
		/* Detect the end of the options. */
		if (iIndex == -1)
			break;
		/* loop through the different option possibilities */
		switch (iIndex) {
             		case 'd':
				bDialogFlag = 1;
               			break;
			case 'p':
				bBestPicFlag = 1;
				break;
             		case 'w':
				config->iWidth = atoi(optarg);
				/* TODO: print error if <0 */
               			break;
             		case 'h':
				config->iHeight = atoi(optarg);
				/* TODO: print error if <0 */
               			break;
             		case 'n':
				config->uiMaxIterations = (unsigned int) atoi(optarg);
				/* TODO: print error if <0 */
				break;
			case 'a':
				config->dImMin = atof(optarg);
				break;
			case 'i':
				config->dReMin = atof(optarg);
				break;
			case 'j':
				config->dReMax = atof(optarg);
				break;
			case 'f':
				strcpy(config->FileName, optarg);
				/* TODO: print error if =="" */
				break;
             		case '?':
               			/* getopt_long already printed an default error message. */
               			break;
             		default:
				/* when does this get called? maybe if no parameter is set? */
				/* TODO: print usage information!? */
               			exit(1);
		}
	}

	if(bDialogFlag == 1 && bBestPicFlag == 1) {
		/* TODO: print error msg or do something */
	}

	if(bDialogFlag == 1){
		dialog(config);
	}
	
	if(bBestPicFlag == 1) {
		bestpictureConfig(config);
	}
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

	/* doing the work */
	if(rank==0) {
		/* master branch */
		/* initialize a default config */
		TImageConfig image;
		initializeConfig(&image);
		/* parse command line parameters and update the config */
		parseCommandLineParameters(argc, argv, &image);
		/* finalize config - automatic computation of dynamic 
   		   configuration values */
		finalizeConfig(&image);
		/* do master tasks */
		master(&iSize, &image);
	} else {
		/* slave branch */
		#ifdef __mandelbrot__debug__
		slave(&rank);
		#else
		slave();
		#endif
	}
	return 0;
}
