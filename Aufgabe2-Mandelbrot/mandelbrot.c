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
	int iWidth, iHeight, iSlices, iRed, iGreen, iBlue, iRedBg, iGreenBg, iBlueBg;
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
		img[iPosition] 	 = (Tuchar)(image->iBlueBg);
		img[++iPosition] = (Tuchar)(image->iGreenBg);
		img[++iPosition] = (Tuchar)(image->iRedBg);
	} else {
		dColor = dColor / image->uiMaxIterations;
		img[iPosition]   = (Tuchar) (dColor * image->iBlue);
		img[++iPosition] = (Tuchar) (dColor * image->iGreen);
		img[++iPosition] = (Tuchar) (dColor * image->iRed);
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

	/* compute average slice height and the rest (will be added to the 
	   last average slice) */
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
	
	/* get iWidth */
	MPI_Recv(&(config->iWidth), 	1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	/* get iHeight */
	MPI_Recv(&(config->iHeight), 	1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	/* get uiMaxIterations */
	MPI_Recv(&(config->uiMaxIterations), 	1, MPI_UNSIGNED, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	/* get iSlices */
	MPI_Recv(&(config->iSlices), 	1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	/* get iRed */
	MPI_Recv(&(config->iRed), 	1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	/* get iGreen */
	MPI_Recv(&(config->iGreen), 	1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	/* get iBlue */
	MPI_Recv(&(config->iBlue), 	1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	/* get iRedBg */
	MPI_Recv(&(config->iRedBg), 	1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	/* get iGreenBg */
	MPI_Recv(&(config->iGreenBg), 	1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	/* get iBlueBg */
	MPI_Recv(&(config->iBlueBg), 	1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	/* get dReMin */
	MPI_Recv(&(config->dReMin), 	1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	/* get dReMax */	
	MPI_Recv(&(config->dReMax), 	1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	/* get dImMin */
	MPI_Recv(&(config->dImMin), 	1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	/* get dImMax */
	 MPI_Recv(&(config->dImMax), 	1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	/* get dReFactor */
	MPI_Recv(&(config->dReFactor), 	1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	/* get dImFactor */
	MPI_Recv(&(config->dImFactor), 	1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

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
		/* get iRed */
		MPI_Ssend(&(image->iRed), 	1, 	MPI_INT, 	iProcess, 	1, 	MPI_COMM_WORLD);
		/* get iGreen */
		MPI_Ssend(&(image->iGreen), 	1, 	MPI_INT, 	iProcess, 	1, 	MPI_COMM_WORLD);
		/* get iBlue */
		MPI_Ssend(&(image->iBlue), 	1, 	MPI_INT, 	iProcess, 	1, 	MPI_COMM_WORLD);
		/* get iRedBg */
		MPI_Ssend(&(image->iRedBg), 	1, 	MPI_INT, 	iProcess, 	1, 	MPI_COMM_WORLD);
		/* get iGreen */
		MPI_Ssend(&(image->iGreenBg), 	1, 	MPI_INT, 	iProcess, 	1, 	MPI_COMM_WORLD);
		/* get iBlue */
		MPI_Ssend(&(image->iBlueBg), 	1, 	MPI_INT, 	iProcess, 	1, 	MPI_COMM_WORLD);
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

	image->iWidth = 1200;
	image->iHeight = 900;
	image->uiMaxIterations = 120;
	image->iSlices = 11;
	image->dReMin = -2.0;
	image->dReMax = 1.0;
	image->dImMin = -1.2;
	image->FileName = "mandelbrot.bmp";
	image->iRed = 255;
	image->iGreen = 255;
	image->iBlue = 0;
	image->iRedBg = 0;
	image->iGreenBg = 0;
	image->iBlueBg = 0;

}

void finalizeConfig(TImageConfig* image) {

	/* auto compute some values */
	image->dImMax = image->dImMin + (image->dReMax - image->dReMin) * image->iHeight / image->iWidth;
	image->dReFactor = (image->dReMax - image->dReMin) / (image->iWidth - 1);
	image->dImFactor = (image->dImMax - image->dImMin) / (image->iHeight - 1);

}

void bestpictureConfig(TImageConfig* config) {

	/* only change, what will be different from the basic config */
	config->iWidth = 1200;
	config->iHeight= 900;
	config->uiMaxIterations = 255;
	config->dReMin = -0.3871904296875;
	config->dReMax = -0.3858232421874879;
	config->dImMin = 0.6238156738281248;
	config->iRedBg = 255;
	config->iGreenBg = 255;
	config->iBlueBg = 255;

}

void dialog(TImageConfig* config) {

	printf("here would the dialog start\n");	
}

void printUsageInformation() {

	printf("mandelbrot generator in mpi-conform C\n"
		"(c) 2012 R. Fruth, S. Kaiser & E. Kuhnt\n\n"
		"usage information:\n"
		"\tmpirun program [options]\n\n"
		"options:\n"
		"long\t\tshort\targ\t\tdefault\t\teffect\n"
		"--dialog\t-d\tnone\t\tnot set\t\tenable the configuration dialog\n"
		"--bestpic\t-p\tnone\t\tnot set\t\tgenerate the best picture\n"
		"--help\t\t-u\tnone\t\tnot set\t\tdisplay this information\n"
		"--width\t\t-w\tint\t\t1200\t\tset the picture width\n"
		"--height\t-h\tint\t\t900\t\tset the picture height\n"
		"--iterations\t-n\tunsigned\t120\t\tset the iterations max\n"
		"--immin\t\t-a\tdouble\t\t-1.2\t\tset the minimal imaginary value\n"
		"--remin\t\t-i\tdouble\t\t-2.0\t\tset the minimal real value\n"
		"--remax\t\t-j\tdouble\t\t1.0\t\tset the maximal real value\n"
		"--red\t\t-r\tint\t\t255\t\t set the red value\n"
		"--green\t\t-g\tint\t\t255\t\t set the green value\n"
		"--blue\t\t-b\tint\t\t0\t\t set the blue value\n"
		"--redBg\t\t-y\tint\t\t0\t\t set the background red value\n"
		"--greenBg\t\t-x\tint\t\t0\t\t set the background green value\n"
		"--blueBg\t\t-z\tint\t\t0\t\t set the background blue value\n"
		"--filename\t-f\tchar[]\t\tmandelbrot.bmp\tset the filename\n"
		"--slices\t-s\tint\t\t11\t\tset the number of slices\n\n"
		"note: values not mentioned above will be computed automatically for you!"
		"\n"
		);
	exit(0);
}

void parseCommandLineParameters(int argc, char** argv, TImageConfig* config) {

	int iIndex, bDialogFlag = 0, bBestPicFlag = 0;

	/* TODO: extract underflow test */
	/* we will write directly to the config in order to prevent more memory allocations */ 
	while (1) {
		/* getopt_long stores the option index here. */
		int iOptionIndex = 0;
		/* options have to be declared here -- otherwise they will not work as expected */
		struct option oLongOptions[] = {
			/* These options set a flag. */
			{"dialog", 	no_argument,       	0, 	'd'},
			{"bestpic",	no_argument,	   	0, 	'p'},
			{"help",	no_argument,		0,	'u'},
			/* These options don't set a flag.
			   We distinguish them by their indices. */
			{"width",     	required_argument,      0, 	'w'},
			{"height",  	required_argument,      0, 	'h'},
			{"iterations",  required_argument, 	0, 	'n'},
			{"immin",  	required_argument, 	0, 	'a'},
			/* skipping one indice in order be able to improve later */
			{"remin",	required_argument,	0, 	'i'},
			{"remax",	required_argument,	0, 	'j'},
			{"slices",	required_argument,	0, 	's'},
			{"filename",    required_argument, 	0,	'f'},
			{"red",		required_argument,	0,	'r'},
			{"green",	required_argument,	0,	'g'},
			{"blue",	required_argument,	0,	'b'},
			{"redBg",	required_argument,	0,	'x'},
			{"greenBg",	required_argument,	0,	'y'},
			{"blueBg",	required_argument,	0,	'z'},			
		       	{0, 		0, 			0,	0}
		};
		/* get the next option */
		iIndex = getopt_long (argc, argv, "udpw:h:n:a:i:j:s:f:r:g:b:x:z:y:", oLongOptions, &iOptionIndex);
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
				if(config->iWidth < 0) {
					fprintf(stderr, "ERROR: You may not set a negative value for --width / -w !\n");
					exit(1);
				}
               			break;
             		case 'h':
				config->iHeight = atoi(optarg);
				if(config->iHeight < 0) {
					fprintf(stderr, "ERROR: You may not set a negative value for --height / -h !\n");
					exit(1);
				}				
               			break;
             		case 'n':
				/* prevent unsinged underflow */
				if(optarg[0]=='-') {
					fprintf(stderr, "ERROR: You may not set a negative value for --iterations / -n !\n");
					exit(1);
				}
				config->uiMaxIterations = (unsigned int) atoi(optarg);
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
				config->FileName = optarg;
				break;
			case 'r':
				config->iRed = atoi(optarg);
				if(config->iRed < 0) {
					fprintf(stderr, "ERROR: You may not set a negative value for --red / -r !\n");
					exit(1);
				}
				break;
			case 'g':
				config->iGreen = atoi(optarg);
				if(config->iGreen < 0) {
					fprintf(stderr, "ERROR: You may not set a negative value for --green / -g !\n");
					exit(1);
				}
				break;
			case 'b':
				config->iBlue = atoi(optarg);
				if(config->iBlue < 0) {
					fprintf(stderr, "ERROR: You may not set a negative value for --blue / -b !\n");
					exit(1);
				}
				break;
			case 'x':
				config->iRedBg = atoi(optarg);
				if(config->iRedBg < 0) {
					fprintf(stderr, "ERROR: You may not set a negative value for --redBg / -x !\n");
					exit(1);
				}
				break;
			case 'y':
				config->iGreenBg = atoi(optarg);
				if(config->iGreenBg < 0) {
					fprintf(stderr, "ERROR: You may not set a negative value for --greenBg / -y !\n");
					exit(1);
				}
				break;
			case 'z':
				config->iBlueBg = atoi(optarg);
				if(config->iBlueBg < 0) {
					fprintf(stderr, "ERROR: You may not set a negative value for --blueBg / -z !\n");
					exit(1);
				}
				break;
			case 's':
				config->iSlices = atoi(optarg);
				if(config->iSlices < 0) {
					fprintf(stderr, "ERROR: You may not set a negative value for --slices / -s !\n");
					exit(1);
				}
				break;
			case 'u':
				printUsageInformation();
				break;
             		case '?':
               			/* getopt_long already printed an default error message. */
				fprintf(stderr, "type 'mpirun program --help' for usage information!\n");
				exit(1);
             		default:
               			exit(1);
		}
	}

	/* configuring bestpic */
	if(bBestPicFlag == 1) {
		bestpictureConfig(config);
	}
	/* using the dialog
	   note: if we have set the bestpic before (e.g. the user sets both flags),
		 the dialog defaults will show the bestpic configuration
		 this enables us to work-around an error message if both flags are 
		 set
	 */
	if(bDialogFlag == 1){
		dialog(config);
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
