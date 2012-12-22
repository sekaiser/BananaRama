/*
	mandelbrot generator in mpi-conform C
	(c) 2012 by R. Fruth, S. Kaiser, E. Kuhnt
	
	for details please see the attached readme
*/
#ifndef __mandelbrot__h__
#define __mandelbrot__h__

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

Tuchar* getBmpFileHeader(int* filesize);

Tuchar* getBmpInfoHeader(int* iImageWidth, int* iImageHeight);

FILE* getBmpFileHandler(const char* filename);

void writeBmp(int* iImageWidth, int* iImageHeight, const char* FileName, const Tuchar* pImage);

void iterateAndStoreAPoint(TImageConfig* image, int* x, int* y, double* dZre, double* dZim, double* dCre, 
			   double* dCim, Tuchar* img);

void computeSlice(TImageConfig* image, TSlice* mySlice, Tuchar* img);

void initializeSlices(int* iHeight, int* iSlices, TSlice* sliceDimensions);

#ifdef __mandelbrot__debug__
void slaveReceiveSlices(int* rank, TImageConfig* image, TSlice* sliceDimensions, Tuchar* buffer);
#else
void slaveReceiveSlices(TImageConfig* image, TSlice* sliceDimensions, Tuchar* buffer);
#endif

void slaveReceiveConfig(TImageConfig* config);

#ifdef __mandelbrot__debug__
void slave(int* rank);
#else
void slave();
#endif

void configureProcesses(int* iSize, TImageConfig* image);

int initializeProcesses(int* iSize, int* iSlices, int* iSliceCache);

void processSlices(TImageConfig* image, int* iSlice, TSlice* sliceDimensions, 
			int* iSliceCache, Tuchar* buffer, Tuchar* out);

void master(int* iSize, TImageConfig* image);

void initializeConfig(TImageConfig* image);

void finalizeConfig(TImageConfig* image);

void bestpictureConfig(TImageConfig* config);

void dialog(TImageConfig* config);
void printUsageInformation();

void underflowCheckInteger(int* iInteger, char* cLParam, char* cSParam);

void parseCommandLineParameters(int argc, char** argv, TImageConfig* config);

int main(int argc, char **argv);

#endif

