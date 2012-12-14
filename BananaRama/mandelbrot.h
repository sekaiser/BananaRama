/*
 *
 * mandelbrot.h                                                       
 *                                                                   
 * Disclaimer:                                                     
 * This implementation of Mandelbrot is inspired by an excellent      
 * tutorial given by Michel Valleries, Professor for Physics at Drexel
 * University. You can find his tutorial here:                        
 * http://www.physics.drexel.edu/~valliere/PHYS405/Content.html      
 *                                                                   
 *                                                                   
 * Inspired by Prof. Velleries we programmed a parallel version     
 * Mandelbrot. Therefore, we used the Master-Slave model. In order 
 * to compute the area, the master process divides  the computation 
 * area into slices. These slices will be distributed among the     
 * available slaves, who perform the compuation. Afterwards the slaves
 * send the results back to the master process, who is will compute  
 * the received data in order to generate a picture.                  
 * The user interaction is realized by implementing a small dialog.  
 * Furthermore we tried to implement a solution, that is as fast as  
 * possible.                                                         
 *                                                                   
 *                                                                                                                                     
 *  author Sebastian Kaiser (743121)                                
 *  author Eric Kuhnt                                               
 *  author Robert                                                   
 *                                                                  
 *  version 1.0  12/06/12                                           
 *                                                                   
 */

#ifndef MANDELBROT__H
#define MANDELBROT__H

/*
 * Constants
 */

/* if |z| goes beyond the threshold, point C is not in the set */
#define THRESHOLD_RADIUS 4.0

/* Maximum number of iterations before declaring a point in the
 * Mandelbrot set */
#define MAX_ITERATIONS 120

/* image size. default is 700x500 pixels */
#define NX 640
#define NY 400

/* number of slices that are going to be distributed among the slave
 * processes. */
#define N_SLICES 64

/* Maximum nuber of processes */
#define MAX_PROCESSES 4

/* #define debug */

/*
 * Functions
 */

/**
   Setup the grid in a complex C plane. Therefore, define the grid
   spanning of the complex plane for the C variable. In addition,
   define the C spacing in the complex plane (dcr, dci).
   
   param crMin The lower bound of real part of complex number C.
   param crMax The upper bound of real part of complex number C.
   param ciMin The lower bound of imaginary  part of complex number C.
   param ciMax The upper bound of imaginary  part of complex number C.
   param dcr
   param dci
 */
void setGrid(double *crMin, double *crMax, double *ciMin,
             double *ciMax, double *dcr, double *dci);


/**
   Setup the slices through which compute the Mandelbrot set.
   
   param iarrBegin  Begin of slice.
   param iarrEnd    End of slice.
   param iarrOffset 
   return Average number of slices.
 */
int setSlices(int *iarrBegin, int *iarrEnd, int *iarrOffset);


/**
   Master node -- control the calculation.

   param numProcs  Number of processes.
   param rank      Rank of current process.
   param iarrBegin 
   param iarrEnd
   param average   Average number of slices.
   param crMin     The lower bound of real part of complex number C.
   param ciMin     The lower bound of imaginary  part of complex
                    number C.
   param dcr
   param dci
   param storage   Allocated memory for a slice.
*/
void master(int mpiSize, int *iarrBegin, int *iarrEnd,
            int average, double *storage, int *offset );

/**
   Iterate the Mandelbrot map. Within each computation step, we
   need to check whether the given input point in complex C plane
   belongs to the set or not. 

   param cReal The real part of the complex number.
   param cImg  The imaginary part of the complex number.
   param count Iteration number for all scanned values.
   return The corresponding color value.
 */
double iterate(double cReal, double cImg);


/**
   Calculate the Mandelbrot set in slice.

   param slice     Number of current slide.
   param iarrBegin
   param iarrEnd
   param count
   param crMin
   param ciMin
   param dcr
   param dci
   param storage   Data of current slide.
 */
void computeSlice(int slice, int *iarrBegin, int *iarrEnd,
                  double crMin, double ciMin,
                  double dcr, double dci, double *storage);


/**
   Slave processes. They calculate a specific slice as
   received from process 0.
   
   param rank       Rank of slave process
   param iarrBegin 
   param iarrEnd
   param average
   param crMin
   param ciMin
   param dcr
   param dci
   param storage
 */
#ifdef debug
void slave(int rank, int *iarrBegin, int *iarrEnd, int average, double crMin, double ciMin,
           	   double dcr, double dci, double *storage);
#else
void slave(int *iarrBegin, int *iarrEnd, int average, double crMin, double ciMin, double dcr, double dci, double* storage);
#endif

/**
   Kill the slave processes.
   param mpiSize
 */
void closeMPI(int mpiSize);
#endif  /* MANDELBROT__H */
