#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <fcntl.h>
#include "mandelbrot.h"
////////////////////////////////////////////////////////////////////////
//                                                                    //
// The Mandelbrot Set                                                 //
//                                                                    //
// Computation: z_{i+1} = z_{i} * z_{i} + c                           //
//                                                                    //
// The Mandelbrot Set results from a very simple map in the complex   //
// plane by following some rules:                                     //
//   1) For a given complex number c, start with z = 0, and iterate   //
//      the map above.                                                //
//                                                                    //
//   2) If z remains finite, even after an infinite number of         //
//      iterations, c belongs to M.                                   //
//                                                                    //
//   3) Repeat the procedure, or scan, for all c in the complex plane,//
//      to find the points belonging to M,                            //
//                                                                    //
//                                                                    //
// Disclaimer:                                                        //
// This implementation of Mandelbrot is inspired by an excellent      //
// tutorial given by Michel Valleries, Professor for Physics at Drexel//
// University. You can find his tutorial here:                        //
// http://www.physics.drexel.edu/~valliere/PHYS405/Content.html       //
//                                                                    //
//                                                                    //
// Inspired by Prof. Velleries we programmed a parallel version       //
// Mandelbrot. Therefore, we used the Master-Slave model. In order    //
// to compute the area, the master process divides  the computation   //
// area into slices. These slices will be distributed among the       //
// available slaves, who perform the compuation. Afterwards the slaves//
// send the results back to the master process, who is will compute   //
// the received data in order to generate a picture.                  //
// The user interaction is realized by implementing a small dialog.   //
// Furthermore we tried to implement a solution, that is as fast as   //
// possible.                                                          //
//                                                                    //
//                                                                    //
//                                                                    //
//  @author Sebastian Kaiser (743121)                                 //
//  @author Eric Kuhnt                                                //
//  @author Robert                                                    //
//                                                                    //
//  @version 1.0  12/06/12                                            //
//                                                                    //
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// Constants                                                          //
////////////////////////////////////////////////////////////////////////
#define test
////////////////////////////////////////////////////////////////////////
// Function declarations                                              //
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// Function implementations                                           //
////////////////////////////////////////////////////////////////////////
/* 
  Set up the grid in complex C plane to generate image.
  Zoom factor.
*/
void setGrid(double *crMin, double *crMax, double *ciMin,
             double *ciMax, double *dcr, double *dci) {
  /* grid spanning for variable C */
  *crMin = -1.7;
  *crMax = 0.8;
  *ciMin = -1.1;
  *ciMax = 1.1;
  /* C spacing */
  *dcr = (*crMax - *crMin)/(NX - 1);
  *dci = (*ciMax - *ciMin)/(NY - 1);
}

/* Setup the slices through which compute the Mandelbrot set  */
int setSlices(int *iarrBegin, int *iarrEnd, int *iarrOffset) {
  int slice, avg_width_slice, leftover, temp, tempOffset;
  avg_width_slice = NY / N_SLICES;
  leftover = NY % N_SLICES;
  temp = -1;
  tempOffset = 0;
  
  /* Initialize picture coordinates for slices. */
  for (slice = 0; slice < N_SLICES; slice++) {
    iarrBegin[slice] = temp + 1;
    iarrEnd[slice] = iarrBegin[slice] + avg_width_slice - 1;

    if (leftover>0) {
      iarrEnd[slice]++;
      leftover--;
    }
    
    temp = iarrEnd[slice];
    iarrOffset[slice] = tempOffset;
    tempOffset = tempOffset + NX * (iarrEnd[slice] - iarrBegin[slice] + 1);
  }
  
  return avg_width_slice;
}

/* Master node -- control the calculation */
void master(int mpiSize, int *iarrBegin, int *iarrEnd, 
            int average, double crMin, double ciMin, double dcr, 
            double dci, double *storage, int *offset ) {

  int    k, k_slice, slice, process;
  int    number, source;
  int i, j;
  static int slices[MAX_PROCESSES];
  int    kill, recv_slice;
  double mandelbrot[NX * NY];
  double rgb[NX][NY][3];
  MPI_Status status;
  char buffer[100];

  
  /* scan over slices */
  slice = 0;
  
  /* send a slice to each slave */
  for (process = 1; process < mpiSize; process++) { 
    printf("1st loop start\n");
    /* generate Mandelbrot set in each process */
    MPI_Ssend(&slice, 1, MPI_INT, process, 325, MPI_COMM_WORLD);
    slices[process] = slice;
    slice++;
    printf("1st loop end\n");
  }
  
  /* scan over received slices */
  for (k_slice = 0; k_slice < N_SLICES; k_slice++) {
    printf("2nd loop start on round %d\n", k_slice);

    /* receive a slice */
    number = NX * (average + 1);
    printf("2nd loop recv\n"); 

    MPI_Recv(storage, number, MPI_DOUBLE, MPI_ANY_SOURCE, 327, MPI_COMM_WORLD, &status);

    /* source of slice */
    source = status.MPI_SOURCE;

    /* received slice */ 
    recv_slice = slices[source];
    
    printf("2nd loop recieved slice %d\n",recv_slice);

    /* send a slice to this slave */
    if (slice < N_SLICES) {   
      printf("2nd loop sending new slice %d\n", slice);
      MPI_Ssend(&slice, 1, MPI_INT, source, 325, MPI_COMM_WORLD);
      slices[source] = slice;
      slice++;
    } else {
      kill = -1;
      
      for (process = 1; process < mpiSize; process++) {
	printf("sending kill to WORLD\n");
        MPI_Ssend(&kill, 1, MPI_INT, process, 325, MPI_COMM_WORLD);
      }
      
      fprintf(stderr, " The Mandelbrot Set will be written out \n");
      for (k = 0; k < NX * NY; k++) {
        printf(" %f \n", mandelbrot[k]);
      }
      
      FILE* file = fopen ("test.tga", "wb");

      if(file==NULL) {
        printf("ERROR: can not write to file");
      } else {

        for (i = 0; i < NX; i++) {
          for (j = 0; j < NY; j++) {
            rgb[i][j][0] = mandelbrot[i * j];
            rgb[i][j][1] = mandelbrot[i * j] / 2;
            rgb[i][j][2] = mandelbrot[i * j] / 2;
          }
        }
        
        buffer[0] = 0;
        buffer[1] = 0;
        buffer[2] = 2;
        buffer[8] = 0; buffer[9] = 0;
        buffer[10] = 0; buffer[11] = 0;
        buffer[12] = (NX & 0x00FF); buffer[13] = (NX & 0xFF00) >> 8;
        buffer[14] = (NY & 0x00FF); buffer[15] = (NY & 0xFF00) >> 8;
        buffer[16] = 24;
        buffer[17] = 0;
        fwrite(buffer,18,1,file);
        fwrite(rgb,NX*NY*3,1,file);
        fclose(file);	
      }
        
      free(storage);
      MPI_Finalize();
      exit(0);
    }
    
    /* actual number */
    number = NX * (iarrEnd[recv_slice] - iarrBegin[recv_slice] + 1);

    /* store set in matrix */
    for (k=0 ; k < number; k++ ) {
      mandelbrot[offset[recv_slice] + k] = storage[k];
    }
    
    printf("2nd loop end\n");
  }
}
                                            
double iterate(double cReal, double cImg, int *count) {
  double zReal, zImg, zCurrentReal, zMagnitude;
  double color;
  int    counter;
  
  /* z = 0 */
  zReal = 0.0;
  zImg  = 0.0;
  counter = 0;
  while (counter < MAX_ITERATIONS) {
    zCurrentReal = zReal;
    zReal = zReal*zReal - zImg * zImg + cReal;
    zImg  = 2.0 * zCurrentReal * zImg + cImg;
    counter++;
    zMagnitude = zReal * zReal + zImg * zImg;

    if (zMagnitude > THRESHOLD_RADIUS) {
      break;
    }
  }
  
  //#ifdef test
  //if (zMagnitude < THRESHOLD_RADIUS) { 
  //  printf (" %f %f  \n ", cReal, cImg );
  //}
  //#endif
  *count++;

  color = (double)(255*counter) / (double)MAX_ITERATIONS;

  return color;
}

/* Calculate the Mandelbrot set in slice  */
void computeSlice(int slice, int *iarrBegin, int *iarrEnd,
                     int *count, double crMin, double ciMin,
                     double dcr, double dci,  double *storage) {
  int i,j,k;
  double cReal, cImg;
  k = 0;
  /* scan over slice */
  for (j = iarrBegin[slice] ; j <= iarrEnd[slice]; j++) {
    for (i = 0; i < NX; i++) {
       cReal = crMin + i * dcr;
       cImg = ciMin + j * dci;

       /* Store resulting color value */
       storage[k++] = iterate(cReal, cImg, count);
    }
  }
}

void slave(int rank, int *iarrBegin, int *iarrEnd,
           int average, double crMin, double ciMin, double dcr, 
           double dci, double *storage) {
  
  static int count;
  int        number, slice;
  MPI_Status status;

  for (;;) {
    /* a new slice to calculate */
    MPI_Recv(&slice, 1, MPI_INT,0, 325, MPI_COMM_WORLD, &status);
    printf("Slave %d will process slice: %d\n.",rank, slice);
    /* suicide signal */
    if (slice < 0) {
      printf("Process %d exiting work loop.\n", rank);
      MPI_Finalize();
      exit(0);
    }
    /* calculate requested slice */
    computeSlice(slice, iarrBegin, iarrEnd, &count, crMin, ciMin, dcr, dci, storage);
    /* send results back to master */
    number = NX * (average + 1);
    MPI_Ssend(storage, number, MPI_DOUBLE, 0, 327, MPI_COMM_WORLD);
    //#ifdef test
    //fprintf( stderr, "slave %d  -->  slice %d is calculated & sent\n", rank, slice );
    //#endif
  }
}

int main(int argc, char *argv[]) {
  double crMin, crMax, ciMin, ciMax; 
  double dcr, dci;
  int    widthAvgSlice;
  int    iarrBeginSlice[1000], iarrEndSlice[1000], offset[1000];
  int    rank, mpiSize;
  double *storage;
   
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  widthAvgSlice = setSlices(iarrBeginSlice, iarrEndSlice, offset);
  setGrid(&crMin, &crMax, &ciMin, &ciMax, &dcr, &dci);
 
  /* memory allocation for color storage */
  storage = (double *)malloc(NX * (widthAvgSlice + 1) * sizeof(double));
  /* master node */
  if (rank == 0) {
    master(mpiSize, iarrBeginSlice, iarrEndSlice, widthAvgSlice,
           crMin, ciMin, dcr, dci, storage, offset);
  }
  /* slave nodes  */
  else {
    slave(rank, iarrBeginSlice, iarrEndSlice, widthAvgSlice,
          crMin, ciMin, dcr, dci, storage);
  }
  
  return 0;
}
  
