#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>

#define MAX_ITERATIONS 120
#define NX 640                 
#define NY 480


struct BMPHeader
{
    char bfType[2];       /* "BM" */
    int bfSize;           /* Size of file in bytes */
    int bfReserved;       /* set to 0 */
    int bfOffBits;        /* Byte offset to actual bitmap data (= 54) */
    int biSize;           /* Size of BITMAPINFOHEADER, in bytes (= 40) */
    int biWidth;          /* Width of image, in pixels */
    int biHeight;         /* Height of images, in pixels */
    short biPlanes;       /* Number of planes in target device (set to 1) */
    short biBitCount;     /* Bits per pixel (24 in this case) */
    int biCompression;    /* Type of compression (0 if no compression) */
    int biSizeImage;      /* Image size, in bytes (0 if no compression) */
    int biXPelsPerMeter;  /* Resolution in pixels/meter of display device */
    int biYPelsPerMeter;  /* Resolution in pixels/meter of display device */
    int biClrUsed;        /* Number of colors in the color table (if 0, use 
                             maximum allowed by biBitCount) */
    int biClrImportant;   /* Number of important colors.  If 0, all colors 
                             are important */
};

int main(int argc, char **argv)
{
    clock_t tstart = clock();
 
    double z_real, z_img, z_magnitude;
    double c_real, c_img, crmin, crmax, cimin, cimax, 
          dcr, dci, z_current_real;

    int inset;
    int i, j;
    int counter;
    double colour;
    int count;
    char rgb[NX][NY][3];
    int fd;
    char buffer[100];
    
    crmin=-1.7;
    crmax = 0.8;
    cimin = -1.0;
    cimax = 1.0;

    dcr = (crmax - crmin)/(NX-1);
    dci = (cimax - cimin)/(NY-1);

    for ( i=0 ; i<NX; i++)
    {

      for ( j=0 ; j<NY; j++)
      {

        c_real = crmin + i*dcr;
        c_img = cimin + j*dci;
        
        z_real = 0.0;
        z_img  = 0.0;
        
        counter = 0;                        
        inset = 1;
        while (counter < MAX_ITERATIONS) 
        {
          z_current_real = z_real;
          z_real = z_real*z_real - z_img*z_img + c_real;
          z_img  = 2.0*z_current_real*z_img + c_img;
          counter = counter + 1;

          z_magnitude = z_real*z_real+z_img*z_img;      
          if (z_magnitude > 4) 
          {
            inset = 0;
            colour = counter;
            count = MAX_ITERATIONS;
          }
        }
        
        if (inset)
        {
            rgb[i][j][0] = 0;
            rgb[i][j][1] = 0;
            rgb[i][j][2] = 0;
        }
        else
        { 
            rgb[i][j][0] = colour / MAX_ITERATIONS * 255;
            rgb[i][j][1] = colour / MAX_ITERATIONS * 255 / 2;
            rgb[i][j][2] = colour / MAX_ITERATIONS * 255 / 2;
        }    
      }
    }

    /* writes the data to a TGA file */
    if ((fd = open("mand.tga", O_RDWR+O_CREAT, 0)) == -1)
    {
      printf("error opening file\n");
      exit(1);
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
    write(fd, buffer, 18);
    write(fd, rgb, NX*NY*3);
    close(fd);
  
    //int dummy = assemble_work(rgb);
    clock_t tend = clock();
    double cpu_time = ((float)(tend-tstart)/CLOCKS_PER_SEC);
    
    printf ( "finished after %f \n", cpu_time);
}
