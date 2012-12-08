#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_ITERATIONS 100
#define NX 320                 
#define NY 320


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



int assemble_work(char *work_result)
{
    // Allocate memory for resulting image
    int bitmap_pixel_count = NX * NY;
    char *bitmap = malloc(bitmap_pixel_count * 3 * sizeof(char));
    if (bitmap == NULL)
        return printf("malloc\n");
    
    // Set pixel colors according to computed limit hit
    int count;
    for (count = 0; count < NX * NY; count++)
    {
        char color = work_result[count] == 0 ? 0 : 255;
        bitmap[count * 3] = (color / 65536) % 256;
        bitmap[count * 3 + 1] = (color / 256) % 256;
        bitmap[count * 3 + 2] = (color % 256);
    }
    
    printf("hier\n");
    
    // Write image into bitmap file
    int exit_code = EXIT_SUCCESS;
    int result = write_bmp("filename.bmp", NX, NY, bitmap);
    if (result != 1)
        exit_code = error("write_bmp", result);
    
    free(bitmap);
    printf("Image creation finished\n");
    
    return(exit_code);
}

int read_bmp(const char *filename, int *width, int *height, unsigned char *rgb)
{
    printf("Sorry, reading of .bmp files isn't supported yet.\n");
    return(0);
}

int write_bmp(const char *filename, int width, int height, char *rgb)
{
    int i, j, ipos;
    int bytesPerLine;
    unsigned char *line;

    FILE *file;
    struct BMPHeader bmph;

    /* The length of each line must be a multiple of 4 bytes */

    bytesPerLine = (3 * (width + 1) / 4) * 4;

    strcpy(bmph.bfType, "BM");
    bmph.bfOffBits = 54;
    bmph.bfSize = bmph.bfOffBits + bytesPerLine * height;
    bmph.bfReserved = 0;
    bmph.biSize = 40;
    bmph.biWidth = width;
    bmph.biHeight = height;
    bmph.biPlanes = 1;
    bmph.biBitCount = 24;
    bmph.biCompression = 0;
    bmph.biSizeImage = bytesPerLine * height;
    bmph.biXPelsPerMeter = 0;
    bmph.biYPelsPerMeter = 0;
    bmph.biClrUsed = 0;       
    bmph.biClrImportant = 0; 

    file = fopen (filename, "wb");
    if (file == NULL) return(0);
  
    fwrite(&bmph.bfType, 2, 1, file);
    fwrite(&bmph.bfSize, 4, 1, file);
    fwrite(&bmph.bfReserved, 4, 1, file);
    fwrite(&bmph.bfOffBits, 4, 1, file);
    fwrite(&bmph.biSize, 4, 1, file);
    fwrite(&bmph.biWidth, 4, 1, file);
    fwrite(&bmph.biHeight, 4, 1, file);
    fwrite(&bmph.biPlanes, 2, 1, file);
    fwrite(&bmph.biBitCount, 2, 1, file);
    fwrite(&bmph.biCompression, 4, 1, file);
    fwrite(&bmph.biSizeImage, 4, 1, file);
    fwrite(&bmph.biXPelsPerMeter, 4, 1, file);
    fwrite(&bmph.biYPelsPerMeter, 4, 1, file);
    fwrite(&bmph.biClrUsed, 4, 1, file);
    fwrite(&bmph.biClrImportant, 4, 1, file);
  
    line = malloc(bytesPerLine);
    if (line == NULL)
    {
        printf("Can't allocate memory for BMP file.\n");
        return(0);
    }

    for (i = height - 1; i >= 0; i--)
    {
        for (j = 0; j < width; j++)
        {
            ipos = 3 * (width * i + j);
            line[3*j] = rgb[ipos + 2];
            line[3*j+1] = rgb[ipos + 1];
            line[3*j+2] = rgb[ipos];
        }
        fwrite(line, bytesPerLine, 1, file);
    }

    free(line);
    fclose(file);

    return(1);
}

int main(int argc, char **argv)
{
    clock_t tstart = clock();
 
    double z_real, z_img, z_magnitude;
    double c_real, c_img, crmin, crmax, cimin, cimax, 
          dcr, dci, z_current_real;

    int i, j;
    int    counter;
    int count;
    char *rgb = malloc(3*NX*NY*sizeof(char));

    crmin=-1.2548;
    crmax = -1.2544;
    cimin = 0.3816;
    cimax = 0.3822;

    dcr = (crmax - crmin)/(NX-1);
    dci = (cimax - cimin)/(NY-1);

    count = 0;
    

    for ( i=0 ; i<NX; i++)
    {

      for ( j=0 ; j<NY; j++)
      {

         c_real = crmin + i*dcr;
         c_img = cimin + j*dci;

         z_real = 0.0;
         z_img  = 0.0;

         counter = 0;
                                     

         while ( counter < MAX_ITERATIONS ) 
         {
            z_current_real = z_real;
            z_real = z_real*z_real - z_img*z_img + c_real;
            z_img  = 2.0*z_current_real*z_img + c_img;
            counter = counter + 1;


            z_magnitude = z_real*z_real+z_img*z_img;
            if ( z_magnitude > 4 ) 
  		{
				//rgb[i*j] = counter;
				rgb[i*j] = (int)((double)(255*counter)/(double)MAX_ITERATIONS);
				break;
			}
				count++;
			}

		}
	}
	
	int dummy = assemble_work(rgb);
	clock_t tend = clock();
	double cpu_time = ((float)(tend-tstart)/CLOCKS_PER_SEC);

	printf ( "finished after %f \n", cpu_time);
}
