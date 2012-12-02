#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv)
{
    int rank, size, i;
    char buf[256], buffer[256];
    MPI_Status status;				/* MPI Status         */
    MPI_Init(&argc, &argv);			/* MPI Initialisieren */
    MPI_Comm_size(MPI_COMM_WORLD, &size);	/* MPI size abfragen  */
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);	/* MPI rank abfragen  */
    if(rank==0)
    {
      sprintf(buffer, "Hello, World! I am a MPI Master Process running on ");
      gethostname(buffer + strlen(buffer), 100);
      printf("From 0 of %d: %s\n", size, buffer);

      for(i=1; i<size; i++)
      {
	/* MPI Nachricht empfangen         */
	MPI_Recv(buf, 256, MPI_CHAR, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
	printf("From %d of %d: %s\n", status.MPI_SOURCE, size, buf);
      }
    }
    else
    {
      sprintf(buf, "Hello, World! I am a MPI Slave Process running on ");
      gethostname(buf + strlen(buf), 100);
      /* MPI Nachricht an Prozess 0 senden */
      MPI_Send(buf, 256, MPI_CHAR, 0, 1, MPI_COMM_WORLD);
    }
    MPI_Finalize();				/* MPI beenden        */

    return 0;
}
