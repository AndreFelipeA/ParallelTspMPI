#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <limits.h>
int main(int argc, char **argv)
{
    int my_rank, n_process;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &n_process);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    int *buffer = malloc(sizeof(int));
    MPI_Request request;
    int min_distance = INT_MAX;
    if(my_rank == 0)
    {
        for(int i = 0; i < 10; i++)
        {
            printf("vou receber:%d\n",i);
            MPI_Recv(buffer,1,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
            if(status.MPI_TAG == 1)
            {
                printf("so envia");
            }
            else
            {
                if(*buffer < min_distance)
                {
                    min_distance = *buffer;
                    printf("min_distance:%d\n",min_distance);
                }
                printf("envia e recebe");
                MPI_Send(&min_distance,1,MPI_INT,status.MPI_SOURCE,0,MPI_COMM_WORLD);
            }

        }
        

    }
    else
    {
        for(int i=10; i > 0; i--)
        {
            min_distance = i;
            printf("vou enviar:%d\n",i);
            MPI_Send(&min_distance,1,MPI_INT,0,1,MPI_COMM_WORLD);
            printf("envieiiii\n");
            MPI_Recv(&min_distance,1,MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);

        }
    }

    MPI_Finalize();
    return 0;
}