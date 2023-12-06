/* WSCAD - 9th Marathon of Parallel Programming 
 * Simple Brute Force Algorithm for the 
 * Traveling-Salesman Problem
 * Author: Emilio Francesquini - francesquini@ic.unicamp.br
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <mpi.h>
#include <math.h>
int min_distance;
int nb_towns;

typedef struct {
    int to_town;
    int dist;
} d_info;

d_info **d_matrix;
int *dist_to_origin;


int present (int town, int depth, int *path) {
    int i;
    for (i = 0; i < depth; i++)
        if (path[i] == town) return 1;
    return 0;
}

void tsp (int depth, int current_length, int *path, int my_rank, int n_process, MPI_Status status) {
    
    if (current_length >= min_distance) 
    {
        return;
    }

    if (depth == nb_towns) 
    {
        current_length += dist_to_origin[path[nb_towns - 1]];
        if(depth % 4 == 0)
        {
            MPI_Recv(&min_distance,1,MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
            MPI_Send(&min_distance,1,MPI_INT,0,0,MPI_COMM_WORLD);

        }
        if (current_length < min_distance)
            min_distance = current_length;
    } 
    if(depth == 1)
    {
        int town, me, dist;
        me = path[depth - 1];
        for (int i =(my_rank-1); i < nb_towns; i+=(n_process-1)) 
        {
            town = d_matrix[me][i].to_town;
            if (!present (town, depth, path)) {
                path[depth] = town;
                dist = d_matrix[me][i].dist;
                tsp (depth + 1, current_length + dist, path,my_rank,n_process,status);
            }
        }
    }
    else {
        int town, me, dist;
        me = path[depth - 1];
        for (int i = 0; i < nb_towns; i++) {
            town = d_matrix[me][i].to_town;
            if (!present (town, depth, path)) {
                path[depth] = town;
                dist = d_matrix[me][i].dist;
                tsp (depth + 1, current_length + dist, path,my_rank,n_process,status);
            }
        }
    }

}

void greedy_shortest_first_heuristic(int *x, int *y) {
    int i, j, k, dist;
    int *tempdist;

    tempdist = (int*) malloc(sizeof(int) * nb_towns);
    //Could be faster, albeit not as didactic.
    //Anyway, for tractable sizes of the problem it
    //runs almost instantaneously.
    for (i = 0; i < nb_towns; i++) {
        for (j = 0; j < nb_towns; j++) {
            int dx = x[i] - x[j];
            int dy = y[i] - y[j];
            tempdist [j] = dx * dx + dy * dy;
        }
        for (j = 0; j < nb_towns; j++) {
            int tmp = INT_MAX;
            int town = 0;
            for (k = 0; k < nb_towns; k++) {
                if (tempdist [k] < tmp) {
                    tmp = tempdist [k];
                    town = k;
                }
            }
            tempdist [town] = INT_MAX;
            d_matrix[i][j].to_town = town;
            dist = (int) sqrt (tmp);
            d_matrix[i][j].dist = dist;
            if (i == 0)
                dist_to_origin[town] = dist;
        }
    }

    free(tempdist);
}

void init_tsp() {
    int i, st;
    int *x, *y;

    min_distance = INT_MAX;

    st = scanf("%u", &nb_towns);
    if (st != 1) exit(1);
    MPI_Bcast(&nb_towns,1,MPI_INT,0,MPI_COMM_WORLD);
    
    d_matrix = (d_info**) malloc (sizeof(d_info*) * nb_towns);
    for (i = 0; i < nb_towns; i++)
        d_matrix[i] = (d_info*) malloc (sizeof(d_info) * nb_towns);
    dist_to_origin = (int*) malloc(sizeof(int) * nb_towns);
   
    x = (int*) malloc(sizeof(int) * nb_towns);
    y = (int*) malloc(sizeof(int) * nb_towns);
    

    for (i = 0; i < nb_towns; i++) {
        st = scanf("%u %u", x + i, y + i);
        if (st != 2) exit(1);
    }
    MPI_Bcast(x,nb_towns,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(y,nb_towns,MPI_INT,0,MPI_COMM_WORLD);


    // greedy_shortest_first_heuristic(x, y);
    
    // free(x);
    // free(y);
}

int run_tsp(int my_rank, int n_process, MPI_Status status) {
    int i, *path;
    if(my_rank == 0)
    {
        path = (int*) malloc(sizeof(int) * nb_towns);
        path[0] = 0;

        init_tsp();
        int end = 0;

        while(end < n_process-1)
        {
            int *buffer = malloc(sizeof(int));
            MPI_Recv(buffer,1,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
            if(*buffer < min_distance)
            {
                min_distance = *buffer;
            }
            if(status.MPI_TAG == 0)
            {
                MPI_Send(&min_distance,1,MPI_INT,status.MPI_SOURCE,0,MPI_COMM_WORLD);
            }
            if(status.MPI_TAG == 1)
            {
                end++;
            }
        }
    }
    else
    {
        min_distance = INT_MAX;
        MPI_Bcast(&nb_towns,1,MPI_INT,0,MPI_COMM_WORLD);

        path = (int*) malloc(sizeof(int) * nb_towns);
        path[0] = 0;

        d_matrix = (d_info**) malloc (sizeof(d_info*) * nb_towns);
        for (i = 0; i < nb_towns; i++)
            d_matrix[i] = (d_info*) malloc (sizeof(d_info) * nb_towns);

        dist_to_origin = (int*) malloc(sizeof(int) * nb_towns);
        
        int *x, *y;
        x = (int*) malloc(sizeof(int) * nb_towns);
        y = (int*) malloc(sizeof(int) * nb_towns);

        MPI_Bcast(x,nb_towns,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(y,nb_towns,MPI_INT,0,MPI_COMM_WORLD);
        greedy_shortest_first_heuristic(x,y);


        
        tsp (1, 0, path,my_rank,n_process,status);
        MPI_Send(&min_distance,1,MPI_INT,0,1,MPI_COMM_WORLD);

    }

    MPI_Barrier(MPI_COMM_WORLD);
    free(path);
    for (i = 0; i < nb_towns; i++)
        free(d_matrix[i]);
    free(d_matrix);

    return min_distance;
}

int main (int argc, char **argv) {
    int num_instances, st;

    int my_rank, n_process;

    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &n_process);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    

     if(my_rank == 0){
        st = scanf("%u", &num_instances);
        if (st != 1) exit(1);
    }


    printf("%d\n", run_tsp(my_rank,n_process,status));
    MPI_Finalize();
    return 0;
}