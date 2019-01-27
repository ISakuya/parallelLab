#include "unistd.h"
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <mpi.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <time.h>
#include <string.h>
#include <algorithm>
#include <unordered_map>
#define BLOCK 1024*1024
#define ERO if(!myRank) 

int          commSize;             /* Number of processes    */
int          myRank;               /* My process rank        */
int64_t*        data;
int64_t*        subData1;
int          todo;
std::unordered_map<int64_t, int64_t> subHash1 = {};
std::unordered_map<int64_t, int64_t> subHash2 = {};

int readData() {
    char pre[10];
    char pos[5] = ".txt";
    char filename[50] = "../../dazuoye/";
    sprintf(pre, "%d", myRank + 1);
    strcat(pre, pos);
    strcat(filename, pre);

    data = (int64_t*)malloc(todo*sizeof(int64_t));
    if (data == NULL) printf("malloc failed\n");

    int handle = open(filename, O_RDWR);
    if (handle == -1)
    {
        printf("Cannot Open file %s\n", filename);
        return -1;
    }
    int readRes;
    readRes = read(handle, data, todo*sizeof(int64_t));
    if (!readRes) return 0;
    //printf("%ld\n", data[1]);
    //printf("%s\n", filename);
    return 1;
}

void initHash(){
    /*for(int i=0;i<todo/2;i++){
        if (subHash1.count(data[i])){
            subHash1[data[i]] += 1;
        }
        else{
            subHash1[data[i]] = 1;
        }
    }*/
    for(int i=todo/2;i<todo;i++){
        if (subHash2.count(data[i])){
            subHash2[data[i]] += 1;
        }
        else{
            subHash2[data[i]] = 1;
        }
    }
}
void initHash1(){
    for(int i=0;i<todo/2;i++){
        if (subHash1.count(data[i])){
            subHash1[data[i]] += 1;
        }
        else{
            subHash1[data[i]] = 1;
        }
    }
}

int main(int argc, char const *argv[])
{
    time_t start, stop;
    start = time(NULL);
    /* Start up MPI */
    MPI_Init(NULL, NULL);

    /* Get the number of processes */
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);

    /* Get my rank among all the processes */
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    todo = 1024*1024*1024/8;

    readData();
    printf("%ld\n", data[1]);
    //subData = (int64_t*)malloc(todo/2*sizeof(int64_t));
    //MPI_Scatter(data, todo/4, MPI_int64_t, subData, todo/4, MPI_int64_t, myRank-myRank%2, MPI_COMM_WORLD);
    //#   pragma omp parallel num_threads(2)
    initHash();
    subHash2.clear();
    //initHash1();
    stop = time(NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    ERO printf("time: %lds\n", stop-start);

    //initHash();
    return 0;
}