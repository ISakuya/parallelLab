#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <mpi.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "unistd.h"
#include <algorithm>
#include <time.h>
#include <string.h>
#define BLOCK 1024*1024
#define ERO if(!myRank) 
class OmpPsrs
{
public:
	OmpPsrs(){
		vdata = NULL;
		sampleList = NULL;
		primarySamples = NULL;
		sortPoint = NULL;
		primaryDividePoint = NULL;
	}

	~OmpPsrs(){};
	//运行OpenMP实现的PSRS排序
	double* run(double* dataSet, int dataSize) {
		double* swapData = NULL;
		thread_count = 32;
		todo = dataSize / thread_count;

		vdata = dataSet;

		sampleList = (double*)malloc(sizeof(double)*thread_count*thread_count);
	#	pragma omp parallel num_threads(thread_count)
		selectSamples();


		//select samples
		std::sort(sampleList, sampleList + thread_count*thread_count);
		primarySamples = (double*)malloc(sizeof(double)*(thread_count - 1));
		for (int i = 0;i < thread_count - 1;++i) primarySamples[i] = sampleList[(i + 1)*thread_count];
		primaryDividePoint = (double**)malloc(sizeof(double*)*(thread_count*thread_count + 1));

	#	pragma omp parallel num_threads(thread_count)
		primaryDivide();

		//global swap
		primaryDividePoint[thread_count*thread_count] = vdata + dataSize;

		int pointer = 0;
		sortPoint = (int*)malloc(sizeof(int)*(thread_count + 1));
		swapData = (double*)malloc(sizeof(double)*dataSize);
		sortPoint[0] = 0;
		for (int processer = 0; processer < thread_count; ++processer) {
			for (int segment = 0; segment < thread_count; ++segment)
			{
				double* head = primaryDividePoint[segment*thread_count + processer];
				double* tail = primaryDividePoint[segment*thread_count + processer + 1];
				for (double* i = head;i < tail;++i) {
					swapData[pointer++] = *i;
				}
			}
			sortPoint[processer + 1] = pointer;
		}
		free(vdata);

	#	pragma omp parallel num_threads(thread_count)
		lastSort(swapData);
		//printf("time to OmpPsrs sort :%fs\n", (double)(stop - start)/CLOCKS_PER_SEC);
		//printf("%d\n", dataSet[todo-1]);
		//printf("%d\n", dataSet[0]);

		/*bool check = true;
		for (int i = 0;i < dataSize - 1;++i) if (swapData[i]>swapData[i + 1]) check = false;
		if (check) printf("pass check\n\n");*/

		return swapData;
	}
private:
	int 			todo; 				/*各线程分配的数据数量	*/
	int 			thread_count;  		/*线程数				*/
	double* 		vdata;				/*各线程数据 			*/
	double* 		sampleList; 		/*采样数组				*/
	double* 		primarySamples;		/*主元数组				*/
	int* 			sortPoint;			/*全局交换点			*/
	double** 		primaryDividePoint;	/*划分点				*/
	
	void selectSamples() {
		//选择样本
		int my_rank = omp_get_thread_num();
		int start = my_rank * todo;
		int end = start + todo;

		std::sort(vdata + start, vdata + end);

		int sampleStride = todo / thread_count;
		for (int i = 0;i < thread_count;++i) sampleList[my_rank*thread_count + i] = vdata[my_rank*todo + sampleStride*i];
	}
	void primaryDivide() {
		//主元划分
		int my_rank = omp_get_thread_num();
		int start = my_rank * todo;

		primaryDividePoint[my_rank*thread_count] = vdata + start;
		for (int i = 0;i < thread_count - 1;++i) {
			double* pos = std::lower_bound(vdata + start, vdata + start + todo, primarySamples[i]);
			while (*pos == primarySamples[i]) {
				if (pos >= vdata + start + todo) break;
				pos++;
				pos = std::lower_bound(pos, vdata + start + todo, primarySamples[i]);
			}
			primaryDividePoint[my_rank*thread_count + 1 + i] = pos;
		}
	}

	void lastSort(double* swapData) {
		//最后阶段排序
		int my_rank = omp_get_thread_num();
		std::sort(swapData + sortPoint[my_rank], swapData + sortPoint[my_rank + 1]);
	}
};

int         commSize;             		/* Number of processes   */
int         myRank;              		/* My process rank       */
double*		data = NULL;		  		/* 汇总数据				*/
double*	    subData = NULL;		 		/* 各进程内数据			*/
double*		newData = NULL;		  		/* 全局交换后数据		*/
int		    todo;						/* 各进程最初分配的数据	*/
long		totalSize;					/* 数据总个数			*/
double*		samples = NULL;				/*采样 					*/
double*		primarySamples = NULL;		/*主元 					*/
double**	primaryDividePoint = NULL;	/*主元划分点 			*/
int*		subSegmentLength = NULL;	/*全局交换长度 			*/
int*		segmentLength = NULL;
int			newSize;					/*全局交换后数据个数 	*/
OmpPsrs 	ompPsrs = OmpPsrs();		/*OpenMp排序 			*/

void selectSamples() {
	subData = ompPsrs.run(subData, todo);
	//std::sort(subData, todo);
	int sampleStride = todo / commSize;
	double* selectedSample = (double*)malloc(sizeof(double)*commSize);
	for (int i = 0;i < commSize;++i) {
		selectedSample[i] = subData[i*sampleStride];
		//printf("%d\n", selectedSample[i]);
	}
	MPI_Gather(selectedSample, commSize, MPI_DOUBLE, samples, commSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void selectPrimarySamples() {
	//select primary samples
	primarySamples = (double*)malloc(sizeof(double)*(commSize - 1));
	ERO {
		samples = ompPsrs.run(samples, commSize*commSize);
		//std::sort(samples, samples + commSize*commSize);
		//for (int i = 0;i < commSize*commSize;i++) printf("%d  ", samples[i]);
		for (int i = 0;i < commSize - 1;++i) {
			primarySamples[i] = samples[(i + 1)*commSize];
			//printf("%d   ", primarySamples[i]);
		}
	}
	MPI_Bcast(primarySamples, commSize - 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void primaryDivide() {
	subSegmentLength = (int*)malloc(sizeof(int) * commSize);
	segmentLength = (int*)malloc(sizeof(int)*commSize*commSize);
	primaryDividePoint = (double**)malloc(sizeof(double*) * (commSize + 1));
	primaryDividePoint[0] = subData;

	for (int i = 0;i < commSize - 1;++i) {
		double* pos = std::lower_bound(subData, subData + todo, primarySamples[i]);
		while (*pos == primarySamples[i]) {
			if (pos >= subData + todo) break;
			pos++;
			pos = std::lower_bound(pos, subData + todo, primarySamples[i]);
		}
		primaryDividePoint[i + 1] = pos;
		//printf("%d %d\n", myRank, *pos);
	}

	primaryDividePoint[commSize] = subData + todo;
	for (int i = 0;i < commSize;++i) subSegmentLength[i] = primaryDividePoint[i + 1] - primaryDividePoint[i];

	MPI_Gather(subSegmentLength, commSize, MPI_INT, segmentLength, commSize, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(segmentLength, commSize*commSize, MPI_INT, 0, MPI_COMM_WORLD);

	//if (myRank==1) for (int i = 0;i < commSize*commSize;i++) printf("%d  ", segmentLength[i]);
}

void globalSwap() {
	int* sendNum = (int*)malloc(sizeof(int)*commSize);
	int* recvNum = (int*)malloc(sizeof(int)*commSize);
	int* sendStride = (int*)malloc(sizeof(int)*commSize);
	int* recvStride = (int*)malloc(sizeof(int)*commSize);
	newSize = 0;
	for (int i = 0;i < commSize;++i) {
		sendNum[i] = segmentLength[myRank*commSize + i];
		recvNum[i] = segmentLength[i*commSize + myRank];
		newSize += recvNum[i];
	}
	sendStride[0] = 0;
	recvStride[0] = 0;
	for (int i = 1;i < commSize;++i) {
		sendStride[i] = sendStride[i - 1] + sendNum[i - 1];
		recvStride[i] = recvStride[i - 1] + recvNum[i - 1];
	}

	newData = (double*)malloc(sizeof(double)*newSize);

	MPI_Alltoallv(subData, sendNum, sendStride, MPI_DOUBLE, newData, recvNum, recvStride, MPI_DOUBLE, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	free(subData);

	/*if (myRank == 0) {
		for (int i = 0;i < newSize;i++) printf("%d\n", newData[i]);
	}*/
}

bool subSortCheck(){
	for (int i = 0;i < newSize - 1;++i) if (newData[i] > newData[i + 1]) return false;
	return true;
}

void finalSort() {
	newData = ompPsrs.run(newData, newSize);
	//std::sort(newData, newData + newSize);
	/*
	ERO{
		data = (double*)malloc(4* totalSize);
		double* ff = (double*)malloc(4* totalSize);
		//double* ff = (double*)malloc(4* totalSize);
		if (data == NULL || ff == NULL) printf("final data malloc failed\n");
	}
	int* sendNum = (int*)malloc(sizeof(int)*(commSize));
	int* sendStride = (int*)malloc(sizeof(int)*(commSize));

	MPI_Gather(&newSize, 1, MPI_INT, sendNum, 1, MPI_INT, 0, MPI_COMM_WORLD);
	sendStride[0] = 0;
	for (int i = 1;i < commSize;++i) sendStride[i] = sendStride[i - 1] + sendNum[i - 1];

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Gatherv(newData, newSize, MPI_DOUBLE, data, sendNum, sendStride, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);*/
	
}

bool sortCheck() {
	for (int i = 0;i < totalSize - 1;++i) if (data[i] > data[i + 1]) return false;
	return true;
}

int readSubData() {
	char pre[10];
	char pos[5] = ".txt";
	char filename[50] = "../../test4/";
	sprintf(pre, "%d", myRank + 1);
	strcat(pre, pos);
	strcat(filename, pre);

	subData = (double*)malloc(todo*sizeof(double));
	if (subData == NULL) printf("malloc failed\n");

	int handle = open(filename, O_RDWR);
	if (handle == -1)
	{
		printf("Cannot Open file %s\n", filename);
		return -1;
	}
	int readRes;
	readRes = read(handle, subData, todo*sizeof(double));
	if (!readRes) return 0;
	//printf("%f\n", subData[1]);
	//printf("%d\n", readRes);
	return 1;
}

int main(int argc, char* argv[]) {
	time_t start, stop, sm, sps, pd, gs, fs;
	time_t ts,te;
	/* Start up MPI */
	MPI_Init(NULL, NULL);

	/* Get the number of processes */
	MPI_Comm_size(MPI_COMM_WORLD, &commSize);

	/* Get my rank among all the processes */
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

	totalSize = 1024ll*1024*1024*10;
	todo = 167772160;//1342177280 / 8;

	ERO{
		samples = (double*)malloc(sizeof(double)*commSize*commSize);
	}

	//makeData();
	readSubData();				//读取数据

	ERO printf("load complete\n");
	ERO ts = time(NULL);
	time_t subdataTime = time(NULL);

	selectSamples();			//选择样本
	ERO{
		sm = time(NULL);
		printf("selectSamples time: %lds\n", (sm - subdataTime));
		sm = time(NULL);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	selectPrimarySamples();		//选择主元
	ERO{
		sps = time(NULL);
		printf("selectPrimarySamples time: %lds\n", (sps - sm));
		sps = time(NULL);
	}

	primaryDivide();			//主元划分
	ERO{
		pd = time(NULL);
		printf("primaryDivide time: %lds\n", (pd - sps));
		pd = time(NULL);
	}

	globalSwap();				//全局交换
	ERO{
		gs = time(NULL);
		printf("globalSwap time: %lds\n", (gs - pd));
		gs = time(NULL);
	}

	finalSort();				//子数组排序
	ERO{
		fs = time(NULL);
		printf("finalSort time: %lds\n", (fs - gs));
	}

	MPI_Barrier(MPI_COMM_WORLD);
	ERO te = time(NULL);
	//计算时间
	ERO printf("time : %lds\n", te - ts);
	if (!subSortCheck()) printf("subSort failed\n");
	free(newData);
	ERO{
		//free(data);
	}

	MPI_Finalize();

	return 0;
}