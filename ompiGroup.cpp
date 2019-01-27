#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <mpi.h>     /* For MPI functions, etc */
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

	int64_t* run(int64_t* dataSet, int dataSize) {
		int64_t* swapData = NULL;
		thread_count = 32;
		todo = dataSize / thread_count;

		vdata = dataSet;

		sampleList = (int64_t*)malloc(sizeof(int64_t)*thread_count*thread_count);
	#	pragma omp parallel num_threads(thread_count)
		selectSamples();

		//select samples
		std::sort(sampleList, sampleList + thread_count*thread_count);
		primarySamples = (int64_t*)malloc(sizeof(int64_t)*(thread_count - 1));
		for (int i = 0;i < thread_count - 1;++i) primarySamples[i] = sampleList[(i + 1)*thread_count];
		primaryDividePoint = (int64_t**)malloc(sizeof(int64_t*)*(thread_count*thread_count + 1));

	#	pragma omp parallel num_threads(thread_count)
		primaryDivide();

		//global swap
		primaryDividePoint[thread_count*thread_count] = vdata + dataSize;

		int pointer = 0;
		sortPoint = (int64_t*)malloc(sizeof(int64_t)*(thread_count + 1));
		swapData = (int64_t*)malloc(sizeof(int64_t)*dataSize);
		sortPoint[0] = 0;
		for (int processer = 0; processer < thread_count; ++processer) {
			for (int segment = 0; segment < thread_count; ++segment)
			{
				int64_t* head = primaryDividePoint[segment*thread_count + processer];
				int64_t* tail = primaryDividePoint[segment*thread_count + processer + 1];
				for (int64_t* i = head;i < tail;++i) {
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
	int 				loop; //read loopMB data
	int 				todo;
	int 				thread_count;
	int64_t* 			vdata;
	int64_t* 			sampleList; //samples array
	int64_t* 			primarySamples;
	int64_t* 			sortPoint;
	int64_t** 			primaryDividePoint;
	
	void selectSamples() {
		int my_rank = omp_get_thread_num();
		int start = my_rank * todo;
		int end = start + todo;

		std::sort(vdata + start, vdata + end);

		int sampleStride = todo / thread_count;
		for (int i = 0;i < thread_count;++i) sampleList[my_rank*thread_count + i] = vdata[my_rank*todo + sampleStride*i];
	}
	void primaryDivide() {
		int my_rank = omp_get_thread_num();
		int start = my_rank * todo;

		primaryDividePoint[my_rank*thread_count] = vdata + start;
		for (int i = 0;i < thread_count - 1;++i) {
			int64_t* pos = std::lower_bound(vdata + start, vdata + start + todo, primarySamples[i]);
			while (*pos == primarySamples[i]) {
				if (pos >= vdata + start + todo) break;
				pos++;
				pos = std::lower_bound(pos, vdata + start + todo, primarySamples[i]);
			}
			primaryDividePoint[my_rank*thread_count + 1 + i] = pos;
		}
	}

	void lastSort(int64_t* swapData) {
		int my_rank = omp_get_thread_num();
		std::sort(swapData + sortPoint[my_rank], swapData + sortPoint[my_rank + 1]);
	}
};

int         	commSize;             			/* Number of processes    	*/
int         	myRank;              		 	/* My process rank        	*/
int		    	todo;							/* 各进程分配的数据个数		*/
int64_t*		data = NULL;
int64_t*	    subData = NULL;		  			/* 各进程内数据				*/
int64_t*		newData = NULL;		  			/* 全局交换后数据				*/
int64_t			totalSize;						/* 总字节数					*/	
int64_t*		samples = NULL;					/* 采样列表					*/
int64_t*		primarySamples = NULL;			/* 主元列表					*/
int64_t**		primaryDividePoint = NULL;		/* 主元划分点				*/
int64_t*		subGroupNumList = NULL;			/* 子数组的分组号列表			*/
int64_t*		groupNumPresum = NULL;			/* 计算分组号用到的前缀和		*/
int64_t*		groupNum = NULL;				/* 各个子数组分组数量			*/
int*			subSegmentLength = NULL;		/* 子数组交换段长度			*/
int*			segmentLength = NULL;			/* 汇总的子数组交换段长度		*/
int				newSize;						/* 全局交换后数据个数			*/
OmpPsrs 		ompPsrs = OmpPsrs();			/* openMP实现的PSRS排序    	*/

void selectSamples() {
	subData = ompPsrs.run(subData, todo);
	//std::sort(subData, todo);
	int sampleStride = todo / commSize;
	int64_t* selectedSample = (int64_t*)malloc(sizeof(int64_t)*commSize);
	for (int i = 0;i < commSize;++i) {
		selectedSample[i] = subData[i*sampleStride];
		//printf("%d\n", selectedSample[i]);
	}
	MPI_Gather(selectedSample, commSize, MPI_LONG, samples, commSize, MPI_LONG, 0, MPI_COMM_WORLD);
}

void selectPrimarySamples() {
	//select primary samples
	primarySamples = (int64_t*)malloc(sizeof(int64_t)*(commSize - 1));
	ERO {
		samples = ompPsrs.run(samples, commSize*commSize);
		//std::sort(samples, samples + commSize*commSize);
		//for (int i = 0;i < commSize*commSize;i++) printf("%d  ", samples[i]);
		for (int i = 0;i < commSize - 1;++i) {
			primarySamples[i] = samples[(i + 1)*commSize];
			//printf("%d   ", primarySamples[i]);
		}
	}
	MPI_Bcast(primarySamples, commSize - 1, MPI_LONG, 0, MPI_COMM_WORLD);
}

void primaryDivide() {
	subSegmentLength = (int*)malloc(sizeof(int) * commSize);
	segmentLength = (int*)malloc(sizeof(int)*commSize*commSize);
	primaryDividePoint = (int64_t**)malloc(sizeof(int64_t*) * (commSize + 1));
	primaryDividePoint[0] = subData;

	for (int i = 0;i < commSize - 1;++i) {
		int64_t* pos = std::lower_bound(subData, subData + todo, primarySamples[i]);
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

	newData = (int64_t*)malloc(sizeof(int64_t)*newSize);

	MPI_Alltoallv(subData, sendNum, sendStride, MPI_LONG, newData, recvNum, recvStride, MPI_LONG, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	free(subData);

	/*if (myRank == 0) {
		for (int i = 0;i < newSize;i++) printf("%d\n", newData[i]);
	}*/
}

bool subSortCheck(){
	for (int64_t i = 0;i < newSize - 1;++i) if (newData[i] > newData[i + 1]) return false;
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
	MPI_Gatherv(newData, newSize, MPI_LONG, data, sendNum, sendStride, MPI_LONG, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);*/
	
}

void allocSubGroupNum(){
	//各个进程内分配分组号
	int64_t subGroupNum = 0;
	subGroupNumList = (int64_t*)malloc(sizeof(int64_t)*newSize);
	if (newSize < 0) return;
	subGroupNumList[0] = subGroupNum;
	for(int64_t i=1;i<newSize;i++){
		if (newData[i]!=newData[i-1]) subGroupNum++;
		subGroupNumList[i] = subGroupNum;
	}
	subGroupNum++;
	MPI_Gather(&subGroupNum, 1, MPI_LONG, groupNum, 1, MPI_LONG, 0, MPI_COMM_WORLD);
}
void allocGloGroupNum(){
	//更改局部分组号，生成全局分组号
	for(int64_t i=1;i<newSize;i++) subGroupNumList[i] += groupNumPresum[0];
}

void allocateGroupNum(){
	//生成分组号
	ERO groupNum = (int64_t*)malloc(sizeof(int64_t)*commSize);
	allocSubGroupNum();
	free(newData);
	groupNumPresum = (int64_t*)malloc(sizeof(int64_t)*commSize);
	ERO{
		groupNumPresum[0] = 0;
		for(int i=1;i<commSize;i++)
			groupNumPresum[i] = groupNumPresum[i-1] + groupNum[i-1];
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(groupNumPresum, commSize, MPI_LONG, 0, MPI_COMM_WORLD);
	allocGloGroupNum();
}

bool sortCheck() {
	for (int64_t i = 0;i < totalSize - 1;++i) if (data[i] > data[i + 1]) return false;
	return true;
}

int readSubData() {
	//读取数据
	char pre[10];
	char pos[5] = ".txt";
	char filename[50] = "../../dazuoye/";
	sprintf(pre, "%d", myRank + 1);
	strcat(pre, pos);
	strcat(filename, pre);

	subData = (int64_t*)malloc(todo*sizeof(int64_t));
	if (subData == NULL) printf("malloc failed\n");

	int handle = open(filename, O_RDWR);
	if (handle == -1)
	{
		printf("Cannot Open file %s\n", filename);
		return -1;
	}
	int readRes;
	readRes = read(handle, subData, todo*sizeof(int64_t));
	if (readRes != todo*sizeof(int64_t)) return 0;
	return 1;
}

int main(int argc, char* argv[]) {
	time_t start, stop, sm, sps, pd, gs, fs, fe, subs, sube;
	time_t ts,te;
	/* Start up MPI */
	MPI_Init(NULL, NULL);

	/* Get the number of processes */
	MPI_Comm_size(MPI_COMM_WORLD, &commSize);

	/* Get my rank among all the processes */
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

	totalSize = 1024ll*1024*1024*8;
	todo = 134217728;

	ERO{
		samples = (int64_t*)malloc(sizeof(int64_t)*commSize*commSize);
	}

	readSubData();		//读取数据
	ERO printf("load complete\n");
	ERO ts = time(NULL);
	time_t subdataTime = time(NULL);

	selectSamples();
	ERO{
		sm = time(NULL);
		printf("selectSamples time: %lds\n", (sm - subdataTime));
		sm = time(NULL);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	selectPrimarySamples();
	ERO{
		sps = time(NULL);
		printf("selectPrimarySamples time: %lds\n", (sps - sm));
		sps = time(NULL);
	}

	primaryDivide();
	ERO{
		pd = time(NULL);
		printf("primaryDivide time: %lds\n", (pd - sps));
		pd = time(NULL);
	}

	globalSwap();
	ERO{
		gs = time(NULL);
		printf("globalSwap time: %lds\n", (gs - pd));
		gs = time(NULL);
	}

	finalSort();
	ERO{
		fs = time(NULL);
		printf("finalSort time: %lds\n", (fs - gs));
	}

	MPI_Barrier(MPI_COMM_WORLD);
	ERO fe = time(NULL);

	ERO printf("PSRS time : %lds\n", fe - ts);

	ERO subs = time(NULL);
	//if (!subSortCheck()) printf("subSort failed\n");
	allocateGroupNum();		//分配分组号
	ERO sube = time(NULL);

	ERO{
		printf("allocate time : %lds\n", sube - subs);
		printf("total time : %lds\n", sube - ts);
	}
	MPI_Finalize();

	return 0;
}