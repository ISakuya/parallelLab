#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include<time.h>
#include <string.h>

#define H 0.001
#define BLOCK 1024*1024

int64_t f(int x){
     return int64_t(rand());
}

int64_t data[BLOCK];
int main(int argc,char* argv[]){
     srand((int)time(0));
     int start,loop;
     int i,j;
     if(argc!=4){
         printf("Usage:%s start_value loop_times data_file\n",argv[0]);
     }
     start=atoi(argv[1]);
     loop=atoi(argv[2]);
     int handle;
     printf("starte_value=%d loop=%d\n",start,loop);
	 
	 char FileName[10];
	 char prefix[5]=".txt";
	 for (int a=1;a<=64;a++)
	 {
		//handle=open(argv[3],O_WRONLY);
		sprintf(FileName,"%d",a);
		strcat(FileName,prefix);
		handle=open(FileName,O_WRONLY);
		if(handle==-1)
		{
			printf("Cannot Open file %s\n",FileName);
			return -1;
		}
		int s=(int)start;
		for(i=0;i<loop;i++)
		{
            for(j=0;j<BLOCK;j++)
		 {
             int64_t s1=f(s);
		     data[j]=s1;
         }
		write(handle,data,sizeof(int64_t)*BLOCK);
		}

     }
	close(handle);
	printf ("ok");
                return 0;
 }
