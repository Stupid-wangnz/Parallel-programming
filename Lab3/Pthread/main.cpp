#include <iostream>
#include<pthread.h>
#include<stdio.h>
#include<ctime>
#include<stdlib.h>
#include<cstdlib>
#include<sys/time.h>
#include<thread>
using namespace std;


int Arr_size=12;
float **A;
int MAXTHREADs=4;
void reset(float** G)
{

    for(int i=0;i<Arr_size;i++)
        for(int j=0;j<Arr_size;j++)
            A[i][j]=G[i][j];
}

typedef struct{
    int k;//消去的轮次
    int t_id;//线程id
}threadParam_t;

typedef struct{
    int k;//消去的轮次
    int t_id;//线程id
    int threads;
}threadParam_t1;
void *threadFunc(void*param){
    threadParam_t*p=(threadParam_t*)param;

    int k=p->k;
    int t_id=p->t_id;
    int i=k+t_id+1;//计算任务

    for(int j=k+1;j<Arr_size;j++)
       A[i][j]=A[i][j]-A[i][k]*A[k][j];

    A[i][k]=0;
    pthread_exit(NULL);

}

void *threadFunc1(void*param){
    threadParam_t1*p=(threadParam_t1*)param;

    int k=p->k;
    int t_id=p->t_id;
    int i=k+t_id+1;//计算任务
    int step=p->threads;

    for(int i=k+t_id+1;i<Arr_size;i+=step){
        for(int j=k+1;j<Arr_size;j++)
            A[i][j]=A[i][j]-A[i][k]*A[k][j];
        A[i][k]=0;
    }

    pthread_exit(NULL);

}

void Serial()
{
    for(int k=0;k<Arr_size;k++)
    {
        for(int j=k+1;j<Arr_size;j++)
        {
            A[k][j]/=A[k][k];
        }

        A[k][k]=1.0;

        for(int i=k+1;i<Arr_size;i++)
        {
            for(int j=k+1;j<Arr_size;j++)
            {
                A[i][j]-=A[i][k]*A[k][j];
            }

            A[i][k]=0.0;
        }

    }

}

void Dynamic_thread()
{

    for(int k=0;k<Arr_size;k++){
        for(int j=k+1;j<Arr_size;j++)
            A[k][j]/=A[k][k];

        A[k][k]=1;

        int worker_count=Arr_size-1-k;
        //动态创建线程
        pthread_t*handles=new pthread_t[worker_count];
        threadParam_t*param=new threadParam_t[worker_count];

        for(int i=0;i<worker_count;i++){
            param[i].k=k;
            param[i].t_id=i;

        }

        for(int i=0;i<worker_count;i++){
            pthread_create(&handles[i],NULL,threadFunc,&param[i]);
        }
        for(int i=0;i<worker_count;i++){
            pthread_join(handles[i],NULL);
        }
    }

}

void Dynamic1_thread()
{
    for(int k=0;k<Arr_size;k++){
        for(int j=k+1;j<Arr_size;j++)
            A[k][j]/=A[k][k];

        A[k][k]=1;

        int rows=Arr_size-(k+1);
        if(rows==0)
            break;

        int threadCount=min(rows,MAXTHREADs);
        pthread_t*handles=new pthread_t[threadCount];
        threadParam_t1 *param=new threadParam_t1[threadCount];
        for(int i=0;i<threadCount;i++){
            param[i].k=k;
            param[i].t_id=i;
            param[i].threads=threadCount;
        }
        for(int i=0;i<threadCount;i++)
            pthread_create(&handles[i],NULL,threadFunc1,&param[i]);
        for(int i=0;i<threadCount;i++)
            pthread_join(handles[i],NULL);

    }

}

int main()
{
    A=new float*[Arr_size];
    for(int i=0;i<Arr_size;i++)
        A[i]=new float[Arr_size];

    float** Gauss_arr;
    Gauss_arr=new float*[Arr_size];
    for(int i=0;i<Arr_size;i++)
        Gauss_arr[i]=new float[Arr_size];

    for(int i=0;i<Arr_size;i++)
    {
        for(int j=0;j<i;j++)
            Gauss_arr[i][j]=0;
        Gauss_arr[i][i]=1.0;
        for(int j=i;j<Arr_size;j++)
            Gauss_arr[i][j]=rand();
    }
    for(int k=0;k<Arr_size;k++)
        for(int i=k+1;i<Arr_size;i++)
            for(int j=0;j<Arr_size;j++)
                Gauss_arr[i][j]+=Gauss_arr[k][j];

    double time=0;
    struct timeval tv_begin,tv_end;

    reset(Gauss_arr);
    gettimeofday(&tv_begin,NULL);
    Serial();
    gettimeofday(&tv_end,NULL);

    time+=(tv_end.tv_sec-tv_begin.tv_sec)*1000.0+(tv_end.tv_usec-tv_begin.tv_usec)/1000.0;
    cout<<time<<"ms"<<endl;

    reset(Gauss_arr);

    time=0;
    gettimeofday(&tv_begin,NULL);
    Dynamic_thread();
    gettimeofday(&tv_end,NULL);

    time+=(tv_end.tv_sec-tv_begin.tv_sec)*1000.0+(tv_end.tv_usec-tv_begin.tv_usec)/1000.0;
    cout<<time<<"ms"<<endl;
    reset(Gauss_arr);

    time=0;
    gettimeofday(&tv_begin,NULL);
    Dynamic1_thread();
    gettimeofday(&tv_end,NULL);

    time+=(tv_end.tv_sec-tv_begin.tv_sec)*1000.0+(tv_end.tv_usec-tv_begin.tv_usec)/1000.0;
    cout<<time<<"ms"<<endl;
    reset(Gauss_arr);


    return 0;
}
