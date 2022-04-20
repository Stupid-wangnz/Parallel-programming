#include <iostream>
#include<pthread.h>
#include<semaphore.h>

using namespace std;

int Arr_size=12;
int NUM_THREADS=4;
float **A;

//信号量
sem_t sem_main;
sem_t sem_workerstart;
sem_t sem_workerend;

void reset(float** G)
{

    for(int i=0;i<Arr_size;i++)
        for(int j=0;j<Arr_size;j++)
            A[i][j]=G[i][j];
}

typedef struct{
    int t_id;//线程id
}threadParam_t;

void *threadFunc(void*param){
    threadParam_t*p=(threadParam_t*)param;

    int t_id=p->t_id;
    for(int k=0;k<Arr_size;k++){
        sem_wait(&sem_workerstart);

        for(int i=k+1+t_id;i<Arr_size;i+=NUM_THREADS){
            for(int j=k+1;j<Arr_size;j++)
                A[i][j]=A[i][j]-A[i][k]*A[k][j];

            A[i][k]=0;
        }
        sem_post(&sem_main);
        sem_wait(&sem_workerend);
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

    for(int i=0;i<Arr_size;i++){
        for(int j=0;j<Arr_size;j++)
            cout<<A[i][j]<<" ";
        cout<<endl;
    }

}

void Static_thread()
{

    sem_init(&sem_main,0,0);
    sem_init(&sem_workerstart,0,0);
    sem_init(&sem_workerend,0,0);

    pthread_t*handles=new pthread_t[NUM_THREADS];
    threadParam_t*param=new threadParam_t[NUM_THREADS];
    for(int i=0;i<NUM_THREADS;i++){
        param[i].t_id=i;
        pthread_create(&handles[i],NULL,threadFunc,&param[i]);
    }
    for(int k=0;k<Arr_size;k++){
        for(int j=k+1;j<Arr_size;j++)
            A[k][j]/=A[k][k];

        A[k][k]=1;

        for(int i=0;i<NUM_THREADS;i++)
            sem_post(&sem_workerstart);

        for(int i=0;i<NUM_THREADS;i++)
            sem_wait(&sem_main);

        for(int i=0;i<NUM_THREADS;i++)
            sem_post(&sem_workerend);
    }

    for(int i=0;i<NUM_THREADS;i++)
        pthread_join(handles[i],NULL);

    sem_destroy(&sem_main);
    sem_destroy(&sem_workerend);
    sem_destroy(&sem_workerstart);

    for(int i=0;i<Arr_size;i++){
        for(int j=0;j<Arr_size;j++)
            cout<<A[i][j]<<" ";
        cout<<endl;
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

    reset(Gauss_arr);
    Serial();
    reset(Gauss_arr);
    Static_thread();
    reset(Gauss_arr);




    return 0;
}
