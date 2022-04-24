#include <iostream>
#include<pthread.h>
#include<semaphore.h>

using namespace std;

int Arr_size=16;
int NUM_THREADS=4;
float **A;

//信号量
//sem_t sem_leader;
sem_t *sem_Divsion;
sem_t *sem_Elimination;

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

    //t_id为0时，线程做除法，让其余线程阻塞
        if(t_id==0){
            for(int j=k+1;j<Arr_size;j++)
                A[k][j]/=A[k][k];
            A[k][k]=1;

         for(int i=0;i<NUM_THREADS-1;i++)
                sem_post(&sem_Divsion[i]);
        }
        else{
            sem_wait(&sem_Divsion[t_id-1]);
        }

        for(int i=k+1+t_id;i<Arr_size;i+=NUM_THREADS){
            for(int j=k+1;j<Arr_size;j++)
                A[i][j]=A[i][j]-A[i][k]*A[k][j];

            A[i][k]=0;
        }

        if(t_id==0){
            for(int i=0;i<NUM_THREADS-1;i++)
                sem_wait(&sem_Elimination[i]);
        }
        else{
            sem_post(&sem_Elimination[t_id-1]);
        }

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

void Static_thread()
{
    sem_Divsion=new sem_t[Arr_size-1];
    sem_Elimination=new sem_t[Arr_size-1];
    for(int i=0;i<Arr_size-1;i++){
        sem_init(&sem_Divsion[i],0,0);
        sem_init(&sem_Elimination[i],0,0);
    }
    pthread_t*handles=new pthread_t[NUM_THREADS];
    threadParam_t*param=new threadParam_t[NUM_THREADS];
    for(int i=0;i<NUM_THREADS;i++){
        param[i].t_id=i;
        pthread_create(&handles[i],NULL,threadFunc,&param[i]);
    }

    for(int i=0;i<NUM_THREADS;i++)
        pthread_join(handles[i],NULL);
    for(int i=0;i<Arr_size-1;i++){
        sem_destroy(&sem_Divsion[i]);
        sem_destroy(&sem_Elimination[i]);
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
    for(int i=0;i<Arr_size;i++){
        for(int j=0;j<Arr_size;j++)
            cout<<A[i][j]<<" ";
        cout<<endl;
    }
    reset(Gauss_arr);

    Static_thread();

    for(int i=0;i<Arr_size;i++){
        for(int j=0;j<Arr_size;j++)
            cout<<A[i][j]<<" ";
        cout<<endl;
    }


    return 0;
}
