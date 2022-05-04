#include <iostream>
#include<pthread.h>
#include<semaphore.h>
#include<ctime>
#include<stdlib.h>
#include<cstdlib>
#include<sys/time.h>
using namespace std;

int Arr_size=16;
int NUM_THREADS=4;
float **A;

sem_t sem_main;
sem_t *sem_workerstart;
sem_t *sem_workerend;

//信号量
//sem_t sem_leader;
sem_t *sem_Divsion;
sem_t *sem_Elimination;

//信号量
pthread_barrier_t barrier_Divsion;
pthread_barrier_t barrier_Elimination;

void reset(float** G)
{

    for(int i=0;i<Arr_size;i++)
        for(int j=0;j<Arr_size;j++)
            A[i][j]=G[i][j];
}

typedef struct{
    int t_id;//Ïß³Ìid
}threadParam_t;

void *threadFunc1(void*param){

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
void *threadFunc2(void*param){

    threadParam_t*p=(threadParam_t*)param;

    int t_id=p->t_id;
    for(int k=0;k<Arr_size;k++){
        int workCount=(Arr_size-k-1)/NUM_THREADS;
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

        for(int i=k+1+t_id*workCount;i<k+1+(t_id+1)*workCount;i++){
            for(int j=k+1;j<Arr_size;j++)
                A[i][j]=A[i][j]-A[i][k]*A[k][j];

            A[i][k]=0;
        }

        if(t_id==0){
            for(int i=k+1+NUM_THREADS*workCount;i<Arr_size;i++){
                for(int j=k+1;j<Arr_size;j++)
                    A[i][j]-=A[i][k]*A[k][j];
                A[i][k]=0;
            }

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

void *threadFunc3(void*param){

    threadParam_t*p=(threadParam_t*)param;

    int t_id=p->t_id;
    for(int k=0;k<Arr_size;k++){
    if(k>0&&t_id==0){
        for(int i=k;i<Arr_size;i++)
            A[i][k-1]=0;
    }
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

        for(int i=k+1;i<Arr_size;i++){
            for(int j=k+1+t_id;j<Arr_size;j+=NUM_THREADS)
                A[i][j]=A[i][j]-A[i][k]*A[k][j];

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
void *threadFunc4(void*param){
    double time=0;
    struct timeval tv_begin,tv_end;
    gettimeofday(&tv_begin,NULL);
    threadParam_t*p=(threadParam_t*)param;

    int t_id=p->t_id;
    for(int k=0;k<Arr_size;k++){
        if(k>0&&t_id==0){
        for(int i=k;i<Arr_size;i++)
            A[i][k-1]=0;
    }
        int workCount=(Arr_size-k-1)/NUM_THREADS;
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

        for(int i=k+1;i<Arr_size;i++){
            for(int j=k+1+t_id*workCount;j<k+1+(t_id+1)*workCount;j++)
                A[i][j]=A[i][j]-A[i][k]*A[k][j];

        }

        if(t_id==0){
            for(int i=k+1;i<Arr_size;i++){
                for(int j=k+1+NUM_THREADS*workCount;j<Arr_size;j++)
                    A[i][j]-=A[i][k]*A[k][j];
            }

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
void Static_thread_1()
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
        pthread_create(&handles[i],NULL,threadFunc1,&param[i]);
    }

    for(int i=0;i<NUM_THREADS;i++)
        pthread_join(handles[i],NULL);
    for(int i=0;i<Arr_size-1;i++){
        sem_destroy(&sem_Divsion[i]);
        sem_destroy(&sem_Elimination[i]);
    }

}
void Static_thread_2()
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
        pthread_create(&handles[i],NULL,threadFunc2,&param[i]);
    }

    for(int i=0;i<NUM_THREADS;i++)
        pthread_join(handles[i],NULL);
    for(int i=0;i<Arr_size-1;i++){
        sem_destroy(&sem_Divsion[i]);
        sem_destroy(&sem_Elimination[i]);
    }
}
void Static_thread_3()
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
        pthread_create(&handles[i],NULL,threadFunc3,&param[i]);
    }

    for(int i=0;i<NUM_THREADS;i++)
        pthread_join(handles[i],NULL);
    for(int i=0;i<Arr_size-1;i++){
        sem_destroy(&sem_Divsion[i]);
        sem_destroy(&sem_Elimination[i]);
    }
}
void Static_thread_4()
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
        pthread_create(&handles[i],NULL,threadFunc4,&param[i]);
    }

    for(int i=0;i<NUM_THREADS;i++)
        pthread_join(handles[i],NULL);
    for(int i=0;i<Arr_size-1;i++){
        sem_destroy(&sem_Divsion[i]);
        sem_destroy(&sem_Elimination[i]);
    }
}
void Run()
{
    srand(Arr_size);

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
            Gauss_arr[i][j]=rand()%10+1;


    }
    reset(Gauss_arr);

    for(int k=0;k<Arr_size;k++)
        for(int i=k+1;i<Arr_size;i++)
            for(int j=0;j<Arr_size;j++)
                Gauss_arr[i][j]+=Gauss_arr[k][j];
    reset(Gauss_arr);

    double time=0;
    struct timeval tv_begin,tv_end;
//    for(int i=0;i<3;i++){
//        reset(Gauss_arr);
//        gettimeofday(&tv_begin,NULL);
//        Serial();
//        gettimeofday(&tv_end,NULL);
//
//        time+=(tv_end.tv_sec-tv_begin.tv_sec)*1000.0+(tv_end.tv_usec-tv_begin.tv_usec)/1000.0;
//    }
//    cout<<time/3<<"ms"<<endl;
//    time=0;

    for(int i=0;i<3;i++){
        reset(Gauss_arr);
        gettimeofday(&tv_begin,NULL);

        Static_thread_1();
        gettimeofday(&tv_end,NULL);

        time+=(tv_end.tv_sec-tv_begin.tv_sec)*1000.0+(tv_end.tv_usec-tv_begin.tv_usec)/1000.0;
    }
    cout<<time/3<<"ms"<<endl;
    time=0;

//    for(int i=0;i<3;i++){
//        reset(Gauss_arr);
//        gettimeofday(&tv_begin,NULL);
//        Static_thread_2();
//        gettimeofday(&tv_end,NULL);
//
//        time+=(tv_end.tv_sec-tv_begin.tv_sec)*1000.0+(tv_end.tv_usec-tv_begin.tv_usec)/1000.0;
//    }
//    cout<<time/3<<"ms"<<endl;
//    time=0;
//
//    reset(Gauss_arr);

//
//    for(int i=0;i<3;i++){
//        reset(Gauss_arr);
//        gettimeofday(&tv_begin,NULL);
//        Static_thread_3();
//        gettimeofday(&tv_end,NULL);
//
//        time+=(tv_end.tv_sec-tv_begin.tv_sec)*1000.0+(tv_end.tv_usec-tv_begin.tv_usec)/1000.0;
//    }
//    cout<<time/3<<"ms"<<endl;
//    time=0;
//    reset(Gauss_arr);

//    for(int i=0;i<3;i++){
//        reset(Gauss_arr);
//        gettimeofday(&tv_begin,NULL);
//        Static_thread_4();
//        gettimeofday(&tv_end,NULL);
//
//        time+=(tv_end.tv_sec-tv_begin.tv_sec)*1000.0+(tv_end.tv_usec-tv_begin.tv_usec)/1000.0;
//    }
//    cout<<time/3<<"ms"<<endl;
//    time=0;
//
//    reset(Gauss_arr);

    return ;
}
int main(){

    Arr_size=1024;
    Run();

    return 0;
}
