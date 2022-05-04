#include <iostream>
#include<pthread.h>
#include<semaphore.h>
#include<ctime>
#include<cstdlib>
#include<windows.h>
#include<time.h>
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

void Simd_AVX()
{
    for(int k=0;k<Arr_size;k++)
    {

        __m256 vt=_mm256_set1_ps(A[k][k]);
        int j;
        for(j=k+1;j+8<=Arr_size;j+=8)
        {
            __m256 va=_mm256_loadu_ps(&A[k][j]);
            va=_mm256_div_ps(va,vt);
            _mm256_storeu_ps(&A[k][j],va);

        }
        for(j;j<Arr_size;j++)
            A[k][j]/=A[k][k];
        A[k][k]=1.0;

        for(int i=k+1;i<Arr_size;i++)
        {
            __m256 vaik=_mm256_set1_ps(A[i][k]);
            for(j=k+1;j+8<=Arr_size;j+=8)
            {
                __m256 vakj=_mm256_loadu_ps(&A[k][j]);
                __m256 vaij=_mm256_loadu_ps(&A[i][j]);
                __m256 vx=_mm256_mul_ps(vakj,vaik);
                vaij=_mm256_sub_ps(vaij,vx);
                _mm256_storeu_ps(&A[i][j],vaij);

            }
            for(j;j<Arr_size;j++)
                A[i][j]-=A[k][j]*A[i][k];
            A[i][k]=0;
        }
    }

}

void Simd_SSE_Aligned1()
{

    for(int k=0;k<Arr_size;k++)
    {
        __m128 vt;
        //float temp[4]__attribute__((aligned(16)))={A[k][k],A[k][k],A[k][k],A[k][k]};
        vt=_mm_set1_ps(A[k][k]);
        int j;
        for(j=k+1;j+4<=Arr_size;j+=4)
        {
            __m128 va=_mm_loadu_ps(&A[k][j]);
            va=_mm_div_ps(va,vt);
            _mm_storeu_ps(&A[k][j],va);

        }
        for(j;j<Arr_size;j++)
            A[k][j]/=A[k][k];
        A[k][k]=1.0;

        for(int i=k+1;i<Arr_size;i++)
        {
            __m128 vaik=_mm_set1_ps(A[i][k]);
            for(j=k+1;j+4<=Arr_size;j+=4)
            {
                __m128 vakj=_mm_loadu_ps(&A[k][j]);
                __m128 vaij=_mm_loadu_ps(&A[i][j]);
                __m128 vx=_mm_mul_ps(vakj,vaik);
                vaij=_mm_sub_ps(vaij,vx);
                _mm_storeu_ps(&A[i][j],vaij);

            }
            for(j;j<Arr_size;j++)
                A[i][j]-=A[k][j]*A[i][k];
            A[i][k]=0;
        }
    }
}


void *threadFunc1(void*param){
    threadParam_t*p=(threadParam_t*)param;

    int t_id=p->t_id;
    for(int k=0;k<Arr_size;k++){
        sem_wait(&sem_workerstart[t_id]);

        for(int i=k+1+t_id;i<Arr_size;i+=NUM_THREADS){
            for(int j=k+1;j<Arr_size;j++)
                A[i][j]=A[i][j]-A[i][k]*A[k][j];

            A[i][k]=0;
        }
        sem_post(&sem_main);
        sem_wait(&sem_workerend[t_id]);
    }
    pthread_exit(NULL);

}

void *threadFunc2(void*param){
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

void *threadFunc3(void*param){
    threadParam_t*p=(threadParam_t*)param;

    int t_id=p->t_id;

    for(int k=0;k<Arr_size;k++){

        //t_id为0时，线程做除法，让其余线程阻塞
        if(t_id==0){
            for(int j=k+1;j<Arr_size;j++)
                A[k][j]/=A[k][k];
            A[k][k]=1;

        }

        pthread_barrier_wait(&barrier_Divsion);

        for(int i=k+1+t_id;i<Arr_size;i+=NUM_THREADS){
            for(int j=k+1;j<Arr_size;j++)
                A[i][j]=A[i][j]-A[i][k]*A[k][j];

            A[i][k]=0;
        }

        pthread_barrier_wait(&barrier_Elimination);

    }
    pthread_exit(NULL);
}

void *threadFunc4(void*param){
    threadParam_t*p=(threadParam_t*)param;

    int t_id=p->t_id;
    for(int k=0;k<Arr_size;k++){
        int j;
        //t_id为0时，线程做除法，让其余线程阻塞
        if(t_id==0){
            __m128 vt=_mm_set1_ps(A[k][k]);
            for(j=k+1;j+4<=Arr_size;j+=4)
            {
                __m128 va=_mm_loadu_ps(&A[k][j]);
                va=_mm_div_ps(va,vt);
                _mm_storeu_ps(&A[k][j],va);

            }
            for(j;j<Arr_size;j++)
                A[k][j]/=A[k][k];
            A[k][k]=1.0;

            for(int i=0;i<NUM_THREADS-1;i++)
                sem_post(&sem_Divsion[i]);
        }
        else{
            sem_wait(&sem_Divsion[t_id-1]);
        }

        for(int i=k+1+t_id;i<Arr_size;i+=NUM_THREADS){
            __m128 vaik=_mm_set1_ps(A[i][k]);
            for(j=k+1;j+4<=Arr_size;j+=4)
            {
                __m128 vakj=_mm_loadu_ps(&A[k][j]);
                __m128 vaij=_mm_loadu_ps(&A[i][j]);
                __m128 vx=_mm_mul_ps(vakj,vaik);
                vaij=_mm_sub_ps(vaij,vx);
                _mm_storeu_ps(&A[i][j],vaij);

            }
            for(j;j<Arr_size;j++)
                A[i][j]-=A[k][j]*A[i][k];
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

void *threadFunc5(void*param){
    threadParam_t*p=(threadParam_t*)param;

    int t_id=p->t_id;
    for(int k=0;k<Arr_size;k++){
        int j;
        //t_id为0时，线程做除法，让其余线程阻塞
        if(t_id==0){
            __m256 vt=_mm256_set1_ps(A[k][k]);
            for(j=k+1;j+8<=Arr_size;j+=8)
            {
                __m256 va=_mm256_loadu_ps(&A[k][j]);
                va=_mm256_div_ps(va,vt);
                _mm256_storeu_ps(&A[k][j],va);

            }
            for(j;j<Arr_size;j++)
                A[k][j]/=A[k][k];
            A[k][k]=1.0;

            for(int i=0;i<NUM_THREADS-1;i++)
                sem_post(&sem_Divsion[i]);
        }
        else{
            sem_wait(&sem_Divsion[t_id-1]);
        }

        for(int i=k+1+t_id;i<Arr_size;i+=NUM_THREADS){
            __m256 vaik=_mm256_set1_ps(A[i][k]);
            for(j=k+1;j+8<=Arr_size;j+=8)
            {
                __m256 vakj=_mm256_loadu_ps(&A[k][j]);
                __m256 vaij=_mm256_loadu_ps(&A[i][j]);
                __m256 vx=_mm256_mul_ps(vakj,vaik);
                vaij=_mm256_sub_ps(vaij,vx);
                _mm256_storeu_ps(&A[i][j],vaij);

            }
            for(j;j<Arr_size;j++)
                A[i][j]-=A[k][j]*A[i][k];
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

void Static_thread_1()
{
    sem_workerend=new sem_t[NUM_THREADS];
    sem_workerstart=new sem_t[NUM_THREADS];

    sem_init(&sem_main,0,0);
    for(int i=0;i<NUM_THREADS;i++){
        sem_init(&sem_workerstart[i],0,0);
        sem_init(&sem_workerend[i],0,0);
    }
    pthread_t*handles=new pthread_t[NUM_THREADS];
    threadParam_t*param=new threadParam_t[NUM_THREADS];
    for(int i=0;i<NUM_THREADS;i++){
        param[i].t_id=i;
        pthread_create(&handles[i],NULL,threadFunc1,&param[i]);
    }
    for(int k=0;k<Arr_size;k++){
        for(int j=k+1;j<Arr_size;j++)
            A[k][j]/=A[k][k];

        A[k][k]=1;

        for(int i=0;i<NUM_THREADS;i++)
            sem_post(&sem_workerstart[i]);

        for(int i=0;i<NUM_THREADS;i++)
            sem_wait(&sem_main);

        for(int i=0;i<NUM_THREADS;i++)
            sem_post(&sem_workerend[i]);
    }

    for(int i=0;i<NUM_THREADS;i++)
        pthread_join(handles[i],NULL);

    sem_destroy(&sem_main);
    for(int i=0;i<NUM_THREADS;i++){
        sem_destroy(&sem_workerend[i]);
        sem_destroy(&sem_workerstart[i]);
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
    pthread_barrier_init(&barrier_Divsion,NULL,NUM_THREADS);
    pthread_barrier_init(&barrier_Elimination,NULL,NUM_THREADS);

    pthread_t*handles=new pthread_t[NUM_THREADS];
    threadParam_t*param=new threadParam_t[NUM_THREADS];
    for(int i=0;i<NUM_THREADS;i++){
        param[i].t_id=i;
        pthread_create(&handles[i],NULL,threadFunc3,&param[i]);
    }

    for(int i=0;i<NUM_THREADS;i++)
        pthread_join(handles[i],NULL);

    pthread_barrier_destroy(&barrier_Divsion);
    pthread_barrier_destroy(&barrier_Elimination);
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

void Static_thread_5()
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
        pthread_create(&handles[i],NULL,threadFunc5,&param[i]);
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

    long long head,tail,freq;
    double time=0;

    for(int i=0;i<1;i++){
        reset(Gauss_arr);
        QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
        QueryPerformanceCounter((LARGE_INTEGER*)&head);
        Serial();
        QueryPerformanceCounter((LARGE_INTEGER*)&tail);
        time+=(tail-head)*1000.0/freq;
    }
    cout<<time/1<<"ms"<<endl;
    time=0;

    for(int i=0;i<1;i++){
        reset(Gauss_arr);
        QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
        QueryPerformanceCounter((LARGE_INTEGER*)&head);
        Simd_SSE_Aligned1();
        QueryPerformanceCounter((LARGE_INTEGER*)&tail);
        time+=(tail-head)*1000.0/freq;
    }
    cout<<time/1<<"ms"<<endl;
    time=0;

    for(int i=0;i<1;i++){
        reset(Gauss_arr);
        QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
        QueryPerformanceCounter((LARGE_INTEGER*)&head);
        Static_thread_4();
        QueryPerformanceCounter((LARGE_INTEGER*)&tail);
        time+=(tail-head)*1000.0/freq;
    }
    cout<<time/1<<"ms"<<endl;
    time=0;

    for(int i=0;i<1;i++){
        reset(Gauss_arr);
        QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
        QueryPerformanceCounter((LARGE_INTEGER*)&head);
        Simd_AVX();
        QueryPerformanceCounter((LARGE_INTEGER*)&tail);
        time+=(tail-head)*1000.0/freq;
    }
    cout<<time/1<<"ms"<<endl;
    time=0;

    for(int i=0;i<1;i++){
        reset(Gauss_arr);
        QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
        QueryPerformanceCounter((LARGE_INTEGER*)&head);
        Static_thread_5();
        QueryPerformanceCounter((LARGE_INTEGER*)&tail);
        time+=(tail-head)*1000.0/freq;
    }
    cout<<time/1<<"ms"<<endl;
    time=0;

    return ;
}
int main(){




    Arr_size=3072;
    Run();

    Arr_size=4096;
    Run();

    return 0;
}
