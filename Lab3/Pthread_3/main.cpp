#include <iostream>
#include<pthread.h>
#include<semaphore.h>
#include<ctime>
#include<stdlib.h>
#include<cstdlib>
#include<sys/time.h>
#include<arm_neon.h>

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


void Simd_Aligned()
{
    for(int k=0;k<Arr_size;k++)
    {
        float32x4_t vt=vmovq_n_f32(A[k][k]);
        int j=k+1;

        while(j%4!=0){
            A[k][j]/=A[k][k];
            j++;
        }

        for(;j+4<=Arr_size;j+=4)
        {
            float32x4_t va=vld1q_f32(&A[k][j]);

            va=vdivq_f32(va,vt);
            vst1q_f32(&A[k][j],va);

        }
        for(j;j<Arr_size;j++)
            A[k][j]/=A[k][k];
        A[k][k]=1.0;

        for(int i=k+1;i<Arr_size;i++)
        {
            float32x4_t vaik=vmovq_n_f32(A[i][k]);
            j=k+1;
            while(j%4!=0){
                A[i][j]-=A[k][j]*A[i][k];
                j++;
            }

            for(;j+4<=Arr_size;j+=4)
            {
                float32x4_t vakj=vld1q_f32(&A[k][j]);
                float32x4_t vaij=vld1q_f32(&A[i][j]);
                float32x4_t vx=vmulq_f32(vakj,vaik);
                vaij=vsubq_f32(vaij,vx);
                vst1q_f32(&A[i][j],vaij);

            }
            for(j;j<Arr_size;j++)
                A[i][j]-=A[k][j]*A[i][k];
            A[i][k]=0;
        }
    }

}

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

void *threadFunc4(void*param){

    threadParam_t*p=(threadParam_t*)param;

    int t_id=p->t_id;
    for(int k=0;k<Arr_size;k++){
        sem_wait(&sem_workerstart[t_id]);

        for(int i=k+1+t_id;i<Arr_size;i+=NUM_THREADS){
            float32x4_t vaik=vmovq_n_f32(A[i][k]);

            int j=k+1;
            while(j%4!=0){
                A[i][j]-=A[k][j]*A[i][k];
                j++;
            }

            for(;j+4<=Arr_size;j+=4)
            {
                float32x4_t vakj=vld1q_f32(&A[k][j]);
                float32x4_t vaij=vld1q_f32(&A[i][j]);
                float32x4_t vx=vmulq_f32(vakj,vaik);
                vaij=vsubq_f32(vaij,vx);
                vst1q_f32(&A[i][j],vaij);
            }
            for(j;j<Arr_size;j++)
                A[i][j]-=A[k][j]*A[i][k];

            A[i][k]=0;

        }
        sem_post(&sem_main);
        sem_wait(&sem_workerend[t_id]);
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

void Static_thread_4()
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
        pthread_create(&handles[i],NULL,threadFunc4,&param[i]);
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

    double time=0;
    struct timeval tv_begin,tv_end;
    for(int i=0;i<3;i++){
        reset(Gauss_arr);
        gettimeofday(&tv_begin,NULL);
        Serial();
        gettimeofday(&tv_end,NULL);

        time+=(tv_end.tv_sec-tv_begin.tv_sec)*1000.0+(tv_end.tv_usec-tv_begin.tv_usec)/1000.0;
    }
    cout<<time/3<<"ms"<<endl;
    time=0;

    reset(Gauss_arr);

    for(int i=0;i<3;i++){
        reset(Gauss_arr);
        gettimeofday(&tv_begin,NULL);

        Static_thread_1();
        gettimeofday(&tv_end,NULL);

        time+=(tv_end.tv_sec-tv_begin.tv_sec)*1000.0+(tv_end.tv_usec-tv_begin.tv_usec)/1000.0;
    }
    cout<<time/3<<"ms"<<endl;
    time=0;

    reset(Gauss_arr);

    for(int i=0;i<3;i++){
        reset(Gauss_arr);
        gettimeofday(&tv_begin,NULL);
        Simd_Aligned();
        gettimeofday(&tv_end,NULL);

        time+=(tv_end.tv_sec-tv_begin.tv_sec)*1000.0+(tv_end.tv_usec-tv_begin.tv_usec)/1000.0;
    }
    cout<<time/3<<"ms"<<endl;
    time=0;

    reset(Gauss_arr);

    for(int i=0;i<3;i++){
        reset(Gauss_arr);
        gettimeofday(&tv_begin,NULL);
        Static_thread_4();
        gettimeofday(&tv_end,NULL);

        time+=(tv_end.tv_sec-tv_begin.tv_sec)*1000.0+(tv_end.tv_usec-tv_begin.tv_usec)/1000.0;
    }
    cout<<time/3<<"ms"<<endl;
    time=0;

    reset(Gauss_arr);



    return ;
}
int main(){
    for(int i=0;i<8;i++){
        cout<<Arr_size<<endl;
        Run();
        Arr_size*=2;
    }
    return 0;
}
