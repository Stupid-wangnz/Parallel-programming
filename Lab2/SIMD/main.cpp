#include <iostream>
#include<arm_neon.h>
#include<stdio.h>
#include<ctime>
#include<stdlib.h>
#include<cstdlib>
#include<sys/time.h>

using namespace std;
int Arr_size=16;

void Serial(float** A)
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

void Simd(float** A)
{

    for(int k=0;k<Arr_size;k++)
    {
        float32x4_t vt=vmovq_n_f32(A[k][k]);
        int j;
        for(j=k+1;j+4<=Arr_size;j+=4)
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
            for(j=k+1;j+4<=Arr_size;j+=4)
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

void Simd_Aligned(float** A)
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

void reset(float **A,float** B)
{

    for(int i=0;i<Arr_size;i++)
        for(int j=0;j<Arr_size;j++)
            B[i][j]=A[i][j];
}
void Run()
{
    cout<<"__________________________________"<<endl;
    cout<<Arr_size<<endl;

    float **Gauss_arr=new float*[Arr_size];
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

    float **Copy_arr=new float*[Arr_size];
    for(int i=0;i<Arr_size;i++)
        Copy_arr[i]=new float[Arr_size];

        float **Aligned_Gauss_arr;

    Aligned_Gauss_arr=new  float*[Arr_size];
    float* testArray[Arr_size];
    for(int i=0;i<Arr_size;i++)
       {
            Aligned_Gauss_arr[i]=(float*)aligned_alloc(128,Arr_size*4);
            testArray[i]=&Aligned_Gauss_arr[i][0];
       }

    for(int i=0;i<Arr_size;i++)
        for(int j=0;j<Arr_size;j++)
            Aligned_Gauss_arr[i][j]=Gauss_arr[i][j];

    double time=0;
    struct timeval tv_begin,tv_end;
    for(int i=0;i<20;i++){
        reset(Gauss_arr,Copy_arr);
       gettimeofday(&tv_begin,NULL);
        Serial(Copy_arr);
        gettimeofday(&tv_end,NULL);
        time+=(tv_end.tv_sec-tv_begin.tv_sec)*1000.0+(tv_end.tv_usec-tv_begin.tv_usec)/1000.0;
    }
    cout<<time/20<<endl;

    time=0;
  for(int i=0;i<20;i++){
        reset(Gauss_arr,Copy_arr);
       gettimeofday(&tv_begin,NULL);
        Simd(Copy_arr);
        gettimeofday(&tv_end,NULL);
        time+=(tv_end.tv_sec-tv_begin.tv_sec)*1000.0+(tv_end.tv_usec-tv_begin.tv_usec)/1000.0;
    }
    cout<<time/20<<endl;


    time=0;
    for(int i=0;i<20;i++){
        reset(Gauss_arr,Copy_arr);
       gettimeofday(&tv_begin,NULL);
        Simd_Aligned(Copy_arr);
        gettimeofday(&tv_end,NULL);
        time+=(tv_end.tv_sec-tv_begin.tv_sec)*1000.0+(tv_end.tv_usec-tv_begin.tv_usec)/1000.0;
    }
    cout<<time/20<<endl;

    for(int i=0;i<Arr_size;i++)
        free((void*)testArray[i]);
}
int main()
{
    for(int i=0;i<40;i++)
    {
        Run();

        Arr_size+=16;
    }
}
