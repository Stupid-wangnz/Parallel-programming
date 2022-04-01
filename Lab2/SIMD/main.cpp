#include <iostream>
#include<arm_neon.h>
#include<stdio.h>
#include<time.h>

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

            va=c(va,vt);
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
        int j;
        for(j=k+1;j+4<=Arr_size;j+=4)
        {
            float32x4_t va=vld1q_f32(&A[k][j]);

            va=c(va,vt);
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
            Aligned_Gauss_arr[i]=(float*)_aligned_malloc(Arr_size*4,128);
            testArray[i]=&Aligned_Gauss_arr[i][0];
       }

    for(int i=0;i<Arr_size;i++)
        for(int j=0;j<Arr_size;j++)
            Aligned_Gauss_arr[i][j]=Gauss_arr[i][j];

    long dnesc=0£»
    struct timespec sts,ets;
    time_t dsec;
    for(int i=0;i<20;i++){
        reset(Gauss_arr,Copy_arr);
        timespec_get(&sts,TIME_UTC);
        Serial(Copy_arr);
        timespec_get(&ets,TIME_UTC);
        dsec+=ets.tv_sec-sts.tv_sec;
        dnesc+=ets.tv_nsec-sts.tv_nsec;
    }
    printf("%llu.%09llus\n",dsec,dnesc);
    cout<<endl;

   for(int i=0;i<20;i++){
        reset(Gauss_arr,Copy_arr);
        timespec_get(&sts,TIME_UTC);
        Simd(Copy_arr);
        timespec_get(&ets,TIME_UTC);
        dsec+=ets.tv_sec-sts.tv_sec;
        dnesc+=ets.tv_nsec-sts.tv_nsec;
    }
    printf("%llu.%09llus\n",dsec,dnesc);
    cout<<endl;

    for(int i=0;i<20;i++){
        reset(Aligned_Gauss_arr,Copy_arr);
        timespec_get(&sts,TIME_UTC);
        Simd_Aligned(Copy_arr);
        timespec_get(&ets,TIME_UTC);
        dsec+=ets.tv_sec-sts.tv_sec;
        dnesc+=ets.tv_nsec-sts.tv_nsec;
    }
    printf("%llu.%09llus\n",dsec,dnesc);
    cout<<endl;
}
int main()
{
    for(int i=0;i<15;i++)
    {
        Run();

        Arr_size+=16;
    }
}
