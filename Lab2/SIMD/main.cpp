#include <iostream>
#include<arm_neon.h>
#include<stdio.h>
#include<time.h>

using namespace std;
const int Arr_size=16;

void Serial(float A[Arr_size][Arr_size])
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

void Simd(float A[Arr_size][Arr_size])
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
void reset(float A[Arr_size][Arr_size],float B[Arr_size][Arr_size])
{

    for(int i=0;i<Arr_size;i++)
        for(int j=0;j<Arr_size;j++)
            B[i][j]=A[i][j];
}
int main()
{
    float Gauss_arr[Arr_size][Arr_size]{};
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

    float Copy_arr[Arr_size][Arr_size]{};
    reset(Gauss_arr,Copy_arr);
    struct timespec sts,ets;
    timespec_get(&sts,TIME_UTC);
    Serial(Copy_arr);
    timespec_get(&ets,TIME_UTC);
    time_t dsec=ets.tv_sec-sts.tv_sec;
    long dnesc =ets.tv_nsec-sts.tv_nsec;
    printf("%llu.%09llus\n",dsec,dnesc);

    timespec_get(&sts,TIME_UTC);
    Simd(Gauss_arr);
    timespec_get(&ets,TIME_UTC);
    dsec=ets.tv_sec-sts.tv_sec;
    dnesc =ets.tv_nsec-sts.tv_nsec;
    printf("%llu.%09llus\n",dsec,dnesc);

}
