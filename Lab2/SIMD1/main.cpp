#include <iostream>
#include<stdio.h>
#include<pmmintrin.h>
#include<stdlib.h>
#include<algorithm>
#include<windows.h>
#include<time.h>
#include<malloc.h>

using namespace std;


const int Arr_size=160;

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

void Simd_SSE(float A[Arr_size][Arr_size])
{

    for(int k=0;k<Arr_size;k++)
    {

        __m128 vt;
        float temp[4]={A[k][k],A[k][k],A[k][k],A[k][k]};
        vt=_mm_loadu_ps(temp);
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
            float temp1[4]={A[i][k],A[i][k],A[i][k],A[i][k]};
            __m128 vaik=_mm_loadu_ps(temp1);
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

void Simd_SSE_Aligned(float A[Arr_size][Arr_size])
{
    /*
    对齐下的SIMD SSE指令操作
    */
    for(int k=0;k<Arr_size;k++)
    {
        __m128 vt;
        float temp[4]__attribute__((aligned(16)))={A[k][k],A[k][k],A[k][k],A[k][k]};
        vt=_mm_load_ps(temp);
        int j;
        for(j=k+1;j+4<=Arr_size;j+=4)
        {
            float temp1[4]__attribute__((aligned(16)))={A[k][j],A[k][j+1],A[k][j+2],A[k][j+3]};
            __m128 va=_mm_load_ps(temp1);
            va=_mm_div_ps(va,vt);
            _mm_storeu_ps(&A[k][j],va);
        }
        for(j;j<Arr_size;j++)
            A[k][j]/=A[k][k];
        A[k][k]=1.0;
        for(int i=k+1;i<Arr_size;i++)
        {
            float temp1[4]__attribute__((aligned(16)))={A[i][k],A[i][k],A[i][k],A[i][k]};
            __m128 vaik=_mm_load_ps(temp1);
            for(j=k+1;j+4<=Arr_size;j+=4)
            {
                float tempkj[4]__attribute__((aligned(16)))={A[k][j],A[k][j+1],A[k][j+2],A[k][j+3]};
                __m128 vakj=_mm_load_ps(tempkj);
                float tempij[4]__attribute__((aligned(16)))={A[i][j],A[i][j+1],A[i][j+2],A[i][j+3]};
                __m128 vaij=_mm_load_ps(tempij);
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

void Simd_SSE_Aligned1(float A[Arr_size][Arr_size])
{

    for(int k=0;k<Arr_size;k++)
    {

        __m128 vt;
        float temp[4]__attribute__((aligned(16)))={A[k][k],A[k][k],A[k][k],A[k][k]};
        vt=_mm_load_ps(temp);
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
            float temp1[4]__attribute__((aligned(16)))={A[i][k],A[i][k],A[i][k],A[i][k]};
            __m128 vaik=_mm_load_ps(temp1);
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

    float **Aligned_Gauss_arr;

//    Aligned_Gauss_arr=new  float*[Arr_size];
//    float* testArray[Arr_size];
//    for(int i=0;i<Arr_size;i++)
//       {
//            Aligned_Gauss_arr[i]=(float*)_aligned_malloc(Arr_size*4,16);
//            testArray[i]=&Aligned_Gauss_arr[i][0];
//       }
//
//
//
//    for(int i=0;i<Arr_size;i++){
//        cout<<testArray[i]<<endl;
//        _aligned_free((void*)testArray[i]);
//    }

    float Copy_arr[Arr_size][Arr_size]{};
    reset(Gauss_arr,Copy_arr);
    long long head,tail,freq;
    double time=0;

    //串行高斯消元
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    QueryPerformanceCounter((LARGE_INTEGER*)&head);
    Serial(Copy_arr);
    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
    time=(tail-head)*1000.0/freq;
    cout<<time<<endl;

    //不对齐的高斯消元
    reset(Gauss_arr,Copy_arr);
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    QueryPerformanceCounter((LARGE_INTEGER*)&head);
    Simd_SSE(Copy_arr);
    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
    time=(tail-head)*1000.0/freq;
    cout<<time<<endl;

    //将数据全部对齐后的高斯消元
    reset(Gauss_arr,Copy_arr);
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    QueryPerformanceCounter((LARGE_INTEGER*)&head);
    Simd_SSE_Aligned(Copy_arr);
    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
    time=(tail-head)*1000.0/freq;
    cout<<time<<endl;

    //部分数据对齐后的高斯消元
    reset(Gauss_arr,Copy_arr);
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    QueryPerformanceCounter((LARGE_INTEGER*)&head);
    Simd_SSE_Aligned1(Copy_arr);
    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
    time=(tail-head)*1000.0/freq;
    cout<<time<<endl;
}
