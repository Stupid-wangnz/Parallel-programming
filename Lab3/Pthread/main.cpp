#include <iostream>
#include<pthread.h>


using namespace std;


int Arr_size;
float **A;

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

void Dynamic_thread(float**A)
{


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




    return 0;
}
