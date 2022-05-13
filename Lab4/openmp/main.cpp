#include <iostream>
//使用openmp完成高斯消元的并行化
#include <omp.h>
#include<cstdlib>
#include<ctime>
#include<stdlib.h>
#include<sys/time.h>
#include<random>

using namespace std;

int Arr_size = 32;
float **A;
int NUM_THREADS = 4;

int cycle = 5;

void reset(float** G)
{

    for(int i=0;i<Arr_size;i++)
        for(int j=0;j<Arr_size;j++)
            A[i][j]=G[i][j];
}
void init() {
    A = new float *[Arr_size];
    for (int i = 0; i < Arr_size; i++) {
        A[i] = new float[Arr_size];
    }
    for (int i = 0; i < Arr_size; i++) {
        for (int j = 0; j < Arr_size; j++) {
            A[i][j] = 0;
        }
    }
}
void Serial() {
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

void Parallel() {
    //利用openmp对A进行并行化高斯消元
    int i,j,k;
    float temp;
    #pragma omp parallel if(Arr_size>=4) num_threads(NUM_THREADS) private(i,j,k,temp)
    {
        for (k = 0; k < Arr_size; k++) {
            #pragma omp single
            {
                temp = A[k][k];
                for (j = k + 1; j < Arr_size; j++) {
                    A[k][j] /= temp;
                }
                A[k][k] = 1.0;
            }
            #pragma omp for
            for (i = k + 1; i < Arr_size; i++) {
                temp = A[i][k];
                for (j = k + 1; j < Arr_size; j++) {
                    A[i][j] -= temp * A[k][j];
                }

                A[i][k] = 0.0;
            }
        }
    }
}


void Parallel_schedule_static_1() {
    //利用openmp对A进行并行化高斯消元,并且使用自定义的schedule
    int i,j,k;
    float temp;
    #pragma omp parallel if(Arr_size>=4) num_threads(NUM_THREADS) private(i,j,k,temp)
    {
        for (k = 0; k < Arr_size; k++) {
            #pragma omp single
            {
                temp = A[k][k];
                for (j = k + 1; j < Arr_size; j++) {
                    A[k][j] /= temp;
                }
                A[k][k] = 1.0;
            }
            #pragma omp for schedule(static, 1)
            for (i = k + 1; i < Arr_size; i++) {
                temp = A[i][k];
                for (j = k + 1; j < Arr_size; j++) {
                    A[i][j] -= temp * A[k][j];
                }

                A[i][k] = 0.0;
            }
        }
    }
}

void Parallel_schedule_static_2() {
    //利用openmp对A进行并行化高斯消元,并且使用自定义的schedule
    int i,j,k;
    float temp;

    #pragma omp parallel if(Arr_size>=4) num_threads(NUM_THREADS) private(i,j,k,temp)
    for(k=0;k<Arr_size;k++)
    {
        #pragma omp single
        {
            temp = A[k][k];
            for (j = k + 1; j < Arr_size; j++) {
                A[k][j] /= temp;
            }
            A[k][k] = 1.0;
        }

        #pragma omp for schedule(static,4)
        for(i=k+1;i<Arr_size;i++)
        {
            temp=A[i][k];
            for(j=k+1;j<Arr_size;j++)
            {
                A[i][j]-=temp*A[k][j];
            }

            A[i][k]=0.0;
        }
    }

}

void Parallel_schedule_dynamic_1(){
    int i,j,k;
    float temp;
    #pragma omp parallel if(Arr_size>=4) num_threads(NUM_THREADS) private(i,j,k,temp)
    {
        for (k = 0; k < Arr_size; k++) {
            #pragma omp single
            {
                temp = A[k][k];
                for (j = k + 1; j < Arr_size; j++) {
                    A[k][j] /= temp;
                }
                A[k][k] = 1.0;
            }
            #pragma omp for schedule(dynamic, 1)
            for (i = k + 1; i < Arr_size; i++) {
                temp = A[i][k];
                for (j = k + 1; j < Arr_size; j++) {
                    A[i][j] -= temp * A[k][j];
                }

                A[i][k] = 0.0;
            }
        }
    }
}

void Parallel_schedule_dynamic_2(){
    int i,j,k;
    float temp;
    #pragma omp parallel if(Arr_size>=4) num_threads(NUM_THREADS) private(i,j,k,temp)
    {
        for (k = 0; k < Arr_size; k++) {
            #pragma omp single
            {
                temp = A[k][k];
                for (j = k + 1; j < Arr_size; j++) {
                    A[k][j] /= temp;
                }
                A[k][k] = 1.0;
            }
            #pragma omp for schedule(dynamic, Arr_size/NUM_THREADS)
            for (i = k + 1; i < Arr_size; i++) {
                temp = A[i][k];
                for (j = k + 1; j < Arr_size; j++) {
                    A[i][j] -= temp * A[k][j];
                }

                A[i][k] = 0.0;
            }
        }
    }
}

void Parallel_schedule_guided(){
    int i,j,k;
    float temp;
    #pragma omp parallel if(Arr_size>=4) num_threads(NUM_THREADS) private(i,j,k,temp)
    {
        for (k = 0; k < Arr_size; k++) {
            #pragma omp single
            {
                temp = A[k][k];
                for (j = k + 1; j < Arr_size; j++) {
                    A[k][j] /= temp;
                }
                A[k][k] = 1.0;
            }
            #pragma omp for schedule(guided, Arr_size/NUM_THREADS)
            for (i = k + 1; i < Arr_size; i++) {
                temp = A[i][k];
                for (j = k + 1; j < Arr_size; j++) {
                    A[i][j] -= temp * A[k][j];
                }

                A[i][k] = 0.0;
            }
        }
    }
}

void Parallel_SIMD() {
    //利用openmp的自动向量化对A进行并行化高斯消元
    int i,j,k;
    float temp;
    #pragma omp parallel if(Arr_size>=4) num_threads(NUM_THREADS) private(i,j,k,temp)
    for(k=0;k<Arr_size;k++)
    {
        #pragma omp single
        {
            temp = A[k][k];
            for (j = k + 1; j < Arr_size; j++) {
                A[k][j] /= temp;
            }
            A[k][k] = 1.0;
        }
        #pragma omp for simd
        for(i=k+1;i<Arr_size;i++)
        {
            temp=A[i][k];
            for(j=k+1;j<Arr_size;j++)
            {
                A[i][j]-=temp*A[k][j];
            }

            A[i][k]=0.0;
        }
    }

}
int main() {
    static default_random_engine gengerator(1337);
    uniform_real_distribution<float> dis(-1.0, 1.0);

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
            Gauss_arr[i][j]=dis(gengerator);
    }
    for(int k=0;k<Arr_size;k++)
        for(int i=k+1;i<Arr_size;i++)
            for(int j=0;j<Arr_size;j++)
                Gauss_arr[i][j]+=Gauss_arr[k][j];

    init();
    reset(Gauss_arr);

    omp_set_num_threads(NUM_THREADS);

    //串行消元
    double time=0;
    struct timeval tv_begin,tv_end;
    for(int i=0;i<cycle;i++){
        reset(Gauss_arr);
        gettimeofday(&tv_begin,NULL);
        Serial();
        gettimeofday(&tv_end,NULL);
        time+=(tv_end.tv_sec-tv_begin.tv_sec)*1000.0+(tv_end.tv_usec-tv_begin.tv_usec)/1000.0;
    }
    cout<<"Serial time:"<<time/cycle<<"ms"<<endl;
    time=0;


    //openmp消元
    for(int i=0;i<cycle;i++){
        reset(Gauss_arr);
        gettimeofday(&tv_begin,NULL);
        Parallel();
        gettimeofday(&tv_end,NULL);
        time+=(tv_end.tv_sec-tv_begin.tv_sec)*1000.0+(tv_end.tv_usec-tv_begin.tv_usec)/1000.0;
    }
    cout<<"Parallel time:"<<time/cycle<<"ms"<<endl;
    time=0;


    //并行化高斯消元
    //openmp消元
    for(int i=0;i<cycle;i++){
        reset(Gauss_arr);
        gettimeofday(&tv_begin,NULL);
        Parallel_schedule_static_1();
        gettimeofday(&tv_end,NULL);
        time+=(tv_end.tv_sec-tv_begin.tv_sec)*1000.0+(tv_end.tv_usec-tv_begin.tv_usec)/1000.0;
    }
    cout<<"Parallel_schedule_static_1 time:"<<time/cycle<<"ms"<<endl;
    time=0;

    //并行化高斯消元,并且使用自定义的schedule
    for(int i=0;i<cycle;i++){
        reset(Gauss_arr);
        gettimeofday(&tv_begin,NULL);
        Parallel_schedule_static_2();
        gettimeofday(&tv_end,NULL);
        time+=(tv_end.tv_sec-tv_begin.tv_sec)*1000.0+(tv_end.tv_usec-tv_begin.tv_usec)/1000.0;
    }
    cout<<"Parallel_schedule_static_2 time:"<<time/cycle<<"ms"<<endl;
    time=0;

    //并行消元，自动向量化
    for(int i=0;i<cycle;i++){
        reset(Gauss_arr);
        gettimeofday(&tv_begin,NULL);
        Parallel_SIMD();
        gettimeofday(&tv_end,NULL);
        time+=(tv_end.tv_sec-tv_begin.tv_sec)*1000.0+(tv_end.tv_usec-tv_begin.tv_usec)/1000.0;
    }
    cout<<"Parallel_SIMD time:"<<time/cycle<<"ms"<<endl;
    time=0;

    for(int i=0;i<cycle;i++){
        reset(Gauss_arr);
        gettimeofday(&tv_begin,NULL);
        Parallel_schedule_dynamic_1();
        gettimeofday(&tv_end,NULL);
        time+=(tv_end.tv_sec-tv_begin.tv_sec)*1000.0+(tv_end.tv_usec-tv_begin.tv_usec)/1000.0;
    }
    cout<<"Parallel_schedule_dynamic_1 time:"<<time/cycle<<"ms"<<endl;
    time=0;

    for(int i=0;i<cycle;i++){
        reset(Gauss_arr);
        gettimeofday(&tv_begin,NULL);
        Parallel_schedule_dynamic_2();
        gettimeofday(&tv_end,NULL);
        time+=(tv_end.tv_sec-tv_begin.tv_sec)*1000.0+(tv_end.tv_usec-tv_begin.tv_usec)/1000.0;
    }
    cout<<"Parallel_schedule_dynamic_2 time:"<<time/cycle<<"ms"<<endl;
    time=0;

    for(int i=0;i<cycle;i++){
        reset(Gauss_arr);
        gettimeofday(&tv_begin,NULL);
        Parallel_schedule_guided();
        gettimeofday(&tv_end,NULL);
        time+=(tv_end.tv_sec-tv_begin.tv_sec)*1000.0+(tv_end.tv_usec-tv_begin.tv_usec)/1000.0;
    }
    cout<<"Parallel_schedule_guided time:"<<time/cycle<<"ms"<<endl;
    time=0;



    return 0;
}
