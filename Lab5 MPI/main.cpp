#include <iostream>
#include <cstdlib>
#include <random>
#include <ctime>
#include <sys/time.h>
#include "mpi.h"
#include <arm_neon.h>
#include <omp.h>

using namespace std;

static const int Arr_size = 32;
float A[Arr_size][Arr_size];
float Gauss_arr[Arr_size][Arr_size];
int cycle = 4;
int thread_count = 4;

int task_per_step=1;

void reset(float A[][Arr_size],float G[][Arr_size])
{
    for(int i=0;i<Arr_size;i++)
        for(int j=0;j<Arr_size;j++)
            A[i][j]=G[i][j];
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

void eliminate(float G[][Arr_size], int my_id,int num_proc) {
    int block_size = Arr_size / num_proc;
    int remain=Arr_size%num_proc;

    int start_row = my_id * block_size;
    int end_row = start_row + block_size;
    //if my_id is the last process, then add the remain rows
    if(my_id==num_proc-1)
        end_row+=remain;

    for(int k=0;k<Arr_size;k++)
    {
        if(k>=start_row && k<end_row)
        {
            for(int j=k+1;j<Arr_size;j++)
            {
                G[k][j]/=G[k][k];
            }
            G[k][k]=1.0;
            for(int p=my_id+1;p<num_proc;p++)
            {
                MPI_Send(&G[k],Arr_size,MPI_FLOAT,p,2,MPI_COMM_WORLD);
            }
        }
        else
        {
            int cur_p=k/block_size;
            if(cur_p<my_id)
            {
                MPI_Recv(&G[k],Arr_size,MPI_FLOAT,cur_p,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
        }
        for(int i= start_row;i<end_row;i++)
        {
            if(i>=k+1) {
                for (int j = k + 1; j < Arr_size; j++) {
                    G[i][j] -= G[i][k] * G[k][j];
                }
                G[i][k] = 0.0;
            }
        }
    }
}

void eliminate_SIMD_OPENMP(float G[][Arr_size], int my_id,int num_proc) {
    int block_size = Arr_size / num_proc;
    int remain = Arr_size % num_proc;

    int start_row = my_id * block_size;
    int end_row = start_row + block_size;
    //if my_id is the last process, then add the remain rows
    if (my_id == num_proc - 1)
        end_row += remain;
    int i, j, k;
    float32x4_t vt, va, vaik, vakj, vaij, vx;
#pragma omp parallel if(Arr_size>=4) num_threads(thread_count)  private(i, j, k, vt, va, vaik, vakj, vaij, vx)
    for (k = 0; k < Arr_size; k++) {
#pragma omp single
        {
            if (k >= start_row && k < end_row) {
                {
                    vt = vmovq_n_f32(G[k][k]);
                    for (j = k + 1; j <= Arr_size - 4; j += 4) {
                        va = vld1q_f32(&G[k][j]);
                        va = vdivq_f32(va, vt);
                        vst1q_f32(&G[k][j], va);
                    }
                    for (j; j < Arr_size; j++)
                        G[k][j] /= G[k][k];
                    G[k][k] = 1.0;

                    for (int p = my_id + 1; p < num_proc; p++) {
                        MPI_Send(&G[k], Arr_size, MPI_FLOAT, p, 2, MPI_COMM_WORLD);
                    }
                }
            } else {
                int cur_p = k / block_size;
                if (cur_p < my_id)
                    MPI_Recv(&G[k], Arr_size, MPI_FLOAT, cur_p, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
#pragma omp for
        for (i = start_row; i < end_row; i++) {
            if (i >= k + 1) {
                vaik = vmovq_n_f32(G[i][k]);
                j = k + 1;
                for (; j + 4 <= Arr_size; j += 4) {
                    vakj = vld1q_f32(&G[k][j]);
                    vaij = vld1q_f32(&G[i][j]);
                    vx = vmulq_f32(vakj, vaik);
                    vaij = vsubq_f32(vaij, vx);
                    vst1q_f32(&G[i][j], vaij);
                }
                for (; j < Arr_size; j++)
                    G[i][j] -= G[k][j] * G[i][k];
                G[i][k] = 0;
            }
        }
    }
}

void eliminate_cycle(float G[][Arr_size], int my_id,int num_proc){
    int steps=num_proc;
    for(int k=0;k<Arr_size;k++){
        if(k%steps==my_id)
        {
            for(int j=k+1;j<Arr_size;j++)
            {
                G[k][j]/=G[k][k];
            }
            G[k][k]=1.0;
            for(int p=0;p<num_proc;p++)
                if(p!=my_id) {
                    MPI_Send(&G[k], Arr_size, MPI_FLOAT, p, 2, MPI_COMM_WORLD);
                }
        }
        else
        {
            int cur_p=k%steps;
            MPI_Recv(&G[k],Arr_size,MPI_FLOAT,cur_p,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        for(int i=k+1;i<Arr_size;i++)
        {
            if(i%steps==my_id)
            {
                for(int j=k+1;j<Arr_size;j++)
                {
                    G[i][j]-=G[i][k]*G[k][j];
                }
                G[i][k]=0.0;
            }
        }
    }
}

void eliminate_cycle_task(float G[][Arr_size], int my_id,int num_proc){
    int steps=num_proc*task_per_step;
    for(int k=0;k<Arr_size;k++){
        if((k%steps)/task_per_step==my_id)
        {
            for(int j=k+1;j<Arr_size;j++)
            {
                G[k][j]/=G[k][k];
            }
            G[k][k]=1.0;

            for(int p=0;p<num_proc;p++)
                if(p!=my_id) {
                    MPI_Send(&G[k], Arr_size, MPI_FLOAT, p, 2, MPI_COMM_WORLD);
                }
        }
        else
        {
            int cur_p=(k%steps)/task_per_step;
            MPI_Recv(&G[k],Arr_size,MPI_FLOAT,cur_p,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        for(int i=k+1;i<Arr_size;i++)
        {
            if((i%steps)/task_per_step==my_id)
            {
                for(int j=k+1;j<Arr_size;j++)
                {
                    G[i][j]-=G[i][k]*G[k][j];
                }
                G[i][k]=0.0;
            }
        }
    }
}

void eliminate_cycle_SIMD_OPENMP(float G[][Arr_size], int my_id,int num_proc) {
    int steps = num_proc;
    int i, j, k;
    float32x4_t vt, va, vaik, vakj, vaij, vx;
    #pragma omp parallel if(Arr_size>=4) num_threads(thread_count)  private(i, j, k, vt, va, vaik, vakj, vaij, vx)
    for (k = 0; k < Arr_size; k++) {
        #pragma omp single
        {
            if (k % steps == my_id)
            {
                vt = vmovq_n_f32(G[k][k]);
                for (j = k + 1; j < Arr_size - 4; j += 4) {
                    va = vld1q_f32(&G[k][j]);
                    va = vdivq_f32(va, vt);
                    vst1q_f32(&G[k][j], va);
                }
                for (; j < Arr_size; j++)
                    G[k][j] /= G[k][k];
                G[k][k] = 1.0;

                for (int p = 0; p < num_proc; p++) {
                    if (p != my_id)
                        MPI_Send(&G[k], Arr_size, MPI_FLOAT, p, 2, MPI_COMM_WORLD);
                }
            }
            else
            {
                int cur_p = k % steps;
                MPI_Recv(&G[k], Arr_size, MPI_FLOAT, cur_p, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        #pragma omp for
        for (i = k + 1; i < Arr_size; i++) {
            if (i % steps == my_id) {
                vaik = vmovq_n_f32(G[i][k]);
                j = k + 1;
                for (; j + 4 <= Arr_size; j += 4) {
                    vakj = vld1q_f32(&G[k][j]);
                    vaij = vld1q_f32(&G[i][j]);
                    vx = vmulq_f32(vakj, vaik);
                    vaij = vsubq_f32(vaij, vx);
                    vst1q_f32(&G[i][j], vaij);
                }
                for (; j < Arr_size; j++)
                    G[i][j] -= G[i][k] * G[k][j];
                G[i][k] = 0.0;
            }
        }
    }
}

void eliminate_pipeline(float G[][Arr_size],int my_id,int num_proc){
    int steps=num_proc*task_per_step;
    int next_proc=(my_id+1)%num_proc;
    int pre_proc=(my_id+(num_proc-1))%num_proc;
    for(int k=0;k<Arr_size;k++){
        if((k%steps)/task_per_step==my_id)
        {
            for(int j=k+1;j<Arr_size;j++)
            {
                G[k][j]/=G[k][k];
            }
            G[k][k]=1.0;

            MPI_Send(&G[k],Arr_size,MPI_FLOAT,next_proc,2,MPI_COMM_WORLD);
        }
        else
            MPI_Recv(&G[k],Arr_size,MPI_FLOAT,pre_proc,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        for(int i=k+1;i<Arr_size;i++)
        {
            if((i%steps)/task_per_step==my_id)
            {
                for(int j=k+1;j<Arr_size;j++)
                {
                    G[i][j]-=G[i][k]*G[k][j];
                }
                G[i][k]=0.0;
            }
        }
    }
}



void MPI_Parallel() {
    timeval t_start, t_end;

    int my_id,num_proc;

    MPI_Comm_rank(MPI_COMM_WORLD,&my_id);
    MPI_Comm_size(MPI_COMM_WORLD,&num_proc);

    int block_size = Arr_size / num_proc;
    int remain=Arr_size%num_proc;

    if(my_id==0){
        reset(A,Gauss_arr);
        gettimeofday(&t_start, NULL);
        for(int i=1;i<num_proc;i++)
        {
            if(i!=num_proc-1){
                for(int j=0;j<block_size;j++)
                    MPI_Send(&A[i*block_size+j],Arr_size,MPI_FLOAT,i,0,MPI_COMM_WORLD);
            }
            else{
                for(int j=0;j<block_size+remain;j++)
                    MPI_Send(&A[i*block_size+j],Arr_size,MPI_FLOAT,i,0,MPI_COMM_WORLD);
            }
        }
        eliminate(A,my_id,num_proc);
        for(int i=1;i<num_proc;i++)
        {
            if(i!=num_proc-1){
                for(int j=0;j<block_size;j++)
                    MPI_Recv(&A[i*block_size+j],Arr_size,MPI_FLOAT,i,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
            else{
                for(int j=0;j<block_size+remain;j++)
                    MPI_Recv(&A[i*block_size+j],Arr_size,MPI_FLOAT,i,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
        }
        gettimeofday(&t_end, NULL);
        cout<<"MPI_BLOCK time: "<<(t_end.tv_sec-t_start.tv_sec)*1000+(t_end.tv_usec-t_start.tv_usec)*0.001<<"ms"<<endl;
        /*for(int i=0;i<Arr_size;i++){
            for(int j=0;j<Arr_size;j++)
                cout<<A[i][j]<<" ";
            cout<<endl;
        }*/
    }
    else
    {
        if(my_id!=num_proc-1){
            for(int j=0;j<block_size;j++)
                MPI_Recv(&A[my_id*block_size+j],Arr_size,MPI_FLOAT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        else{
            for(int j=0;j<block_size+remain;j++)
                MPI_Recv(&A[my_id*block_size+j],Arr_size,MPI_FLOAT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }

        eliminate(A,my_id,num_proc);

        if(my_id!=num_proc-1){
            for(int j=0;j<block_size;j++)
                MPI_Send(&A[my_id*block_size+j],Arr_size,MPI_FLOAT,0,1,MPI_COMM_WORLD);
        }
        else{
            for(int j=0;j<block_size+remain;j++)
                MPI_Send(&A[my_id*block_size+j],Arr_size,MPI_FLOAT,0,1,MPI_COMM_WORLD);
        }
    }

}

void MPI_Parallel_SIMD_OPENMP() {
    timeval t_start, t_end;

    int my_id,num_proc;

    MPI_Comm_rank(MPI_COMM_WORLD,&my_id);
    MPI_Comm_size(MPI_COMM_WORLD,&num_proc);

    int block_size = Arr_size / num_proc;
    int remain=Arr_size%num_proc;

    if(my_id==0){
        reset(A,Gauss_arr);
        gettimeofday(&t_start, NULL);
        for(int i=1;i<num_proc;i++)
        {
            if(i!=num_proc-1){
                for(int j=0;j<block_size;j++)
                    MPI_Send(&A[i*block_size+j],Arr_size,MPI_FLOAT,i,0,MPI_COMM_WORLD);
            }
            else{
                for(int j=0;j<block_size+remain;j++)
                    MPI_Send(&A[i*block_size+j],Arr_size,MPI_FLOAT,i,0,MPI_COMM_WORLD);
            }
        }

        eliminate_SIMD_OPENMP(A,my_id,num_proc);
        for(int i=1;i<num_proc;i++)
        {
            if(i!=num_proc-1){
                for(int j=0;j<block_size;j++)
                    MPI_Recv(&A[i*block_size+j],Arr_size,MPI_FLOAT,i,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
            else{
                for(int j=0;j<block_size+remain;j++)
                    MPI_Recv(&A[i*block_size+j],Arr_size,MPI_FLOAT,i,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
        }
        gettimeofday(&t_end, NULL);
        cout<<"MPI_BLOCK_SIMD_OPENMP time: "<<(t_end.tv_sec-t_start.tv_sec)*1000+(t_end.tv_usec-t_start.tv_usec)*0.001<<"ms"<<endl;
        /*for(int i=0;i<Arr_size;i++){
            for(int j=0;j<Arr_size;j++)
                cout<<A[i][j]<<" ";
            cout<<endl;
        }*/
    }
    else
    {
        if(my_id!=num_proc-1){
            for(int j=0;j<block_size;j++)
                MPI_Recv(&A[my_id*block_size+j],Arr_size,MPI_FLOAT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        else{
            for(int j=0;j<block_size+remain;j++)
                MPI_Recv(&A[my_id*block_size+j],Arr_size,MPI_FLOAT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }

        eliminate_SIMD_OPENMP(A,my_id,num_proc);

        if(my_id!=num_proc-1){
            for(int j=0;j<block_size;j++)
                MPI_Send(&A[my_id*block_size+j],Arr_size,MPI_FLOAT,0,1,MPI_COMM_WORLD);
        }
        else{
            for(int j=0;j<block_size+remain;j++)
                MPI_Send(&A[my_id*block_size+j],Arr_size,MPI_FLOAT,0,1,MPI_COMM_WORLD);
        }
    }
}

void MPI_Parallel_cycleSchedule(){
    timeval t_start, t_end;

    int my_id,num_proc;

    MPI_Comm_rank(MPI_COMM_WORLD,&my_id);
    MPI_Comm_size(MPI_COMM_WORLD,&num_proc);

    int steps= num_proc;

    if(my_id==0){
        reset(A,Gauss_arr);
        gettimeofday(&t_start, NULL);
        for(int i=0;i<Arr_size;i++)
        {
            int flag=i%steps;
            if(flag==my_id)
                continue;
            MPI_Send(&A[i],Arr_size,MPI_FLOAT,flag,0,MPI_COMM_WORLD);
        }
        eliminate_cycle(A,my_id,num_proc);
        for(int i=0;i<Arr_size;i++)
        {
            int flag=i%steps;
            if(flag==my_id)
                continue;
            MPI_Recv(&A[i],Arr_size,MPI_FLOAT,flag,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        gettimeofday(&t_end, NULL);
        cout<<"MPI_CYCLE time: "<<(t_end.tv_sec-t_start.tv_sec)*1000+(t_end.tv_usec-t_start.tv_usec)*0.001<<"ms"<<endl;
        /*for(int i=0;i<Arr_size;i++){
            for(int j=0;j<Arr_size;j++)
                cout<<A[i][j]<<" ";
            cout<<endl;
        }*/
    }
    else
    {
        for(int i=my_id;i<Arr_size;i+=steps)
            MPI_Recv(&A[i],Arr_size,MPI_FLOAT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        eliminate_cycle(A,my_id,num_proc);

        for(int i=my_id;i<Arr_size;i+=steps)
            MPI_Send(&A[i],Arr_size,MPI_FLOAT,0,1,MPI_COMM_WORLD);
    }
}

void MPI_Parallel_cycleSchedule_task(){
    timeval t_start, t_end;

    int my_id,num_proc;

    MPI_Comm_rank(MPI_COMM_WORLD,&my_id);
    MPI_Comm_size(MPI_COMM_WORLD,&num_proc);

    int steps= num_proc*task_per_step;

    if(my_id==0){
        reset(A,Gauss_arr);
        gettimeofday(&t_start, NULL);
        for(int i=0;i<Arr_size;i++)
        {
            int flag=i%steps;
            if(flag==my_id)
                continue;
            MPI_Send(&A[i],Arr_size,MPI_FLOAT,flag,0,MPI_COMM_WORLD);
        }
        eliminate_cycle_task(A,my_id,num_proc);
        for(int i=0;i<Arr_size;i++)
        {
            int flag=i%steps;
            if(flag==my_id)
                continue;
            MPI_Recv(&A[i],Arr_size,MPI_FLOAT,flag,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        gettimeofday(&t_end, NULL);
        cout<<"MPI_CYCLE time: "<<(t_end.tv_sec-t_start.tv_sec)*1000+(t_end.tv_usec-t_start.tv_usec)*0.001<<"ms"<<endl;
        /*for(int i=0;i<Arr_size;i++){
            for(int j=0;j<Arr_size;j++)
                cout<<A[i][j]<<" ";
            cout<<endl;
        }*/
    }
    else
    {
        for(int i=my_id;i<Arr_size;i+=steps)
            MPI_Recv(&A[i],Arr_size,MPI_FLOAT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        eliminate_cycle_task(A,my_id,num_proc);

        for(int i=my_id;i<Arr_size;i+=steps)
            MPI_Send(&A[i],Arr_size,MPI_FLOAT,0,1,MPI_COMM_WORLD);
    }
}

void MPI_Parallel_cycleSchedule_SIMD_OPENMP(){
    timeval t_start, t_end;

    int my_id,num_proc;

    MPI_Comm_rank(MPI_COMM_WORLD,&my_id);
    MPI_Comm_size(MPI_COMM_WORLD,&num_proc);

    int steps= num_proc;

    if(my_id==0){
        reset(A,Gauss_arr);
        gettimeofday(&t_start, NULL);
        for(int i=0;i<Arr_size;i++)
        {
            int flag=i%steps;
            if(flag==my_id)
                continue;
            MPI_Send(&A[i],Arr_size,MPI_FLOAT,flag,0,MPI_COMM_WORLD);
        }

        eliminate_cycle_SIMD_OPENMP(A,my_id,num_proc);
        for(int i=0;i<Arr_size;i++)
        {
            int flag=i%steps;
            if(flag==my_id)
                continue;
            MPI_Recv(&A[i],Arr_size,MPI_FLOAT,flag,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        gettimeofday(&t_end, NULL);
        cout<<"MPI_CYCLE_SIMD_OPENMP time: "<<(t_end.tv_sec-t_start.tv_sec)*1000+(t_end.tv_usec-t_start.tv_usec)*0.001<<"ms"<<endl;
       /* for(int i=0;i<Arr_size;i++){
            for(int j=0;j<Arr_size;j++)
                cout<<A[i][j]<<" ";
            cout<<endl;
        }*/
    }
    else
    {
        for(int i=my_id;i<Arr_size;i+=steps)
            MPI_Recv(&A[i],Arr_size,MPI_FLOAT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        eliminate_cycle_SIMD_OPENMP(A,my_id,num_proc);

        for(int i=my_id;i<Arr_size;i+=steps)
            MPI_Send(&A[i],Arr_size,MPI_FLOAT,0,1,MPI_COMM_WORLD);
    }

}

void MPI_Parallel_pipeline(){
    timeval t_start, t_end;

    int my_id,num_proc;

    MPI_Comm_rank(MPI_COMM_WORLD,&my_id);
    MPI_Comm_size(MPI_COMM_WORLD,&num_proc);

    int steps= num_proc;

    if(my_id==0){
        reset(A,Gauss_arr);
        gettimeofday(&t_start, NULL);
        for(int i=0;i<Arr_size;i++)
        {
            int flag=i%steps;
            if(flag==my_id)
                continue;
            MPI_Send(&A[i],Arr_size,MPI_FLOAT,flag,0,MPI_COMM_WORLD);
        }
        eliminate_cycle(A,my_id,num_proc);
        for(int i=0;i<Arr_size;i++)
        {
            int flag=i%steps;
            if(flag==my_id)
                continue;
            MPI_Recv(&A[i],Arr_size,MPI_FLOAT,flag,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        gettimeofday(&t_end, NULL);
        cout<<"MPI_PIPELINE time: "<<(t_end.tv_sec-t_start.tv_sec)*1000+(t_end.tv_usec-t_start.tv_usec)*0.001<<"ms"<<endl;
        for(int i=0;i<Arr_size;i++){
            for(int j=0;j<Arr_size;j++)
                cout<<A[i][j]<<" ";
            cout<<endl;
        }
    }
    else
    {
        for(int i=my_id;i<Arr_size;i+=steps)
            MPI_Recv(&A[i],Arr_size,MPI_FLOAT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        eliminate_cycle(A,my_id,num_proc);

        for(int i=my_id;i<Arr_size;i+=steps)
            MPI_Send(&A[i],Arr_size,MPI_FLOAT,0,1,MPI_COMM_WORLD);
    }
}

int main() {
    //MPI优化高斯消去

    static default_random_engine gengerator(1337);
    uniform_real_distribution<float> dis(-1.0, 1.0);

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
    reset(A,Gauss_arr);
    /*timeval t_start,t_end;
    gettimeofday(&t_start,NULL);
    Serial();
    gettimeofday(&t_end,NULL);
    cout<<"Serial time: "<<(t_end.tv_sec-t_start.tv_sec)*1000+(t_end.tv_usec-t_start.tv_usec)*0.001<<"ms"<<endl;
    for(int i=0;i<Arr_size;i++){
        for(int j=0;j<Arr_size;j++)
            cout<<A[i][j]<<" ";
        cout<<endl;
    }*/

    reset(A,Gauss_arr);

    MPI_Init(NULL,NULL);

    MPI_Parallel();

    MPI_Parallel_cycleSchedule();

    MPI_Parallel_SIMD_OPENMP();

    MPI_Parallel_cycleSchedule_SIMD_OPENMP();

    MPI_Parallel_pipeline();

    MPI_Finalize();

    return 0;
}


