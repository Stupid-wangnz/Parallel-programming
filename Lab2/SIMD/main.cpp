#include <iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<stdio.h>
#include<ctime>
#include<stdlib.h>
#include<cstdlib>
#include<sys/time.h>
#include<arm_neon.h>
using namespace std;

bool Finish(int*A,int n2)
{
    for(int i=0;i<n2;i++)
        if(A[i])
            return 0;

    return 1;
}

int FindTheMax(short*A,int n)
{

    for(int i=0;i<n;i++)
        if(A[i])
            return i;
    return -1;
}

void XY(short*A,short*B,int n)
{
    for(int i=0;i<n;i++)
        B[i]=A[i]^B[i];

}
void XY1(short*A,short*B,int n)
{
    int i;
    for(i=0;i+4<=n;i+=4)
        {
            int16x4_t a=vld1_s16(&A[i]);
            int16x4_t b=vld1_s16(&B[i]);
            b=veor_s16(a,b);
            vst1_s16(&B[i],b);
        }
    while(i%4!=0){
        B[i]=A[i]^B[i];
        i++;
    }

}
void XY2(short*A,short*B,int n)
{
    int i;
    for(i=0;i+8<=n;i+=8)
        {
            int16x8_t a=vld1q_s16(&A[i]);
            int16x8_t b=vld1q_s16(&B[i]);
            b=veorq_s16(a,b);
            vst1q_s16(&B[i],b);
        }
    while(i%8!=0){
        B[i]=A[i]^B[i];
        i++;
    }
}

void reset(short**a,short**A,short**b,short**B,int n,int n2,int*record,int*r,int*p,int*pc)
{
    for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
            A[i][j]=a[i][j];

    for(int i=0;i<n2;i++)
        for(int j=0;j<n;j++)
            B[i][j]=b[i][j];

    for(int i=0;i<n;i++)
        r[i]=record[i];
    for(int i=0;i<n2;i++)
        pc[i]=p[i];
}
int main()
{


    ifstream infile1("/home/data/Groebner/2_254_106_53/1.txt",ios::in);
    if(!infile1){
        cout<<"not open!";
        return 0;
    }
    ifstream infile2("/home/data/Groebner/2_254_106_53/2.txt",ios::in);
    if(!infile2){
        cout<<"not open!";
        return 0;
    }


    int n=254,n1=106,n2=53;

    //n列数,n1消元子的数量，n2消元子的数量

    int *p=new int[n2];
    int *pc=new int[n2];
    for(int i=0;i<n2;i++)
        p[i]=1;
    int *record=new int[n];
    int *r=new int[n];
    for(int i=0;i<n;i++)
        record[i]=0;

    short ** arr1=new short*[n];
    for(int i=0;i<n;i++)
        arr1[i]=new short[n];

    short ** arr2=new short*[n2];
    for(int i=0;i<n2;i++)
        arr2[i]=new short[n];


    short ** A=new short*[n];
    for(int i=0;i<n;i++)
        A[i]=new short[n];

    short ** B=new short*[n2];
    for(int i=0;i<n2;i++)
        B[i]=new short[n];

    string temp;
    int l=0;
    while(getline(infile1,temp))
    {
        stringstream ss(temp);
        int num;
        ss>>num;
        num=n-num-1;
        record[num]=1;
        arr1[num][num]=1;
        int t;
        while(!ss.eof())
        {
            ss>>t;
            if(t<=0)
                break;
            arr1[num][n-t-1]=1;
        }
    }
    infile1.close();

    temp="";
    while(getline(infile2,temp))
    {
        stringstream ss(temp);
        int num;
        while(!ss.eof())
        {
            ss>>num;
            if(num<=0)
                break;
            arr2[l][n-num-1]=1;
        }
        l++;
    }
    infile2.close();


    double time=0;
    struct timeval tv_begin,tv_end;

    for(int i=0;i<20;i++){

    reset(arr1,A,arr2,B,n,n2,record,r,p,pc);
    gettimeofday(&tv_begin,NULL);
    for(int i=0;i<n2;i++)
    {
        if(pc[i]==0)
            continue;

        int k=FindTheMax(B[i],n);
        while(r[k])
        {
            XY(A[k],B[i],n);
            k=FindTheMax(B[i],n);
            if(k<0)
	break;
        }
        pc[i]=0;
        if(k<0)
            continue;
        r[k]=1;

        for(int x=0;x<n;x++)
            A[k][x]=B[i][x];

    }

    gettimeofday(&tv_end,NULL);
    time+=(tv_end.tv_sec-tv_begin.tv_sec)*1000.0+(tv_end.tv_usec-tv_begin.tv_usec)/1000.0;
    }
    cout<<time/20<<endl;
    time=0;
    for(int i=0;i<20;i++){

    reset(arr1,A,arr2,B,n,n2,record,r,p,pc);
    gettimeofday(&tv_begin,NULL);
    for(int i=0;i<n2;i++)
    {
        if(pc[i]==0)
            continue;

        int k=FindTheMax(B[i],n);
        while(r[k])
        {
            XY1(A[k],B[i],n);
            k=FindTheMax(B[i],n);
            if(k<0)
	break;
        }
        pc[i]=0;
        if(k<0)
            continue;
        r[k]=1;

        for(int x=0;x<n;x++)
            A[k][x]=B[i][x];
    }

    gettimeofday(&tv_end,NULL);
    time+=(tv_end.tv_sec-tv_begin.tv_sec)*1000.0+(tv_end.tv_usec-tv_begin.tv_usec)/1000.0;
    }
    cout<<time/20<<endl;


    for(int i=0;i<20;i++){

    reset(arr1,A,arr2,B,n,n2,record,r,p,pc);
    gettimeofday(&tv_begin,NULL);
    for(int i=0;i<n2;i++)
    {
        if(pc[i]==0)
            continue;

        int k=FindTheMax(B[i],n);
        while(r[k])
        {
            XY2(A[k],B[i],n);
            k=FindTheMax(B[i],n);
            if(k<0)
	break;
        }
        pc[i]=0;
        if(k<0)
            continue;
        r[k]=1;

        for(int x=0;x<n;x++)
            A[k][x]=B[i][x];

    }

    gettimeofday(&tv_end,NULL);
    time+=(tv_end.tv_sec-tv_begin.tv_sec)*1000.0+(tv_end.tv_usec-tv_begin.tv_usec)/1000.0;
    }
    cout<<time/20<<endl;
    time=0;
    return 0;
}
