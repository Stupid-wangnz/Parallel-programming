#include <iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<stdio.h>
#include<ctime>
#include<stdlib.h>
#include<cstdlib>
#include<sys/time.h>
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
    for(int i=0;i+4<=n;i+=4)
        {
            int16x4_t a=vld1_s16(&A[i]);
            int16x4_t b=vld1_s16(&B[i]);
            b=veor_s16(a,b);
            vst1_s16(&B[i],b);
        }

}
int main()
{

    ifstream infile1("/home/data/Groebner/1_130_22_8/1.txt",ios::in);
    if(!infile1){
        cout<<"not open!";
        return 0;
    }
    ifstream infile2("/home/data/Groebner/1_130_22_8/2.txt",ios::in);
    if(!infile2){
        cout<<"not open!";
        return 0;
    }


    int n=130,n1=22,n2=8;

    //n列数,n1消元子的数量，n2消元子的数量

    int *p=new int[n2];
    for(int i=0;i<n2;i++)
        p[i]=1;

    int *record=new int[n];
    for(int i=0;i<n;i++)
        record[i]=0;

    short ** arr1=new short*[n];
    for(int i=0;i<n;i++)
        arr1[i]=new short[n];

    short ** arr2=new short*[n2];
    for(int i=0;i<n2;i++)
        arr2[i]=new short[n];

    for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
            arr1[i][j]=0;
    for(int i=0;i<n2;i++)
        for(int j=0;j<n;j++)
            arr2[i][j]=0;



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
    gettimeofday(&tv_begin,NULL);
    for(int i=0;i<n2;i++)
    {
        if(p[i]==0)
            continue;

        int k=FindTheMax(arr2[i],n);
        while(record[k])
        {
            XY(arr1[k],arr2[i],n);
            k=FindTheMax(arr2[i],n);
        }

        record[k]=1;
        p[i]=0;

        for(int x=0;x<n;x++)
            arr1[k][x]=arr2[i][x];
    }

    gettimeofday(&tv_end,NULL);
    time+=(tv_end.tv_sec-tv_begin.tv_sec)*1000.0+(tv_end.tv_usec-tv_begin.tv_usec)/1000.0;
    }
    cout<<time/20<<endl;

    time=0;
    for(int i=0;i<20;i++){
    gettimeofday(&tv_begin,NULL);
    for(int i=0;i<n2;i++)
    {
        if(p[i]==0)
            continue;

        int k=FindTheMax(arr2[i],n);
        while(record[k])
        {
            XY1(arr1[k],arr2[i],n);
            k=FindTheMax(arr2[i],n);
        }

        record[k]=1;
        p[i]=0;

        for(int x=0;x<n;x++)
            arr1[k][x]=arr2[i][x];
    }

    gettimeofday(&tv_end,NULL);
    time+=(tv_end.tv_sec-tv_begin.tv_sec)*1000.0+(tv_end.tv_usec-tv_begin.tv_usec)/1000.0;
    }
    cout<<time/20<<endl;
    return 0;
}
