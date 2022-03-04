#include <iostream>
#include<windows.h>
#include<stdlib.h>
using namespace std;

int main()
{
    int n=1000;
    int**A=new int*[n];
    for(int i=0;i<n;i++)
        A[i]=new int[n];

    int *a=new int[n];

    for(int i=0;i<n;i++)
    {
        a[i]=i*i;

        for(int j=0;j<n;j++)
            {
                A[i][j]=i+j;
            }
    }
    int *sum=new int[n];

    long long head ,tail,freq;

    QueryPerformanceFrequency((LARGE_INTEGER *)&freq);

    QueryPerformanceCounter((LARGE_INTEGER *)&head);
    for(int i=0;i<n;i++)
    {
        sum[i]=0;
        for(int j=0;j<n;j++)
            sum[i]+=A[j][i]*a[j];

    }
    QueryPerformanceCounter((LARGE_INTEGER *)&tail );


    cout<<"Col:"<<(tail-head)*1000.0/freq<<"ms"<<endl;


    QueryPerformanceFrequency((LARGE_INTEGER *)&freq );
    QueryPerformanceCounter((LARGE_INTEGER *)&head);
    for(int i=0;i<n;i++)
        sum[i]=0;
    for(int j=0;j<n;j++)
        for(int i=0;i<n;i++)
            sum[i]=A[j][i]*a[j];

    QueryPerformanceCounter((LARGE_INTEGER *)&tail );
    cout<<"Row:"<<(tail-head)*1000.0/freq<<"ms";

    return 0;
}
