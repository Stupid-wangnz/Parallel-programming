#include <iostream>

#include<sys/time.h>

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
    long head,tail;
    struct timeval tv_begin,tv_end;
    gettimeofday(&tv_begin, NULL);
    head=tv_begin.tv_usec;
    for(int i=0;i<n;i++)
    {
        sum[i]=0;
        for(int j=0;j<n;j++)
            sum[i]+=A[j][i]*a[j];

    }
    gettimeofday(&tv_end,NULL);
    tail=tv_end.tv_usec;

    cout<<"Col:"<<(tail-head)/1000.0<<"ms"<<endl;


    gettimeofday(&tv_begin, NULL);
    head=tv_begin.tv_usec;
    for(int i=0;i<n;i++)
        sum[i]=0;
    for(int j=0;j<n;j++)
        for(int i=0;i<n;i++)
            sum[i]=A[j][i]*a[j];
    gettimeofday(&tv_end,NULL);
    tail=tv_end.tv_usec;

    cout<<"Row:"<<(tail-head)/1000.0<<"ms";

    return 0;
}
