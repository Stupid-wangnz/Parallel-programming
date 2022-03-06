#include <iostream>

#include<sys/time.h>
#include<cstdlib>
#include<stdlib.h>
#include<ctime>
using namespace std;

int main()
{
    int n=500;


    cout<<"数组规模:"<<n<<endl;
    int**A=new int*[n];
    for(int i=0;i<n;i++)
        A[i]=new int[n];

    int *a=new int[n];

    for(int i=0;i<n;i++)
    {
        srand(time(NULL));
        a[i]=rand()%10+1;

        for(int j=0;j<n;j++)
            {
                srand(time(NULL));
                A[i][j]=rand()%100+1;
            }
    }
    int *sum=new int[n];
    double time=0;
    struct timeval tv_begin,tv_end;

    int c=100;
    for(int k=c;k>=0;k--)
    {

        gettimeofday(&tv_begin,NULL);
        for(int i=0;i<n;i++)
        {
            sum[i]=0;
            for(int j=0;j<n;j++)
                sum[i]+=A[j][i]*a[j];

        }
        gettimeofday(&tv_end,NULL);
        time+=(tv_end.tv_sec-tv_begin.tv_sec)*1000.0+(tv_end.tv_usec-tv_begin.tv_usec)/1000.0;
    }


    cout<<"Col:"<<time/c<<"ms"<<endl;


    time=0;
    for(int k=c;k>=0;k--)
    {
        gettimeofday(&tv_begin,NULL);
        for(int i=0;i<n;i++)
            sum[i]=0;
        for(int j=0;j<n;j++)
            for(int i=0;i<n;i++)
                sum[i]=A[j][i]*a[j];

        gettimeofday(&tv_end,NULL);
        time+=(tv_end.tv_sec-tv_begin.tv_sec)*1000.0+(tv_end.tv_usec-tv_begin.tv_usec)/1000.0;

    }


    cout<<"Row:"<<time/c<<"ms";

    return 0;
}
