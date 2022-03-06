#include <iostream>
#include<windows.h>
#include<stdlib.h>
#include<cstdlib>
#include<ctime>
using namespace std;

int main()
{
    int n=2000;

    cout<<"数组规模："<<n<<endl;
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

    long long head ,tail,freq,time=0;

    QueryPerformanceFrequency((LARGE_INTEGER *)&freq);




    int c=10;//重复执行的次数，减小误差


    for(int k=c;k>=0;k--)
    {

        QueryPerformanceCounter((LARGE_INTEGER *)&head);
        for(int i=0;i<n;i++)
        {
            sum[i]=0;
            for(int j=0;j<n;j++)
                sum[i]+=A[j][i]*a[j];

        }
        QueryPerformanceCounter((LARGE_INTEGER *)&tail );
        time+=tail-head;
    }


    cout<<"Col:"<<time*1000.0/freq/c<<"ms"<<endl;
    QueryPerformanceFrequency((LARGE_INTEGER *)&freq );
    time=0;
    for(int k=c;k>=0;k--)
    {

        QueryPerformanceCounter((LARGE_INTEGER *)&head);

        for(int i=0;i<n;i++)
            sum[i]=0;
        for(int j=0;j<n;j++)
            for(int i=0;i<n;i++)
                sum[i]=A[j][i]*a[j];

        QueryPerformanceCounter((LARGE_INTEGER *)&tail );
        time+=tail-head;
    }

    cout<<"Row:"<<time*1000.0/freq/c<<"ms";

    return 0;
}
