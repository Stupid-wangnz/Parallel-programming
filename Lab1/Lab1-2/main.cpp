#include <iostream>

#include<sys/time.h>
#include<cstdlib>
#include<stdlib.h>
#include<ctime>
using namespace std;

int main()
{
    int n=5000;


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
    long head,tail;
    struct timeval tv_begin,tv_end;
    //gettimeofday(&tv_begin, NULL);
    //head=tv_begin.tv_usec;
    clock_t start,ed;

    start=clock();
    int c=100;//循环规模
    for(int k=c;k>=0;k--)
    {


        for(int i=0;i<n;i++)
        {
            sum[i]=0;
            for(int j=0;j<n;j++)
                sum[i]+=A[j][i]*a[j];

        }
    }
    //gettimeofday(&tv_end,NULL);

    //tail=tv_end.tv_usec;
    ed=clock();
    cout<<"Col:"<<(ed-start)*1000.0/c/CLOCKS_PER_SEC<<"ms"<<endl;

    //gettimeofday(&tv_begin, NULL);
    //head=tv_begin.tv_usec;
    start=clock();
    for(int k=c;k>=0;k--)
    {

        for(int i=0;i<n;i++)
            sum[i]=0;
        for(int j=0;j<n;j++)
            for(int i=0;i<n;i++)
                sum[i]=A[j][i]*a[j];

    }
    ed=clock();
    //gettimeofday(&tv_end,NULL);
    //tail=tv_end.tv_usec;
    cout<<"Row:"<<(ed-start)*1000.0/c/CLOCKS_PER_SEC<<"ms";

    return 0;
}
