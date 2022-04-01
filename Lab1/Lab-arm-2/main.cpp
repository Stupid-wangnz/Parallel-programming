#include <iostream>
#include<sys/time.h>

#include <stdio.h>
#include <ctime>
#include<stdlib.h>
using namespace std;

void recursion(int n,int*a)

 {
    if (n == 1)
        return;
else
    {
        for (int i = 0; i < n / 2; i++)
            a[i] += a[n - i - 1];
        n = n / 2;
    recursion(n,a);
    }
}
void reset(int*a,int*b,int n)
{
    for(int i=0;i<n;i++)
        b[i]=a[i];

}

void Proc(int n)
{
    int*a=new int[n];
    int sum=0;

    int c=1;


     for(int i=0;i<n;i++)
    {
        srand(time(NULL));

        a[i]=rand()%100+1;
    }

    double time=0;
    struct timeval tv_begin,tv_end;
    for(int k=c;k>0;k--)
    {
        gettimeofday(&tv_begin,NULL);

        for (int i = 0; i < n; i++)
            sum += a[i];

        gettimeofday(&tv_end,NULL);
        time+=(tv_end.tv_sec-tv_begin.tv_sec)*1000.0+(tv_end.tv_usec-tv_begin.tv_usec)/1000.0;
    }
    cout<<"sum:"<<sum/c<<endl;
    cout<<"平凡算法时间："<<time/c<<"ms"<<endl;
    return;
    int sum1 = 0, sum2 = 0;
    time=0;
    for(int k=c;k>0;k--)
    {

        gettimeofday(&tv_begin,NULL);
        for (int i = 0;i < n; i += 2)
        {
            sum1 += a[i];
            sum2 += a[i + 1];
        }
        gettimeofday(&tv_end,NULL);
        time+=(tv_end.tv_sec-tv_begin.tv_sec)*1000.0+(tv_end.tv_usec-tv_begin.tv_usec)/1000.0;
    }

    sum = sum1 + sum2;
    cout<<"sum:"<<sum/c<<endl;
    cout<<"多链路式算法时间："<<time/c<<"ms"<<endl;

    int*b=new int[n];
     time=0;
     for(int k=c;k>0;k--)
    {
        reset(a,b,n);
        gettimeofday(&tv_begin,NULL);
        recursion(n,b);

        gettimeofday(&tv_end,NULL);
        time+=(tv_end.tv_sec-tv_begin.tv_sec)*1000.0+(tv_end.tv_usec-tv_begin.tv_usec)/1000.0;
    }

    cout<<"sum:"<<b[0]<<endl;
    cout<<"递归函数算法时间："<<time/c<<"ms"<<endl;

     time=0;
    for(int k=c;k>0;k--)
    {
        reset(a,b,n);
        gettimeofday(&tv_begin,NULL);
        for (int m = n; m > 1; m /= 2)
            for (int i = 0; i < m / 2; i++)
               b[i] = b[i * 2] + b[i * 2 + 1];

        gettimeofday(&tv_end,NULL);
        time+=(tv_end.tv_sec-tv_begin.tv_sec)*1000.0+(tv_end.tv_usec-tv_begin.tv_usec)/1000.0;
    }

    cout<<"sum:"<<b[0]<<endl;
    cout<<"二重循环算法时间："<<time/c<<"ms"<<endl;

}
int main()
{

    int n=1024*4;

    int*a=new int[n];
    int sum=0;

    int c=500;


     for(int i=0;i<n;i++)
    {
        srand(time(NULL));

        a[i]=rand()%100+1;
    }

    double time=0;
    struct timeval tv_begin,tv_end;

    int sum1 = 0, sum2 = 0;
    time=0;
    for(int k=c;k>0;k--)
    {
        gettimeofday(&tv_begin,NULL);

        for (int i = 0; i < n; i++)
            sum += a[i];

        gettimeofday(&tv_end,NULL);
        time+=(tv_end.tv_sec-tv_begin.tv_sec)*1000.0+(tv_end.tv_usec-tv_begin.tv_usec)/1000.0;
    }
    cout<<"sum:"<<sum/c<<endl;
    cout<<"平凡算法时间："<<time/c<<"ms"<<endl;


    return 0;
}
