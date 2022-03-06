#include <iostream>
#include<sys/time.h>

#include <stdio.h>
#include <ctime>
#include<stdlib.h>
using namespace std;

int main()
{

    int n=1024;

    int*a=new int[n];
    int sum=0;

    int c=100;


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

    int sum1 = 0, sum2 = 0;

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


    for(int k=c;k>0;k--)
    {
        gettimeofday(&tv_begin,NULL);
        for (int m = n; m > 1; m /= 2)
            for (int i = 0; i < m / 2; i++)
                a[i] = a[i * 2] + a[i * 2 + 1];

        gettimeofday(&tv_end,NULL);
        time+=(tv_end.tv_sec-tv_begin.tv_sec)*1000.0+(tv_end.tv_usec-tv_begin.tv_usec)/1000.0;
    }

    cout<<"sum:"<<a[0]<<endl;
    cout<<"二重循环算法时间："<<time/c<<"ms"<<endl;
    return 0;
}
