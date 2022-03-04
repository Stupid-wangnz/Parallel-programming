#include <iostream>
#include<sys/time.h>

#include <stdio.h>
#include <unistd.h>
using namespace std;

int main()
{

    int n=1024;

    int*a=new int[n];
    int sum=0;

    struct timeval tv_begin,tv_end;


    for(int i=0;i<n;i++)
        a[i]=i*i;
    gettimeofday(&tv_begin, NULL);

    long head=tv_begin.tv_usec;
    for (int i = 0; i < n; i++)
        sum += a[i];

    gettimeofday(&tv_end, NULL);
    long tail=tv_end.tv_usec;

    cout<<"sum:"<<sum<<endl;
    cout<<"平凡算法时间："<<(tail-head)/1000.0<<"ms"<<endl;

    int sum1 = 0, sum2 = 0;
    gettimeofday(&tv_begin, NULL);

    head=tv_begin.tv_usec;
    for (int i = 0;i < n; i += 2) {
        sum1 += a[i];
        sum2 += a[i + 1];
    }
    gettimeofday(&tv_end, NULL);
    tail=tv_end.tv_usec;
    sum = sum1 + sum2;
    cout<<"sum:"<<sum<<endl;
    cout<<"多链路式算法时间："<<(tail-head)/1000.0<<"ms"<<endl;



    gettimeofday(&tv_begin, NULL);

    head=tv_begin.tv_usec;
    for (int m = n; m > 1; m /= 2) // log(n)个步骤
        for (int i = 0; i < m / 2; i++)
            a[i] = a[i * 2] + a[i * 2 + 1];
    gettimeofday(&tv_end, NULL);
    tail=tv_end.tv_usec;
    cout<<"sum:"<<a[0]<<endl;
    cout<<"二重循环算法时间："<<(tail-head)/1000.0<<"ms"<<endl;
    return 0;
}
