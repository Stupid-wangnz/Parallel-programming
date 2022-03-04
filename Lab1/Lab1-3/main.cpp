#include <iostream>
#include<windows.h>
#include<stdlib.h>

using namespace std;

int main()
{

    int n=1024;

    int*a=new int[n];
    int sum=0;

    long long head,tail,freq;
    QueryPerformanceFrequency((LARGE_INTEGER *)&freq );
    QueryPerformanceCounter((LARGE_INTEGER *)&head);

    for(int i=0;i<n;i++)
        a[i]=i*i;

    for (int i = 0; i < n; i++)
        sum += a[i];
    QueryPerformanceCounter((LARGE_INTEGER *)&tail );
    cout<<"sum:"<<sum<<endl;
    cout<<"平凡算法时间："<<(tail-head)*1000.0/freq<<"ms"<<endl;

    int sum1 = 0, sum2 = 0;
    QueryPerformanceCounter((LARGE_INTEGER *)&head);
    for (int i = 0;i < n; i += 2) {
        sum1 += a[i];
        sum2 += a[i + 1];
    }
    QueryPerformanceCounter((LARGE_INTEGER *)&tail );
    sum = sum1 + sum2;
    cout<<"sum:"<<sum<<endl;
    cout<<"多链路式算法时间："<<(tail-head)*1000.0/freq<<"ms"<<endl;



    QueryPerformanceCounter((LARGE_INTEGER *)&head);
    for (int m = n; m > 1; m /= 2) // log(n)个步骤
        for (int i = 0; i < m / 2; i++)
            a[i] = a[i * 2] + a[i * 2 + 1];
    QueryPerformanceCounter((LARGE_INTEGER *)&tail );
    cout<<"sum:"<<a[0]<<endl;
    cout<<"二重循环算法时间："<<(tail-head)*1000.0/freq<<"ms"<<endl;
    return 0;
}
