#include<iostream>
#include<fstream>
#include<sstream>
#include<pthread.h>
#include<ctime>
#include<cstdlib>
#include<sys/time.h>
#include<arm_neno.h>

using namespace std;


//
// Created by LEGION on 2022-04-26.
//
#include <vector>

class Grobner_Matrix {
    /*
     * Grobner_Matrix实现了对Grobner基问题求解中对矩阵的压缩
     * 原矩阵为n*n，因为只在（0,1）域上计算，所以可以压缩到n*n/32大小的矩阵.
     */
public:
    int n;//总行数
    int m;//总列数

    int m_;//压缩后的列数

    vector<int> row_index;//行索引,用于记录消元子或被消元子是否存在

    int** matrix;//原矩阵
    Grobner_Matrix();
    explicit Grobner_Matrix(int n,int m);
    ~Grobner_Matrix();
    void init();
    void set_bit(int i,int j);
    int xor_line(Grobner_Matrix&,int i,int j);
    int Simd_xor_line(Grobner_Matrix &grobnerMatrix, int i, int j)
    int get_max_bit(int i);
    void input_line(int i,vector<int>&);

    vector<int> get_5_line(int s);
    vector<int> get_line(int s,int num);
    void print_line(int i);
};

Grobner_Matrix::Grobner_Matrix() {
    n = 0;
    m = 0;
}

Grobner_Matrix::Grobner_Matrix(int n,int m) {
    this->n = n;
    this->m = m;
    init();
}

Grobner_Matrix::~Grobner_Matrix() {
    delete[] matrix;
}

void Grobner_Matrix::init() {
    matrix = new int*[n];
    if(m%32==0)
        m_ = m/32;
    else
        m_ = m/32+1;

    for (int i = 0; i <n ; i++) {
        matrix[i] = new int[m_];
    }
    row_index.resize(n,-1);
}

void Grobner_Matrix::set_bit(int i, int j) {
    //将第i行的第j个bit置1
    matrix[i][j/32] |= (1<<(j%32));
}

int Grobner_Matrix::xor_line(Grobner_Matrix &grobnerMatrix,int i,int j) {
    //返回异或后最大的非零位
    int max_bit=0;
    for (int k = 0; k < m_; k++) {
        matrix[j][k] ^= grobnerMatrix.matrix[i][k];
    }
    max_bit = get_max_bit(j);
    return max_bit;
}

void Grobner_Matrix::input_line(int i, vector<int> &line) {
    //将第i行的数据输入到矩阵中
    for(int j : line){
        set_bit(i,j);
    }
    row_index[i] = i;
}

vector<int> Grobner_Matrix::get_5_line(int s) {
    //从s开始，取得前五行
    vector<int> top5_line;
    top5_line.resize(row_index.size()-5*s>=5?5:row_index.size()-5*s);
    for(int i=0;i<top5_line.size();i++){
        top5_line[i] = row_index[i+5*s];
    }
    return top5_line;
}

vector<int> Grobner_Matrix::get_line(int s,int num) {
    //从s开始，取得前num行
    vector<int> top_line;
    top_line.resize(row_index.size()-num*s>=num?num:row_index.size()-num*s);
    for(int i=0;i<top_line.size();i++){
        top_line[i] = row_index[i+num*s];
    }
    return top_line;
}

int Grobner_Matrix::get_max_bit(int i) {
    //获取第i行最大的非零位
    int max_bit=0;
    for (int j = 0; j < m_; j++) {
        if(matrix[i][j]!=0) {
            for (int k = 0; k < 32; k++) {
                if (matrix[i][j] & (1 << k))
                    max_bit= j * 32 + k;
            }
        }
    }
    return max_bit;
}

void Grobner_Matrix::print_line(int i) {
    //将第i行各个为1的位，从大到小输出
    for (int j = 0; j < m_; j++) {
        if(matrix[i][j]!=0) {
            for (int k = 0; k < 32; k++) {
                if (matrix[i][j] & (1 << k))
                    cout << j * 32 + k << " ";
            }
        }
    }
    cout << endl;
}

int Grobner_Matrix::Simd_xor_line(Grobner_Matrix &grobnerMatrix, int i, int j) {
    //返回异或后最大的非零位
    int max_bit=0;
    int k=0;
    for (k = 0; k+4<=m_; k+=4) {
        int32x4_t mjk=vld1q_s32(&matrix[j][k]);
        int32x4_t mik=vld1q_s32(&grobnerMatrix.matrix[i][k]);
        mjk=veorq_s32(mjk,mik);
        vst1q_s32(&matrix[j][k],mjk);
        //matrix[j][k] ^= grobnerMatrix.matrix[i][k];
    }
    while(k<m_){
        matrix[j][k] ^= grobnerMatrix.matrix[i][k];
        k++;
    }
    max_bit = get_max_bit(j);
    return max_bit;
}

typedef struct {
    int t_id;
    int n;
    vector<int> v;
}threadParam_t;


int NUM_THREADS=4;

int n=2362;
int eliminator_count=170,elminated_element=453;

Grobner_Matrix eliminator(n,n);
Grobner_Matrix eliminatedElement(elminated_element,n);


void * threadfunc(void*param){
    threadParam_t * p=(threadParam_t*)param;
    int t_id=p->t_id;
    int step=NUM_THREADS;
    vector<int> row_index=p->v;
    for(int i=t_id;i<row_index.size();i+=step){
        int old_max_bit=eliminatedElement.get_max_bit(row_index[i]);
        int new_max_bit;
        while(eliminator.row_index[old_max_bit]!=-1){
            //消元子还没消元完
            new_max_bit=eliminatedElement.xor_line(eliminator,old_max_bit,row_index[i]);
            //如果new_max_bit不在eliminator的row_index中，则被消元子消元完毕
            if(eliminator.row_index[new_max_bit]==-1){
                break;
            }
            old_max_bit=new_max_bit;
        }
    }

}
void func_pthread(){
    fstream infile_eliminator(R"(F:\CL_WorkSpace\Grobner\input1.txt)");
    fstream infile_eliminated_element(R"(F:\CL_WorkSpace\Grobner\input2.txt)");

    if(!infile_eliminator.is_open()||!infile_eliminated_element.is_open()){
        cout<<"err open";
        return ;
    }

    string line;
    while(getline(infile_eliminator,line)){
        stringstream ss(line);
        vector<int>input_line;
        int temp;
        while(ss>>temp){
            input_line.push_back(temp);
        }
        eliminator.input_line(input_line[0],input_line);
    }

    line="";
    int record=0;
    while(getline(infile_eliminated_element,line)){
        stringstream ss(line);
        vector<int>input_line;

        int temp;
        while(ss>>temp){
            input_line.push_back(temp);
        }

        eliminatedElement.input_line(record,input_line);
        record++;
    }
    infile_eliminated_element.close();
    infile_eliminator.close();

    threadParam_t *param=new threadParam_t[NUM_THREADS];
    for(int i=0;i<NUM_THREADS;i++){
        param[i].t_id=i;
    }
    pthread_t *pthreads=new pthread_t[NUM_THREADS];

    double time=0;
    struct timeval tv_begin,tv_end;
    gettimeofday(&tv_begin,NULL);

    while(eliminatedElement.row_index.size()>0){
        //仍有被消元子没被消元完
        vector<int> row_index=eliminatedElement.row_index;
        for(int j=0;j<NUM_THREADS;j++){
            param[j].v=row_index;
        }
        for(int j=0;j<NUM_THREADS;j++){
            pthread_create(&pthreads[j],NULL,threadfunc,&param[j]);
        }
        for(int j=0;j<NUM_THREADS;j++){
            pthread_join(pthreads[j],NULL);
        }
        //线程计算完毕，更新eliminator
        for(int i=0;i<row_index.size();i++){
            int new_max_bit=eliminatedElement.get_max_bit(row_index[i]);
            if(new_max_bit==0){
                //消元子消元完毕
                eliminatedElement.row_index[i]=-1;
                continue;
            }
            if(eliminator.row_index[new_max_bit]==-1){
                eliminator.row_index[new_max_bit]=new_max_bit;
                eliminatedElement.row_index[i]=-1;
                for(int k=0;k<eliminatedElement.m_;k++){
                    eliminator.matrix[new_max_bit][k]=eliminatedElement.matrix[row_index[i]][k];
                }
            }
        }
        vector<int>new_row_index;
        for(int i=0;i<eliminatedElement.row_index.size();i++){
            if(eliminatedElement.row_index[i]!=-1){
                new_row_index.push_back(eliminatedElement.row_index[i]);
            }
        }
        eliminatedElement.row_index=new_row_index;
    }
    gettimeofday(&tv_end,NULL);
    time+=(tv_end.tv_sec-tv_begin.tv_sec)*1000.0+(tv_end.tv_usec-tv_begin.tv_usec)/1000.0;
    cout<<"time:"<<time<<endl;

    delete []param;
    delete []pthreads;
}

int main(){

    /*
     * n是总行数
     *
     * eliminator是初始时消元子的数量
     *
     * eliminated_element是初始时被消元子的数量
     */

    fstream infile_eliminator(R"(F:\CL_WorkSpace\Grobner\input1.txt)");
    fstream infile_eliminated_element(R"(F:\CL_WorkSpace\Grobner\input2.txt)");

    if(!infile_eliminator.is_open()||!infile_eliminated_element.is_open()){
        cout<<"err open";
        return 0;
    }

    string line;
    while(getline(infile_eliminator,line)){
        stringstream ss(line);
        vector<int>input_line;
        int temp;
        while(ss>>temp){
            input_line.push_back(temp);
        }
        eliminator.input_line(input_line[0],input_line);
    }

    line="";
    int record=0;
    while(getline(infile_eliminated_element,line)){
        stringstream ss(line);
        vector<int>input_line;

        int temp;
        while(ss>>temp){
            input_line.push_back(temp);
        }

        eliminatedElement.input_line(record,input_line);
        record++;
    }
    infile_eliminated_element.close();
    infile_eliminator.close();

    //串行的版本
    double time=0;
    struct timeval tv_begin,tv_end;
    gettimeofday(&tv_begin,NULL);
    int num=0;
    while(num*4<elminated_element){
        //仍有被消元子没被消元完
        vector<int> row_index=eliminatedElement.get_line(num,4);
        for(int i=0;i<row_index.size();i++){
            int old_max_bit=eliminatedElement.get_max_bit(row_index[i]);
            int new_max_bit;
            while(eliminator.row_index[old_max_bit]!=-1){
                //消元子还没消元完
                new_max_bit=eliminatedElement.xor_line(eliminator,old_max_bit,row_index[i]);
                if(new_max_bit==0){
                    //消元子消元完了
                    //cout<<row_index[i]<<": all0"<<endl;
                    break;
                }
                //如果new_max_bit不在eliminator的row_index中，则被消元子消元完毕
                if(eliminator.row_index[new_max_bit]==-1){
                    for(int k=0;k<eliminatedElement.m_;k++){
                        eliminator.matrix[new_max_bit][k]=eliminatedElement.matrix[row_index[i]][k];
                    }
                    eliminator.row_index[new_max_bit]=new_max_bit;
                    break;
                }
                old_max_bit=new_max_bit;
            }
        }
        num++;
    }

    gettimeofday(&tv_end,NULL);
    time+=(tv_end.tv_sec-tv_begin.tv_sec)*1000.0+(tv_end.tv_usec-tv_begin.tv_usec)/1000.0;
    cout<<"time:"<<time<<endl;

    func_pthread();

    return 0;
}

