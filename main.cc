#include <cstring>
#include <string>
#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <sys/time.h>
#include <omp.h>

#include<bits/stdc++.h>//直接加个万能头省事

#include<arm_neon.h>

// 可以自行添加需要的头文件

void fRead(int *a, int *b, int *n, int *p, int input_id){
    // 数据输入函数
    std::string str1 = "/nttdata/";
    std::string str2 = std::to_string(input_id);
    std::string strin = str1 + str2 + ".in";
    char data_path[strin.size() + 1];
    std::copy(strin.begin(), strin.end(), data_path);
    data_path[strin.size()] = '\0';
    std::ifstream fin;
    fin.open(data_path, std::ios::in);
    fin>>*n>>*p;
    for (int i = 0; i < *n; i++){
        fin>>a[i];
    }
    for (int i = 0; i < *n; i++){   
        fin>>b[i];
    }
}

void fCheck(int *ab, int n, int input_id){
    // 判断多项式乘法结果是否正确
    std::string str1 = "/nttdata/";
    std::string str2 = std::to_string(input_id);
    std::string strout = str1 + str2 + ".out";
    char data_path[strout.size() + 1];
    std::copy(strout.begin(), strout.end(), data_path);
    data_path[strout.size()] = '\0';
    std::ifstream fin;
    fin.open(data_path, std::ios::in);
    for (int i = 0; i < n * 2 - 1; i++){
        int x;
        fin>>x;
        if(x != ab[i]){
            std::cout<<"多项式乘法结果错误"<<std::endl;
            return;
        }
    }
    std::cout<<"多项式乘法结果正确"<<std::endl;
    return;
}

void fWrite(int *ab, int n, int input_id){
    // 数据输出函数, 可以用来输出最终结果, 也可用于调试时输出中间数组
    std::string str1 = "files/";
    std::string str2 = std::to_string(input_id);
    std::string strout = str1 + str2 + ".out";
    char output_path[strout.size() + 1];
    std::copy(strout.begin(), strout.end(), output_path);
    output_path[strout.size()] = '\0';
    std::ofstream fout;
    fout.open(output_path, std::ios::out);
    for (int i = 0; i < n * 2 - 1; i++){
        fout<<ab[i]<<'\n';
    }
}

void poly_multiply(int *a, int *b, int *ab, int n, int p){
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            ab[i+j]=(1LL * a[i] * b[j] % p + ab[i+j]) % p;
        }
    }
}




int ksm(int a,int b,int p){//快速幂modp
    if(!b)  return 1;
    int tmp=1;
    while(b>1){
        if(b&1) tmp=1ll*tmp*a%p;
        a=1ll*a*a%p;
        b>>=1;
    }
    return 1ll*a*tmp%p;
}


int k=0,lim=1;//k表示最大幂次，lim表示（大于等于2n的）2^k
int qwq[400000];//辅助数组，用于存储蝴蝶变换下标
void realntt(int *a,int p,bool opt){//算法主体,opt表示正运算还是逆运算
    
    for(int i=0;i<lim;++i){
        if(i<qwq[i])    std::swap(a[i],a[qwq[i]]);
    }

    for(int mid=1;mid<lim;mid<<=1){
        int wn=ksm(3,(p-1)/(mid<<1),p);
        if(!opt){//负运算，需要求逆元
            wn=ksm(wn,p-2,p);
        }

        for(int j=0;j<lim;j+=(mid<<1)){//主体
            int w=1;
            for(int k=0;k<mid;++k,w=1ll*w*wn%p){
                int x=a[j+k];
                int y=1ll*a[j+k+mid]*w%p;
                a[j+k]=(x+y)%p;
                a[j+k+mid]=(x-y+p)%p;
            }
        }
    }

    if(!opt){//负运算除以长度
        int awa=ksm(lim,p-2,p);
        for(int i=0;i<lim;++i){
            a[i]=1ll*a[i]*awa%p;
        }
    }
    return;
}
void ntt(int *a,int *b,int *ab,int n,int p){//优化算法

    while(lim<n*2){//初始化幂次
        lim<<=1;
        ++k;
    }

    for(int i=0;i<lim;i++){//对每个位置下标的二进制位操作 
		qwq[i]=((qwq[i>>1]>>1)|((i&1)<<(k-1)));
	}

    realntt(a,p,true);
    realntt(b,p,true);

    for(int i=0;i<lim;++i){
        ab[i]=1ll*a[i]*b[i]%p;
    }
    realntt(ab,p,false);//INTT
    return;
}






int a[500000], b[500000], ab[500000];
int main(int argc, char *argv[])
{
    
    // 保证输入的所有模数的原根均为 3, 且模数都能表示为 a \times 4 ^ k + 1 的形式
    // 输入模数分别为 7340033 104857601 469762049 263882790666241
    // 第四个模数超过了整型表示范围, 如果实现此模数意义下的多项式乘法需要修改框架
    // 对第四个模数的输入数据不做必要要求, 如果要自行探索大模数 NTT, 请在完成前三个模数的基础代码及优化后实现大模数 NTT
    // 输入文件共五个, 第一个输入文件 n = 4, 其余四个文件分别对应四个模数, n = 131072
    // 在实现快速数论变化前, 后四个测试样例运行时间较久, 推荐调试正确性时只使用输入文件 1
    int test_begin = 0;
    int test_end = 1;
    for(int i = test_begin; i <= test_end; ++i){
        long double ans = 0;
        int n_, p_;

        memset(a,0,sizeof(a));
        memset(b,0,sizeof(b));//防止数据污染而加的清零

        fRead(a, b, &n_, &p_, i);
        memset(ab,0,sizeof(ab));
        auto Start = std::chrono::high_resolution_clock::now();
        // TODO : 将 poly_multiply 函数替换成你写的 ntt
    //    poly_multiply(a, b, ab, n_, p_);


        ntt(a,b,ab,n_,p_);//执行优化算法


        auto End = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double,std::ratio<1,1000>>elapsed = End - Start;
        ans += elapsed.count();
        fCheck(ab, n_, i);
        std::cout<<"average latency for n = "<<n_<<" p = "<<p_<<" : "<<ans<<" (us) "<<std::endl;
        // 可以使用 fWrite 函数将 ab 的输出结果打印到 files 文件夹下
        // 禁止使用 cout 一次性输出大量文件内容
        fWrite(ab, n_, i);
    }
    return 0;
}
