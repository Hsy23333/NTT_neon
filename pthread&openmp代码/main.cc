#include <cstring>
#include <string>
#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <sys/time.h>
#include <omp.h>

#include<bits/stdc++.h>//直接加个万能头省事

#include<cstdint>
#include<pthread.h>

// 可以自行添加需要的头文件

void fRead(long long *a,long long *b,long long *n,long long *p,long long input_id){
    // 数据输入函数
    std::string str1="/nttdata/";
    std::string str2=std::to_string(input_id);
    std::string strin=str1+str2+".in";
    char data_path[strin.size()+1];
    std::copy(strin.begin(),strin.end(),data_path);
    data_path[strin.size()]='\0';
    std::ifstream fin;
    fin.open(data_path,std::ios::in);
    fin>>*n>>*p;
    for (long long i=0; i < *n; i++){
        fin>>a[i];
    }
    for (long long i=0; i < *n; i++){   
        fin>>b[i];
    }
}

void fCheck(long long *ab,long long n,long long input_id){
    // 判断多项式乘法结果是否正确
    std::string str1="/nttdata/";
    std::string str2=std::to_string(input_id);
    std::string strout=str1+str2+".out";
    char data_path[strout.size()+1];
    std::copy(strout.begin(),strout.end(),data_path);
    data_path[strout.size()]='\0';
    std::ifstream fin;
    fin.open(data_path,std::ios::in);
    for (long long i=0; i < n*2-1; i++){
        long long x;
        fin>>x;
        if(x != ab[i]){
            std::cout<<"多项式乘法结果错误"<<std::endl;
            return;
        }
    }
    std::cout<<"多项式乘法结果正确"<<std::endl;
    return;
}

void fWrite(long long *ab,long long n,long long input_id){
    // 数据输出函数,可以用来输出最终结果,也可用于调试时输出中间数组
    std::string str1="files/";
    std::string str2=std::to_string(input_id);
    std::string strout=str1+str2+".out";
    char output_path[strout.size()+1];
    std::copy(strout.begin(),strout.end(),output_path);
    output_path[strout.size()]='\0';
    std::ofstream fout;
    fout.open(output_path,std::ios::out);
    for (long long i=0; i < n*2-1; i++){
        fout<<ab[i]<<'\n';
    }
}

void poly_multiply(long long *a,long long *b,long long *ab,long long n,long long p){
    for(long long i=0; i < n; ++i){
        for(long long j=0; j < n; ++j){
            ab[i+j]=(1ll*a[i]*b[j]%p+ab[i+j])%p;
        }
    }
}




long long ksm(long long a,long long b,long long p){//快速幂modp
    if(!b)  return 1;
    __int128 tmp=1,aa=a;

    while(b>1){
        if(b&1) tmp=tmp*aa%p;
        aa=aa*aa%p;
        b>>=1;
    }
    return (long long)(tmp*aa%p);
}



long long mod[4]={1004535809,1998585857,998244353,469762049};//多模数

long long k=0,lim=1;//k表示最大幂次，lim表示（大于等于2n的）2^k
long long qwq[400000];//辅助数组，用于存储蝴蝶变换下标

/*pthread_t ths[4];//朴素优化，使用4个线程
struct qwqqq{//封装参数，传给线程
    long long mid,wn,j,p;
    int KK;//KK=0,1,2,3
    long long *a;
}tmp4[4];//创建4个暂时参数结构体
void* pth(void* tmpw){
    qwqqq* tmpq=(qwqqq*) tmpw;
    long long mid=tmpq->mid;
    long long wn=tmpq->wn;
    long long j=tmpq->j;
    long long p=tmpq->p;
    long long w=ksm(wn,tmpq->KK,p);//不同的KK值对应起点不同
    long long wn4=ksm(wn,4,p);//改为转4下
    long long *a=tmpq->a;
    for(long long k=tmpq->KK;k<mid;k+=4,w=1ll*w*wn4%p){
        long long x=a[j+k];
        long long y=1ll*a[j+k+mid]*w%p;
        a[j+k]=(x+y)%p;
        a[j+k+mid]=(x-y+p)%p;
    }
    return NULL;
}*/
void realntt(long long *a,long long p,bool opt){//算法主体,opt表示正运算还是逆运算
    
    for(long long i=0;i<lim;++i){
        if(i<qwq[i])    std::swap(a[i],a[qwq[i]]);
    }

/*    for(int v=0;v<4;++v){//有些参数可以在外面初始化加速(朴素多线程)
        tmp4[v].p=p;
        tmp4[v].a=a;
    }*/
    for(long long mid=1;mid<lim;mid<<=1){
        long long wn=ksm(3,(p-1)/(mid<<1),p);
        if(!opt){//负运算，需要求逆元
            wn=ksm(wn,p-2,p);
        }


    /*    for(int v=0;v<4;++v){//有些参数可以在外面初始化加速（朴素多线程）
            tmp4[v].wn=wn;
            tmp4[v].mid=mid;
        }*/
        for(long long j=0;j<lim;j+=(mid<<1)){//主体
            long long w=1;
            for(long long k=0;k<mid;++k,w=1ll*w*wn%p){
                long long x=a[j+k];
                long long y=a[j+k+mid]*w%p;
                a[j+k]=(x+y)%p;
                a[j+k+mid]=(x-y+p)%p;
            }
            
/*
            //下面是朴素优化
            for(int v=0;v<4;++v){
                tmp4[v].j=j;
                tmp4[v].KK=v;
//                tmp4[v].mid=mid;
//                tmp4[v].p=p;
//                tmp4[v].wn=wn;
                int errcode=pthread_create(&ths[v],NULL,pth,&tmp4[v]);
            }
            for(int v=0;v<4;++v)    pthread_join(ths[v],NULL);//等待这层所有线程结束
           */ 
        }
    }
    if(!opt){//负运算除以长度
        long long awa=ksm(lim,p-2,p);
        for(long long i=0;i<lim;++i){
            a[i]=a[i]*awa%p;
        }
    }
    return;
}
void ntt(long long *a,long long *b,long long *ab,long long n,long long p){//优化算法
    realntt(a,p,true);
    realntt(b,p,true);
    for(long long i=0;i<lim;++i){
        ab[i]=a[i]*b[i]%p;
    }
    realntt(ab,p,false);//INTT
    return;
}




long long ac[4][500000];
long long bc[4][500000];//CRT的4个拆分
long long a[500000],b[500000],ab[500000];

pthread_t ths[4];//考虑对4模数NTT并行优化
struct qwqqq{//封装参数，传给线程
    long long n_;
    int j;//0,1,2,3
}tmp4[4];//创建4个暂时参数结构体
struct awaaa{
    __int128 M1M2,M1M2M3,invM1_modM2,invM1M2_modM3,invM1M2M3_modM4;
    int K;
    long long p_;
}tmpp4[4];
void* NTT_parallel(void* tmpw){
    qwqqq* tmpq=(qwqqq*) tmpw;
    int j=tmpq->j;
    long long n_=tmpq->n_;
    for(int J=0;J<=lim;++J){
        ac[j][J]=a[J]%mod[j];
        bc[j][J]=b[J]%mod[j];
    }
    ntt(ac[j],bc[j],ac[j],n_,mod[j]);
    return NULL;
}
void* CRT_parallel(void* tmpw){
    awaaa* tmpq=(awaaa*) tmpw;
    __int128 M1M2=tmpq->M1M2;
    __int128 M1M2M3=tmpq->M1M2M3;
    __int128 invM1_modM2=tmpq->invM1_modM2;
    __int128 invM1M2_modM3=tmpq->invM1M2_modM3;
    __int128 invM1M2M3_modM4=tmpq->invM1M2M3_modM4;
    long long p_=tmpq->p_;
    int K=tmpq->K;
    
    for (int j=K; j <= lim; j+=4) {
        __int128 x1=((ac[1][j]-ac[0][j])%mod[1]+mod[1])%mod[1];
        x1=x1*invM1_modM2%mod[1];
        __int128 r1=(x1*mod[0]%M1M2+ac[0][j])%M1M2;

        __int128 x2=((ac[2][j]-r1%mod[2])%mod[2]+mod[2])%mod[2];
        x2=x2*invM1M2_modM3%mod[2];
        __int128 r2=(x2*M1M2%M1M2M3+r1)%M1M2M3;

        __int128 x3=((ac[3][j]-r2%mod[3])%mod[3]+mod[3])%mod[3];
        x3=x3*invM1M2M3_modM4%mod[3];

        __int128 temp=x3%p_;
        temp=temp*(mod[0]%p_)%p_;
        temp=temp*(mod[1]%p_)%p_;
        temp=temp*(mod[2]%p_)%p_;
        temp=(temp+r2%p_)%p_;

        ab[j]=(long long)temp;
    }
    return NULL;
}




int main(int argc,char *argv[])
{
    
    // 保证输入的所有模数的原根均为 3,且模数都能表示为 a \times 4 ^ k+1 的形式
    // 输入模数分别为 7340033 104857601 469762049 1337006139375617
    // 第四个模数超过了整型表示范围,如果实现此模数意义下的多项式乘法需要修改框架
    // 对第四个模数的输入数据不做必要要求,如果要自行探索大模数 NTT,请在完成前三个模数的基础代码及优化后实现大模数 NTT
    // 输入文件共五个,第一个输入文件 n=4,其余四个文件分别对应四个模数,n=131072
    // 在实现快速数论变化前,后四个测试样例运行时间较久,推荐调试正确性时只使用输入文件 1
    long long test_begin=0;
    long long test_end=4;
    for(long long i=test_begin; i <= test_end; ++i){
        long double ans=0;
        long long n_,p_;

        memset(a,0,sizeof(a));
        memset(b,0,sizeof(b));
        memset(ac,0,sizeof(ac));
        memset(bc,0,sizeof(bc));
        //防止数据污染而加的清零

        fRead(a,b,&n_,&p_,i);
        memset(ab,0,sizeof(ab));
        auto Start=std::chrono::high_resolution_clock::now();
        // TODO : 将 poly_multiply 函数替换成你写的 ntt
     //   poly_multiply(a,b,ab,n_,p_);

        lim=1;
        k=0;
        while(lim<n_*2){//初始化幂次
            lim<<=1;
            ++k;
        }
        for(long long j=0;j<lim;j++){//对每个位置下标的二进制位操作 
            qwq[j]=((qwq[j>>1]>>1)|((j&1)<<(k-1)));
        }
 //       ntt(a,b,ab,n_,p_);//执行优化算法

 //下面是多模数ntt
    for(int j=0;j<4;++j){
        tmp4[j].j=j;
        tmp4[j].n_=n_;
        int errcode=pthread_create(&ths[j],NULL,NTT_parallel,&tmp4[j]);
    }//分别执行子任务
    for(int j=0;j<4;++j)    pthread_join(ths[j],NULL);
    __int128 M1M2=mod[0]*mod[1];
    __int128 M1M2M3=M1M2*mod[2];

    __int128 invM1_modM2=ksm(mod[0], mod[1]-2, mod[1]);
    __int128 invM1M2_modM3=ksm(M1M2%mod[2], mod[2]-2, mod[2]);
    __int128 invM1M2M3_modM4=ksm(M1M2M3%mod[3], mod[3]-2, mod[3]);

/*    for(int j=0;j<4;++j){
        tmpp4[j].invM1_modM2=invM1_modM2;
        tmpp4[j].invM1M2_modM3=invM1M2_modM3;
        tmpp4[j].invM1M2M3_modM4=invM1M2M3_modM4;
        tmpp4[j].M1M2=M1M2;
        tmpp4[j].M1M2M3=M1M2M3;
        tmpp4[j].K=j;
        tmpp4[j].p_=p_;
        int errcode=pthread_create(&ths[j],NULL,CRT_parallel,&tmpp4[j]);
    }
    for(int j=0;j<4;++j)    pthread_join(ths[j],NULL);*/
    for (int j=0; j <= lim; ++j) {
        __int128 x1=((ac[1][j]-ac[0][j])%mod[1]+mod[1])%mod[1];
        x1=x1*invM1_modM2%mod[1];
        __int128 r1=(x1*mod[0]%M1M2+ac[0][j])%M1M2;

        __int128 x2=((ac[2][j]-r1%mod[2])%mod[2]+mod[2])%mod[2];
        x2=x2*invM1M2_modM3%mod[2];
        __int128 r2=(x2*M1M2%M1M2M3+r1)%M1M2M3;

        __int128 x3=((ac[3][j]-r2%mod[3])%mod[3]+mod[3])%mod[3];
        x3=x3*invM1M2M3_modM4%mod[3];

        __int128 temp=x3%p_;
        temp=temp*(mod[0]%p_)%p_;
        temp=temp*(mod[1]%p_)%p_;
        temp=temp*(mod[2]%p_)%p_;
        temp=(temp+r2%p_)%p_;

        ab[j]=(long long)temp;
    }





        auto End=std::chrono::high_resolution_clock::now();
        std::chrono::duration<double,std::ratio<1,1000>>elapsed=End-Start;
        ans += elapsed.count();
        fCheck(ab,n_,i);
        std::cout<<"average latency for n="<<n_<<" p="<<p_<<" : "<<ans<<" (us) "<<std::endl;
        // 可以使用 fWrite 函数将 ab 的输出结果打印到 files 文件夹下
        // 禁止使用 cout 一次性输出大量文件内容
        fWrite(ab,n_,i);
    }
    return 0;
}
