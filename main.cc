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
#include<cstdint>

// 可以自行添加需要的头文件

void fRead(int *a, int *b, int *n, int *p, int input_id){
    // 数据输入函数
    std::string str1="/nttdata/";
    std::string str2=std::to_string(input_id);
    std::string strin=str1 + str2 + ".in";
    char data_path[strin.size() + 1];
    std::copy(strin.begin(), strin.end(), data_path);
    data_path[strin.size()]='\0';
    std::ifstream fin;
    fin.open(data_path, std::ios::in);
    fin>>*n>>*p;
    for (int i=0; i < *n; i++){
        fin>>a[i];
    }
    for (int i=0; i < *n; i++){   
        fin>>b[i];
    }
}

void fCheck(int *ab, int n, int input_id){
    // 判断多项式乘法结果是否正确
    std::string str1="/nttdata/";
    std::string str2=std::to_string(input_id);
    std::string strout=str1 + str2 + ".out";
    char data_path[strout.size() + 1];
    std::copy(strout.begin(), strout.end(), data_path);
    data_path[strout.size()]='\0';
    std::ifstream fin;
    fin.open(data_path, std::ios::in);
    for (int i=0; i < n * 2 - 1; i++){
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
    std::string str1="files/";
    std::string str2=std::to_string(input_id);
    std::string strout=str1 + str2 + ".out";
    char output_path[strout.size() + 1];
    std::copy(strout.begin(), strout.end(), output_path);
    output_path[strout.size()]='\0';
    std::ofstream fout;
    fout.open(output_path, std::ios::out);
    for (int i=0; i < n * 2 - 1; i++){
        fout<<ab[i]<<'\n';
    }
}

void poly_multiply(int *a, int *b, int *ab, int n, int p){
    for(int i=0; i < n; ++i){
        for(int j=0; j < n; ++j){
            ab[i+j]=(1ll * a[i] * b[j] % p + ab[i+j]) % p;
        }
    }
}






uint32x4_t mc(uint32x4_t A,uint32x4_t B,int mod){
    uint32x2_t al = vget_low_u32(A);
    uint32x2_t ah = vget_high_u32(A); 
    uint32x2_t bl = vget_low_u32(B);
    uint32x2_t bh = vget_high_u32(B); 

    uint64x2_t awa1=vmull_u32(al,bl);
    uint64x2_t awa2=vmull_u32(ah,bh);
    int tmpa=(uint32_t)(vgetq_lane_u64(awa1,0)%mod);
    int tmpb=(uint32_t)(vgetq_lane_u64(awa1,1)%mod);
    int tmpc=(uint32_t)(vgetq_lane_u64(awa2,0)%mod);
    int tmpd=(uint32_t)(vgetq_lane_u64(awa2,1)%mod);
    uint32x4_t tmpans={tmpa,tmpb,tmpc,tmpd};
    return tmpans;
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

  /*      for(int j=0;j<lim;j+=(mid<<1)){//主体
            int w=1;
            for(int k=0;k<mid;++k,w=1ll*w*wn%p){
                int x=a[j+k];
                int y=1ll*a[j+k+mid]*w%p;
                a[j+k]=(x+y)%p;
                a[j+k+mid]=(x-y+p)%p;
            }
        }*/
 
        uint32x4_t modv=vdupq_n_u32(p);//{p,p,p,p}，方便内部调用
        for(int j=0;j<lim;j+=(mid<<1)){//向量化的主体
            int w=1;
            int k=0;
            for(;k+3<mid;k+=4){
                int32x4_t x={a[j+k],a[j+k+1],a[j+k+2],a[j+k+3]};//加载
                int32x4_t b={a[j+k+mid],a[j+k+mid+1],a[j+k+mid+2],a[j+k+mid+3]};
    
                int w1=w;
                int w2=1ll*w1*wn%p;
                int w3=1ll*w2*wn%p;
                int w4=1ll*w3*wn%p;
                uint32x4_t wv={w1,w2,w3,w4};//预处理4个一组的单位根
                w=1ll*w4*wn%p;//将幅角转4倍单位根

                uint32x4_t y=mc(vreinterpretq_u32_s32(b),wv,p);//先做乘法
    
                int32x4_t sum=vaddq_s32(x,vreinterpretq_s32_u32(y));//x+y->a[j+k]
                int32x4_t diff=vsubq_s32(x,vreinterpretq_s32_u32(y));//x-y->a[j+k+mid]

            //从这里开始
                uint32x4_t mask1=vcgeq_s32(sum,vreinterpretq_s32_u32(modv));//比较sum的每一个元素和p
                sum=vreinterpretq_s32_u32(vsubq_u32(vreinterpretq_u32_s32(sum),vandq_u32(mask1,modv)));//如果大于p（对应位置为0xFFFFFFFF），就-p
    
                uint32x4_t mask2=vcltq_s32(diff,vdupq_n_s32(0));//clt和cge相反
                diff=vreinterpretq_s32_u32(vaddq_u32(vreinterpretq_u32_s32(diff),vandq_u32(mask2,modv)));//类似
            //到这里为止
            //注意到我们的sum=x+y,又有diff=x-y,而x和y都小于p，于是结果在-p到2p-1间，只需要判断溢出再做一次加减就行
                
                vst1q_s32(a+j+k, sum);
                vst1q_s32(a+j+k+mid, diff);//存回内存
            }
            for(;k<mid;++k,w=1ll*w*wn%p){//这里处理最后剩下的部分
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

/*         int i=0;
        for(i=0;i+3<lim;i+=4){
            uint32x4_t tmpp={a[i],a[i+1],a[i+2],a[i+3]};
            uint32x4_t vec_awa=vdupq_n_u32(awa);
            uint32x4_t tmpans=vmulq_u32(tmpp,vec_awa);
            a[i]=vgetq_lane_u32(tmpans,0);
            a[i+1]=vgetq_lane_u32(tmpans,1);
            a[i+2]=vgetq_lane_u32(tmpans,2);
            a[i+3]=vgetq_lane_u32(tmpans,3);
        }
        for(;i<lim;++i){
            a[i]=1ll*a[i]*awa%p;
        }
*/
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
    // 输入文件共五个, 第一个输入文件 n=4, 其余四个文件分别对应四个模数, n=131072
    // 在实现快速数论变化前, 后四个测试样例运行时间较久, 推荐调试正确性时只使用输入文件 1
    int test_begin=0;
    int test_end=1;
    for(int i=test_begin; i <= test_end; ++i){
        long double ans=0;
        int n_, p_;

        memset(a,0,sizeof(a));
        memset(b,0,sizeof(b));//防止数据污染而加的清零

        fRead(a, b, &n_, &p_, i);
        memset(ab,0,sizeof(ab));
        auto Start=std::chrono::high_resolution_clock::now();
        // TODO : 将 poly_multiply 函数替换成你写的 ntt
    //    poly_multiply(a, b, ab, n_, p_);


        ntt(a,b,ab,n_,p_);//执行优化算法


        auto End=std::chrono::high_resolution_clock::now();
        std::chrono::duration<double,std::ratio<1,1000>>elapsed=End - Start;
        ans += elapsed.count();
        fCheck(ab, n_, i);
        std::cout<<"average latency for n="<<n_<<" p="<<p_<<" : "<<ans<<" (us) "<<std::endl;
        // 可以使用 fWrite 函数将 ab 的输出结果打印到 files 文件夹下
        // 禁止使用 cout 一次性输出大量文件内容
        fWrite(ab, n_, i);
    }
    return 0;
}
