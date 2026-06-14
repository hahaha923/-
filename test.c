#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define RAND_MAX_C 32767
#define SAMPLE_NUM 1000000   // 总采样样本量
#define N_TEST 10            // 测试区间 [0,N_TEST-1]
#define PI_SAMPLE 5000000    // 蒙特卡洛π采样量
#define ALPHA 0.05           // 显著性水平

// 全局存储原始随机序列，所有工具共用同一组数据
int* rand_seq = NULL;

// ===================== 工具1：频数+概率误差统计 =====================
void freq_analysis(int N, long long *cnt, long long total)
{
    double p_theo = 1.0 / N;
    printf("==========【工具1：频数概率误差分析】==========\n");
    printf("理论单值概率 p=%.6f\n", p_theo);
    double max_err = 0.0, mse = 0.0;
    for (int i = 0; i < N; i++)
    {
        double p_obs = (double)cnt[i] / total;
        double err = fabs(p_obs - p_theo);
        if (err > max_err) max_err = err;
        mse += pow(err, 2);
        printf("数字%2d | 计数=%6lld | 实测概率=%.6f | 误差=%.6f\n",
               i, cnt[i], p_obs, err);
    }
    mse /= N;
    printf("最大概率绝对误差：%.6f | 均方误差MSE：%.8f\n\n", max_err, mse);
}

// ===================== 工具2：卡方均匀性检验 =====================
int chi2_test(int N, long long *cnt, long long total)
{
    printf("==========【工具2：卡方均匀性检验】==========\n");
    double E = (double)total / N;
    double chi2 = 0.0;
    for (int i = 0; i < N; i++)
    {
        chi2 += pow((double)cnt[i] - E, 2) / E;
    }
    double chi2_critical[] = {3.841, 5.991, 7.815, 9.488, 11.070,
                              12.592, 14.067, 15.507, 16.919, 18.307};
    int df = N - 1;
    double cri = (df < 10) ? chi2_critical[df - 1] : 20.483;
    printf("卡方统计量=%.4f | 自由度df=%d | 临界值=%.4f\n", chi2, df, cri);
    if (chi2 < cri)
    {
        printf("检验结果：通过均匀分布假设（P>0.05）\n\n");
        return 1;
    }
    else
    {
        printf("检验结果：拒绝均匀分布假设（存在显著分布偏差）\n\n");
        return 0;
    }
}

// ===================== 工具3：Gap重复间隔分布检验（读取预存序列） =====================
void gap_analysis(int N, long long total)
{
    printf("==========【工具3：Gap重复间隔分布检验】==========\n");
    long long last_pos[N];
    long long gap_cnt[20] = {0};
    memset(last_pos, -1, sizeof(last_pos));
    for (long long pos = 0; pos < total; pos++)
    {
        int x = rand_seq[pos]; // 读取预存序列，不再调用rand()
        if (last_pos[x] != -1)
        {
            long long g = pos - last_pos[x] - 1;
            if (g < 20) gap_cnt[g]++;
        }
        last_pos[x] = pos;
    }
    double p_g;
    for (int g = 0; g < 10; g++)
    {
        p_g = pow((N - 1.0) / N, g) * (1.0 / N);
        printf("间隔g=%d | 观测频次=%6lld | 理论概率=%.6f\n",
               g, gap_cnt[g], p_g);
    }
    printf("\n");
}

// ===================== 工具4：滞后1阶自相关检验 =====================
void auto_corr_test(int N, long long total)
{
    printf("==========【工具4：滞后1阶自相关检验】==========\n");
    double mean = 0.0, var = 0.0, cov = 0.0;
    int prev = rand_seq[0];
    mean += prev;
    for (long long i = 1; i < total; i++)
    {
        int curr = rand_seq[i]; // 读取预存数组，无新rand调用
        mean += curr;
        cov += (prev - N/2.0) * (curr - N/2.0);
        var += pow(curr - N/2.0, 2);
        prev = curr;
    }
    mean /= total;
    var /= total;
    double r = cov / (total * var);
    printf("序列均值=%.4f | 方差=%.4f | 滞后1阶自相关系数r=%.6f\n", mean, var, r);
    if (fabs(r) > 0.01)
        printf("结论：自相关系数偏离0，序列存在线性相关性缺陷\n\n");
    else
        printf("结论：无明显一阶自相关\n\n");
}

// ===================== 工具5：扑克随机性检验 =====================
// 分组长度k=5，每组5个[0,N-1]整数，统计重复模式
int poker_test(int N, long long total, int* seq)
{
    printf("==========【工具5：扑克随机性检验】==========\n");
    const int k = 5;  // 每组取5个随机数
    long long group_num = total / k;
    if(group_num < 1000)
    {
        printf("样本量不足，无法执行扑克检验\n\n");
        return 0;
    }

    // 5分组下的5种组合类型：0=全不同,1=一对,2=两对,3=三条,4=四/五条
    long long count_type[5] = {0};

    for(long long g=0; g<group_num; g++)
    {
        int bucket[10] = {0};
        // 取出一组5个数
        for(int i=0; i<k; i++)
        {
            int val = seq[g*k + i];
            bucket[val]++;
        }
        // 统计组内重复数字的出现次数
        int repeat[6] = {0};
        for(int i=0; i<N; i++)
        {
            if(bucket[i]>0) repeat[bucket[i]]++;
        }
        // 分类判定
        if(repeat[5]==1) count_type[4]++;    // 五条相同
        else if(repeat[4]==1) count_type[4]++; // 四条相同
        else if(repeat[3]==1) count_type[3]++; // 三条
        else if(repeat[2]==2) count_type[2]++; // 两对
        else if(repeat[2]==1) count_type[1]++; // 一对
        else count_type[0]++;                  // 全部数字不重复
    }

    // N=10、k=5时，理论概率（预计算）
    double prob[5] = {0.3024, 0.5040, 0.1080, 0.0720, 0.0136};
    double chi2 = 0.0;
    for(int t=0; t<5; t++)
    {
        double exp = group_num * prob[t];
        chi2 += pow(count_type[t] - exp, 2) / exp;
    }

    double critical = 9.488; // df=4，α=0.05 卡方临界值
    printf("总组数=%lld | 卡方统计量=%.4f | df=4 临界值=%.4f\n", group_num, chi2, critical);
    printf("全不重复:%lld 一对:%lld 两对:%lld 三条:%lld 四/五条:%lld\n",
           count_type[0],count_type[1],count_type[2],count_type[3],count_type[4]);
    if(chi2 < critical)
    {
        printf("扑克检验通过：序列无异常重复扎堆，随机性合格\n\n");
        return 1;
    }
    else
    {
        printf("扑克检验不通过：数字重复分布偏离随机预期\n\n");
        return 0;
    }
}
// ===================== 无偏均匀随机数（拒绝采样） =====================
int uniform_unbias(int N)
{
    int limit = RAND_MAX_C - (RAND_MAX_C % N);
    int x;
    do
    {
        x = rand();
    } while (x >= limit);
    return x % N;
}

// ===================== Box-Muller正态分布随机数 =====================
double normal_rand(double mu, double sigma)
{
    static int flag = 0;
    static double z1;
    if (flag == 0)
    {
        double u1 = (rand() + 1.0) / (RAND_MAX_C + 2.0);
        double u2 = (rand() + 1.0) / (RAND_MAX_C + 2.0);
        double z0 = sqrt(-2 * log(u1)) * cos(2 * M_PI * u2);
        z1 = sqrt(-2 * log(u1)) * sin(2 * M_PI * u2);
        flag = 1;
        return mu + sigma * z0;
    }
    else
    {
        flag = 0;
        return mu + sigma * z1;
    }
}

// ===================== 蒙特卡洛法求π =====================
void monte_carlo_pi()
{
    printf("==========【蒙特卡洛求π对比实验】==========\n");
    long long inside1 = 0;
    for (long long i = 0; i < PI_SAMPLE; i++)
    {
        double x = (rand() % 10000) / 1000.0;
        double y = (rand() % 10000) / 1000.0;
        x /= 10.0; y /= 10.0;
        if (x*x + y*y <= 1.0) inside1++;
    }
    double pi1 = 4.0 * inside1 / PI_SAMPLE;

    long long inside2 = 0;
    for (long long i = 0; i < PI_SAMPLE; i++)
    {
        double x = uniform_unbias(10000) / 10000.0;
        double y = uniform_unbias(10000) / 10000.0;
        if (x*x + y*y <= 1.0) inside2++;
    }
    double pi2 = 4.0 * inside2 / PI_SAMPLE;

    printf("样本总量：%lld\n", PI_SAMPLE);
    printf("原生rand()%%N计算π = %.8f | 与真值误差=%.8f\n", pi1, fabs(pi1 - M_PI));
    printf("无偏均匀随机数计算π = %.8f | 与真值误差=%.8f\n\n", pi2, fabs(pi2 - M_PI));
}

// 导出二进制比特流（用于中经国昱密评工具箱）
void export_bitstream(int N, long long num)
{
    FILE *fp = fopen("rand_bit.bin", "wb");
    unsigned char buf = 0;
    int bit_cnt = 0;
    for(long long i=0; i<num; i++)
    {
        int x = rand_seq[i];
        for(int b=0;b<4;b++)
        {
            buf <<= 1;
            buf |= (x >> b) & 1;
            bit_cnt++;
            if(bit_cnt == 8)
            {
                fputc(buf, fp);
                buf = 0; bit_cnt = 0;
            }
        }
    }
    if(bit_cnt > 0) fputc(buf << (8-bit_cnt), fp);
    fclose(fp);
    printf("比特流文件 rand_bit.bin 导出完成\n\n");
}

// ===================== 主函数：预先生成统一序列 =====================
int main()
{
    srand((unsigned)time(NULL));
    int N = N_TEST;
    long long *count = calloc(N, sizeof(long long));

    // 步骤1：一次性生成全部随机序列存入数组，所有工具共用
    rand_seq = (int*)malloc(sizeof(int) * SAMPLE_NUM);
    for (long long i = 0; i < SAMPLE_NUM; i++)
    {
        int val = rand() % N;
        rand_seq[i] = val;
        count[val]++; // 同步统计频数
    }

    printf("===== 随机数分布综合测试实验 N=%d 总采样%lld（全工具共用同一序列）=====\n\n", N, SAMPLE_NUM);

    // 执行5套测评工具，全部读取预存rand_seq
    freq_analysis(N, count, SAMPLE_NUM);
    chi2_test(N, count, SAMPLE_NUM);
    gap_analysis(N, SAMPLE_NUM);
    auto_corr_test(N, SAMPLE_NUM);
    poker_test(N, SAMPLE_NUM, rand_seq);
    // 无偏均匀随机数演示
    printf("==========【无偏均匀随机数示例】==========\n");
    printf("连续10个无偏[0,%d)随机数：", N-1);
    for (int i = 0; i < 10; i++)
        printf("%d ", uniform_unbias(N));
    printf("\n\n");

    // 正态随机数演示
    printf("==========【Box-Muller正态随机数 N(0,1)示例】==========\n");
    printf("连续10个标准正态样本：");
    for (int i = 0; i < 10; i++)
        printf("%.3f ", normal_rand(0, 1));
    printf("\n\n");

    // 蒙特卡洛求π
    monte_carlo_pi();

    // 导出密评比特流文件
    export_bitstream(N, 300000);

    // 释放内存
    free(count);
    free(rand_seq);
    rand_seq = NULL;
    return 0;
}