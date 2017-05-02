//*H29年度/dsp2-1/5J06*//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define MAX_SIZE 16384
//#define M_PI 3.14159265359

typedef struct {
    double re; //実部
    double im; //虚部
} complex_t;

FILE * openFile();
FILE * openFile_2();
complex_t comp_add(complex_t in1, complex_t in2);
complex_t comp_sub(complex_t in1, complex_t in2);
complex_t comp_multiply(complex_t in1, complex_t in2);
complex_t comp_division(complex_t in1, complex_t in2);
complex_t comp_conjugate(complex_t in);
void twiddle_factor(complex_t *wnk, int n, double a);
void bit_r(int *bit, int N);
void convert_r(complex_t *in, int *bit, int N);
void fft(complex_t *in, complex_t *wnk, int N);
void power_spectral(complex_t *in1, complex_t *in2, complex_t *ps, int N);
int power(int n);
void cc(double *data[],int size1);
void correlation_function(complex_t *x1, complex_t *x2, complex_t *ps, int size, int N);


int main(){

    int i, j, N, size1, size2, p;
    complex_t x1[MAX_SIZE], x2[MAX_SIZE], ps[MAX_SIZE];
    FILE *fp1, *fp2, *fq;

    printf("\nH29 課題1 出席番号06\n");
    printf("~使い方~\n");
    printf("・読み込む2つのテキストファイルを同じディレクトリに置いて下さい\n");
    printf("・自己相関または相互相関を求めテキストファイルを出力します\n");
    printf("~実行環境~\n");
    printf("OS : OS El Capitan (Version 10.11.6)\nプロセッサ : 1.7 GHz Intel Core i7\nメモリ : 8 GB 1600 MHz DDR3\n");

    printf("1つ目のファイル名を入力: ");
    fp1 = openFile();
    printf("2つ目のファイル名を入力: ");
    fp2 = openFile();
    fq = openFile_2();

    size1 = 0;
    size2 = 0;
    while((j = getc(fp1)) != EOF) if(j == '\n') size1++;
    while((j = getc(fp2)) != EOF) if(j == '\n') size2++;
    if(size1 != size2){
        printf("error\n");
        return 0;
    }

    //点数決定
    p = power(size1*2);
    N = pow(2, p);

    // 構造体初期化
    for(i=0;i<MAX_SIZE;i++){
        x1[i].re = 0;
        x1[i].im = 0;
        x2[i].re = 0;
        x2[i].im = 0;
        ps[i].re = 0;
        ps[i].im = 0;
    }


    //読み込みデータを構造体に代入
    fseek(fp1,  0L, SEEK_SET);
    fseek(fp2,  0L, SEEK_SET);

    for(i=0;i<size1;i++) fscanf(fp1,"%lf", &x1[i].re);
    for(i=0;i<size2;i++) fscanf(fp2,"%lf", &x2[i].re);

    fclose(fp1);
    fclose(fp2);

    correlation_function(x1,x2,ps,size1,N);
    for(i=0;i<size1;i++){
        fprintf(fq, "%f\n", ps[i].re);
    }
    fclose(fq);

    return 0;
}

FILE* openFile(){
    char file[128];
    // printf("Enter the input file name: ");
    scanf("%s",file);

    FILE *fp = fopen(file,"r");
    if (fp==NULL) {
        printf("can't open the input file\n");
        exit(0);
    }
    return fp;
}

FILE* openFile_2(){
    char file[128];
    printf("出力ファイル名を入力: ");
    scanf("%s",file);

    FILE *fp = fopen(file,"w");
    if (fp==NULL) {
        printf("can't open the output file\n");
        exit(0);
    }
    return fp;
}

complex_t comp_add(complex_t in1, complex_t in2){
    complex_t res;
    res.re = in1.re + in2.re;
    res.im = in1.im + in2.im;
    return res;
}

complex_t comp_sub(complex_t in1, complex_t in2){
    complex_t res;
    res.re = in1.re - in2.re;
    res.im = in1.im - in2.im;
    return res;
}

complex_t comp_multiply(complex_t in1, complex_t in2){
    complex_t res;
    res.re = (in1.re * in2.re) - (in1.im * in2.im);
    res.im = (in1.re * in2.im) + (in1.im * in2.re);
    return res;
}

complex_t comp_division(complex_t in1, complex_t in2){
    complex_t res;
    complex_t in2_copy = comp_conjugate(in2);
    res = comp_multiply(in1, in2_copy);
    res.re /= comp_multiply(in2, in2_copy).re;
    res.im /= comp_multiply(in2, in2_copy).re;
    return res;
}

complex_t comp_conjugate(complex_t in){
    complex_t res;
    res.re = in.re;
    res.im = -1 * in.im;
    return res;
}

void twiddle_factor(complex_t *wnk, int n, double a){
    int i;
    for(i=0;i<n;i++){
        wnk[i].re = cos(a*2*M_PI/n*i);//回転子の値
        wnk[i].im = sin(a*-2*M_PI/n*i);
    }
}

void bit_r(int *bit, int N){
    int r,i,j;
    r = (int)(log(N)/log(2.0) + 0.5);

    for(i=0;i<N;i++){
        bit[i] = 0;
        for(j=0;j<r;j++){
            bit[i] += ((i >> j) & 1) << (r-j-1);
        }
    }
}

void convert_r(complex_t *in, int *bit, int N){
    int i;
    complex_t comp[N];
    for(i=0;i<N;i++){
        comp[i].re = in[i].re;
        comp[i].im = in[i].im;
    }
    for(i=0;i<N;i++){
        in[i].re = comp[bit[i]].re;
        in[i].im = comp[bit[i]].im;
    }
}

void fft(complex_t *in, complex_t *wnk, int N){
    int r_big = 1, r_sma = N / 2; //初期値を与えておく
    int i, j, k, in1, in2, nk, r;
    complex_t dummy;
    clock_t start, end;

    r = (int)(log(N) / log(2.0) + 0.5);

    start = clock();
    for(i=0;i<r;i++){// ビット回繰り返し（段数回）
        for(j=0;j<r_big;j++){//繰り返しが段々増える
            for(k=0;k<r_sma;k++){//繰り返しが段々減る：(j+k)の合計=N/2回の繰り返し
                in1 = r_big * 2 * k + j;    //バタ上段の入力順番
                in2 = in1 + r_big;      //バタ下段の入力順番
                nk = j * r_sma;     //回転子の番号

                dummy = comp_multiply(in[in2], wnk[nk]);//バタ入力下段×重みWnk
                in[in2] = comp_sub(in[in1], dummy);//バタ演算の下段出力(そのまま次段の入力)
                in[in1] = comp_add(in[in1], dummy); //バタ演算の上段出力(そのまま次段の入力)
            }
        }
        r_big *= 2; // 2^i      :1,2,4...N/2まで増加
        r_sma /= 2; // 2^(r-1-i):N/2...4.2.1まで減少
    }
    end = clock();
    printf( "処理時間:%lf[s]\n", (double)(end - start)/CLOCKS_PER_SEC);
}

void power_spectral(complex_t *in1, complex_t *in2, complex_t *ps, int N){
    int i;
    for(i=0;i<N;i++) ps[i].re = (in1[i].re * in2[i].re) + (in1[i].im * in2[i].im);
}

int power(int n){
    return (int)ceil(log10(n) / log10(2));
}

void correlation_function(complex_t *x1, complex_t *x2, complex_t *ps, int size, int N){
    int i;
    int bit[MAX_SIZE];
    complex_t r[MAX_SIZE];

    for(i=0;i<MAX_SIZE;i++){
        r[i].re = 0;
        r[i].im = 0;
    }
    twiddle_factor(r, N, 1); //回転子生成
    bit_r(bit,N); //ビットリバーサル

    //データ入れ替え
    convert_r(x1,bit,N);
    convert_r(x2,bit,N);

    //FFT
    fft(x1,r,N);
    fft(x2,r,N);

    power_spectral(x1,x2,ps,N);

    // IFFT
    twiddle_factor(r, N, -1); //回転子生成
    bit_r(bit,N); //ビットリバーサル
    convert_r(ps,bit,N);
    fft(ps,r,N);
    for(i=0;i<N;i++) {
        ps[i].re /= N * size;
        ps[i].im /= N * size;
    }
}