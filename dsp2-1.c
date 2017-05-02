//H29 dsp2-1 5J04
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
// コメントgithub

//#define debug printf("Check\n");

/*************構造体、グローバル変数*************/
typedef struct
{
	double re;//実部
	double im;//虚部
}comp_t;

/****************プロトタイプ宣言****************/
void print_index();							//説明文を表示
FILE* open_rfile();							//ファイルを読み込み形式でオープン
FILE* open_wfile();							//ファイルを書き込み方式でオープン
int count_data(FILE*);						//データ数を返す
int* int_calloc(int);						//int個分のint型配列を動的生成
double* double_calloc(int);					//int個分のdouble型配列を動的生成
comp_t* comp_calloc(int);					//int個分のcomp_t型配列を動的生成
comp_t comp_add(comp_t,comp_t);				//複素数の和を返す
comp_t comp_sub(comp_t,comp_t);				//複素数の差を返す
comp_t comp_multi(comp_t,comp_t);			//複素数の積を返す
comp_t comp_div(comp_t,comp_t);				//複素数の商を返す
comp_t conjugate(comp_t);					//共役複素数を返す
void twid(comp_t*,double,int);				//回転子を計算
void reverse_b(int*,int);					//ビットリバーサルの順番を計算
void reverse_c(comp_t*,int*,int);			//ビットリバーサルの並び替え
void fft(comp_t*,int);						//FFTを実行
void ifft(comp_t*,int);						//IFFTを実行
comp_t* pow_spectrum(comp_t*,comp_t*,int);	//パワースペクトルを求める

/**********************MAIN**********************/
int main()
{
	FILE *fr1,*fr2,*fw;
	int i,j;
	int N;							//FFTの点数
	int quant,size;
	comp_t *corre1 = NULL,*corre2 = NULL,*power = NULL;

	// 説明文を表示
	print_index();

	/***********相関係数を求める***********/
	// ファイルオープン
	fr1 = open_rfile();
	fr2 = open_rfile();

	// データ数をカウント，fr1とfr2のデータ数が異なる場合はエラー
	quant = count_data(fr1);
	if(quant != count_data(fr2) $ )
	{
		printf("failed.");
		exit(1);
	}

	// データ数を元に，2のべき乗個の配列を作成
	N = pow(2,(int)ceil(log10(quant*2) / log10(2)));
	corre1 = comp_calloc(N);
	corre2 = comp_calloc(N);

	// データを構造体に読み込む
	for(i=0; i<quant; ++i)
	{
		fscanf(fr1,"%lf",&corre1[i].re);
		fscanf(fr2,"%lf",&corre2[i].re);
	}

	// FFTを行う
	fft(corre1,N);
	fft(corre2,N);

	// パワースペクトルを求める
	power = comp_calloc(N);
	power = pow_spectrum(corre1,corre2,N);

	// IFFT
	ifft(power,N);

	printf("\ncompleted\n");
	fw = open_wfile();

	for(i=0; i<quant; ++i) fprintf(fw,"%lf\n",power[i].re / quant);

	fclose(fr1);
	fclose(fr2);
	fclose(fw);
	free(corre1);
	free(corre2);
	free(power);
}

/**********************関数**********************/
void print_index()
{
	printf("H29 dsp2-1 5J04\n");
	printf("入力した二つのファイル内データの相関係数を求めます。\n");
	printf("相関係数を求めたいデータが入った2つのファイルの名前を入力してください。\n");
}

FILE* open_rfile()
{
	FILE *fr;
	char name[64];
	
	printf("Please tell me Readfile name．\n> ");
	scanf("%s",name);

	sprintf(name,"%s.txt",name);
	
	/*ファイルの読み込み開始*/
	if((fr = fopen(name,"r")) == NULL)
	{
		perror(name);
		exit(1);
	}
	return fr;
}

FILE* open_wfile()
{
	FILE *fw;
	char name[64];

	printf("Please tell me Writefile name．\n> ");
	scanf("%s",name);

	sprintf(name,"%s.csv",name);

	if((fw = fopen(name,"w")) == NULL)
	{
		perror(name);
		exit(1);
	}
	return fw;
}

int count_data(FILE *fr)
{
	int i = 0;
	double a;
	while(fscanf(fr,"%lf",&a) != EOF) i++;
	fseek(fr,0L,SEEK_SET);
	return i;
}

int* int_calloc(int size)
{
	int *array = NULL;
	array = (int*)calloc(size,sizeof(int));
	if(array == NULL)
	{
	    printf("Can't allocate memory.\n");
	    exit(1);
	}
	return array;
}

double* double_calloc(int size)
{
	double *array = NULL;
	array = (double*)calloc(size,sizeof(double));
	if(array == NULL)
	{
	    printf("Can't allocate memory.\n");
	    exit(1);
	}
	return array;
}

comp_t* comp_calloc(int size)
{
	comp_t *array = NULL;
	array = (comp_t*)calloc(size,sizeof(comp_t));
	if(array == NULL)
	{
	    printf("Can't allocate memory.\n");
	    exit(1);
	}
	return array;
}

comp_t comp_add(comp_t x1,comp_t x2)
{
	comp_t ans;
	ans.re=x1.re+x2.re;
	ans.im=x1.im+x2.im;
	return ans;
}

comp_t comp_sub(comp_t x1,comp_t x2)
{
	comp_t ans;
	ans.re=x1.re-x2.re;
	ans.im=x1.im-x2.im;
	return ans;
}

comp_t comp_multi(comp_t x1,comp_t x2)
{
	comp_t ans;
	ans.re=(x1.re*x2.re)-(x1.im*x2.im);
	ans.im=(x1.im*x2.re)+(x1.re*x2.im);
	return ans;
}

comp_t comp_div(comp_t x1,comp_t x2)
{
	comp_t ans;
	ans.re=(x1.re*x2.re+x1.im*x2.im)/(x2.re*x2.re+x2.im*x2.im);
	ans.im=(x1.im*x2.re-x1.re*x2.im)/(x2.re*x2.re+x2.im*x2.im);
	return ans;
}

comp_t conjugate(comp_t x)
{
	comp_t comp;
	comp.re=x.re;
	comp.im=-x.im;
	return comp;
}

void twid(comp_t *fact,double a,int N)
{
	int i; 
	for(i=0; i<N; ++i)
	{
		fact[i].re = cos(2*M_PI*i/N*a);
		fact[i].im = sin(-2*M_PI*i/N*a);
	}
}

void reverse_b(int *bit,int N)
{
	int i,j;
	int r = (int)(log(N) / log(2.0));;
	
	for(i=0; i<N; ++i)
	{
		bit[i] = 0;
		for(j=0; j<r; ++j)
		{
			bit[i] += ((i >> j) & 1) << (r - j - 1);
		}
	}
}

void reverse_c(comp_t *in, int *bit,int N)
{
	int i;
	comp_t *comp = NULL;
	comp=(comp_t*)malloc(sizeof(comp_t)*N);
	for(i=0; i<N; ++i) comp[i] = in[i];
	for(i=0; i<N; ++i) in[i] = comp[bit[i]];
	free(comp);
}

void fft(comp_t *in,int N)
{
	int r_big=1,r_sma=N/2;
	int i,j,k;
	int in1,in2,nk;
	int r=(int)(log(N)/log(2.0));
	int *bit = NULL;
	comp_t *fact = NULL;

	bit = int_calloc(N);
	fact = comp_calloc(N);

	//ビットリバーサル
	reverse_b(bit,N);
	reverse_c(in,bit,N);

	//回転子を計算
	twid(fact,1,N);

	comp_t dummy;
	for(i=0; i<r; ++i)
	{
		for(j=0; j<r_big; ++j)
		{
			for(k=0; k<r_sma; ++k)
			{
				in1=r_big*2*k+j;
				in2=in1+r_big;
				nk=j*r_sma;
				dummy = comp_multi(in[in2], fact[nk]);
				in[in2] = comp_sub(in[in1], dummy);
				in[in1] = comp_add(in[in1], dummy);
			}
		}
		r_big*=2;
		r_sma/=2;
	}
	free(bit);
	free(fact);
}

void ifft(comp_t *in,int N)
{
	int r_big=1,r_sma=N/2;
	int i,j,k;
	int in1,in2,nk;
	int r=(int)(log(N)/log(2.0));
	int *bit = NULL;
	comp_t *fact = NULL;

	bit = int_calloc(N);
	fact = comp_calloc(N);

	//ビットリバーサル
	reverse_b(bit,N);
	reverse_c(in,bit,N);

	//回転子を計算
	twid(fact,-1,N);

	comp_t dummy;
	for(i=0; i<r; ++i)
	{
		for(j=0; j<r_big; ++j)
		{
			for(k=0; k<r_sma; ++k)
			{
				in1=r_big*2*k+j;
				in2=in1+r_big;
				nk=j*r_sma;
				dummy = comp_multi(in[in2], fact[nk]);
				in[in2] = comp_sub(in[in1], dummy);
				in[in1] = comp_add(in[in1], dummy);
			}
		}
		r_big*=2;
		r_sma/=2;
	}
	for(i=0; i<N; ++i) in[i].re /= N;
	free(bit);
	free(fact);
}

comp_t* pow_spectrum(comp_t *in1,comp_t *in2,int N)
{
	comp_t *power = NULL;
	power = comp_calloc(N);
	for(int i=0; i<N; ++i) power[i] = comp_multi(conjugate(in1[i]),in2[i]);
	return power;
}