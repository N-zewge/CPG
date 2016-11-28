#include	<stdio.h>
#include	<stdlib.h>
#include	<math.h>
#include	<time.h>
#include "Neural_Oscillator.h"   // added
 
#define		MaxGens				50
#define		NumParameter		27 // modified.
#define		ParentsOffspring	100
#define		Parents				50
#define		True				1
#define		False				0
////////////////////////////////////////
#define		iteration			600
//#define		Neural_N			3
#define     PI                  3.141592
///////////////////////////////////////
  

double		compare[ParentsOffspring],
			score[ParentsOffspring],
			save_score[ParentsOffspring],
			vector[ParentsOffspring][NumParameter],
			save_vector[ParentsOffspring][NumParameter],
			smallest;
int			wins[ParentsOffspring], K, maximum,
			count, flag, BestSolution, generation;
FILE		*outfile, *outfile2, *sinfile;

void InitialVectors();
void ScoreAllVectors();
void GenerateOffspring();
void CountWins();
void RankPopulation();
void Statistics();
double rndnor(double mean,double var);


double ue[3] = {-0.0042, -0.0042, -0.0042};
double uf[3] = {0.8372, 0.8372, 0.8372};
double ve[3] = {0.1834, 0.1834, 0.1834};
double vf[3] = {0.6501, 0.6501, 0.6501}; // hw에 있는것 그대로 사용

double ye[3] = {0, 0, 0};
double yf[3] = {0, 0, 0};


double wee[9];					//입력하는것
double wef[9];
double wfe[9];
double wff[9];

double score_result[30000] = {0};



double uinput;			
double tc[2];  


double sum1, sum2, sum3;
double	s1, s2, s3;




//참고: wef[i * Neural_N + j]
//////////////////////

void main(void)
{

	int i, k;
	double *result = 0;

	FILE *fout1, *fout2, *fout3;
	fopen_s(&fout1, "output1.txt", "w");
	fopen_s(&fout2, "output2.txt", "w");
	fopen_s(&fout3, "output3.txt", "w");
	fopen_s(&sinfile, "sin_output.txt", "w");
	count = 0;

	if ((outfile = fopen("mean.txt","w"))==NULL) {
		exit(1);
	}
	if ((outfile2 = fopen("best.txt", "w")) == NULL) {
		exit(1);
	}
	fprintf(outfile,"Outfile for mean/var of population score...\n");
	fprintf(outfile2,"Outfile for best score at each generation...\n");

	K = Parents;
	maximum = Parents;
	generation = 0;
	srand( (unsigned)time(NULL) );
	InitialVectors();
	ScoreAllVectors();

	do{
		GenerateOffspring();
		ScoreAllVectors();
		CountWins();
		RankPopulation();
		Statistics();
		maximum = Parents;
		generation++;
	} while ( generation <= MaxGens );

	for(i=0; i< NumParameter ; i++){
		//printf("X[%d] = %lf \n", i, vector[BestSolution][i]);
		fprintf(outfile2,"X[%d] = %lf \n", i, vector[BestSolution][i]);


	}


	
		wee[1] = vector[BestSolution][3];	wee[2] = vector[BestSolution][4];
		wee[3] = vector[BestSolution][5];	wee[5] = vector[BestSolution][6];
		wee[6] = vector[BestSolution][7];	wee[7] = vector[BestSolution][8];

		wef[1] = vector[BestSolution][9];	wef[2] = vector[BestSolution][10];
		wef[3] = vector[BestSolution][11];	wef[5] = vector[BestSolution][12];
		wef[6] = vector[BestSolution][13];	wef[7] = vector[BestSolution][14];

		wfe[1] = vector[BestSolution][15];	wfe[2] = vector[BestSolution][16];
		wfe[3] = vector[BestSolution][17];	wfe[5] = vector[BestSolution][18];
		wfe[6] = vector[BestSolution][19];	wfe[7] = vector[BestSolution][20];

		wff[1] = vector[BestSolution][21];	wff[2] = vector[BestSolution][22];
		wff[3] = vector[BestSolution][23];	wff[5] = vector[BestSolution][24];
		wff[6] = vector[BestSolution][25];	wff[7] = vector[BestSolution][26];
		




	

		
	for(k=0;k<iteration;k++)
	{
		result= Neural_Oscillator(ue, uf, ve, vf, ye, yf, wee, wef, wfe, wff, uinput, tc);

		fprintf(fout1,"%lf\n",result[0]);
		fprintf(fout2, "%lf\n", result[1]);
		fprintf(fout3, "%lf\n",  result[2]);

		//printf("%d: %lf  %lf  %lf \n", i, result[0], result[1], result[2]);
	}

	fclose(fout1);
	fclose(fout2);
	fclose(fout3);


	printf("Score = %lf\n", score[BestSolution]);
	fprintf(outfile,"Score = %lf", score[BestSolution]);
	fclose(outfile);
	fclose(outfile2);

	
}

void InitialVectors()
{
	int i, j;
	
	for(i=0; i< K ; i++)
		for(j=0; j< NumParameter; j++) 
			vector[i][j] = ((double)rand()/RAND_MAX) * 4. - 2.;
	count = K;
}






/////////////////////////
void ScoreAllVectors()
{
	int i, j, h, k, t1, t2;
	int flag = 0;
	double *result =0;
	double score_ev1, score_ev2, score_ev3; 

	for(i=0; i< K; i++){

	sum1 = 0, sum2 = 0, sum3 = 0;
	s1 = 0, s2 = 0, s3 = 0;



	
	uinput = vector[i][0];
	tc[0]  = 0.468;
	tc[1]  = 0.309; // value 3개 할당
	


	
		wee[1] = vector[i][3];	wee[2] = vector[i][4];
		wee[3] = vector[i][5];	wee[5] = vector[i][6];
		wee[6] = vector[i][7];	wee[7] = vector[i][8];

		wef[1] = vector[i][9];	wef[2] = vector[i][10];
		wef[3] = vector[i][11];	wef[5] = vector[i][12];
		wef[6] = vector[i][13];	wef[7] = vector[i][14];

		wfe[1] = vector[i][15];	wfe[2] = vector[i][16];
		wfe[3] = vector[i][17];	wfe[5] = vector[i][18];
		wfe[6] = vector[i][19];	wfe[7] = vector[i][20];

		wff[1] = vector[i][21];	wff[2] = vector[i][22];
		wff[3] = vector[i][23];	wff[5] = vector[i][24];
		wff[6] = vector[i][25];	wff[7] = vector[i][26];



		for(int j=0;j<iteration;j++)
		{	

			result= Neural_Oscillator(ue, uf, ve, vf, ye, yf, wee, wef, wfe, wff, uinput, tc);

			/*
			score_result[j+0] = result[0];
			score_result[j+1] = result[0];
			score_result[j+2] = result[0];	
			*/


			printf("%lf \n", result[1]);		

			
			
			
			if(j>500)
			{
				s1 = PI * j*0.005 / 2.0;
				s2 = PI * j*0.005 / 2.0 + (PI / 3.0);
				s3 = PI * j*0.005 / 2.0 + (2* PI / 3.0);

				sum1 += pow((result[0] - sin(s1)), 2);
				sum2 +=	pow((result[1] - sin(s2)), 2); 
				sum3 +=	pow((result[2] - sin(s3)), 2);
			}

			
		}

	
			score[i] = sum1 + sum2 + sum3;
		

	}

		
		
}



/////////////////////////
void GenerateOffspring() // offspring 생성. 다음번째 vector에 저장된다.
{
	int i, j;
	double Perturb;

	for(i=0; i< K ; i++)
        for(j=0; j< NumParameter; j++){
			Perturb = rndnor( 0., 1.224*1.224*score[i]/(NumParameter*NumParameter) );
			vector[i+K][j] = vector[i][j] + Perturb;
		}
	count+= K;
	K*= 2;
}

double rndnor(double mean, double var)
{

	//double mean, var;
	{
		double a, b, r;
		do{
			a = (double)rand() / RAND_MAX;
			b = (double)rand() / RAND_MAX;
		} while (a <= 0 || b <= 0);

		r = sqrt(var / 2.);
		r *= log(a / b);
		r += mean;
		return r;
	}
}

void CountWins()
{
	int opponent, i, j, Competitions;
	
	Competitions = 10;
	for(i=0; i< K ; i++){
		wins[i] = 0;
        for(j=0; j< Competitions; j++){
			do{
				opponent = rand() % K + 1;
			} while ( opponent == i );
			if( score[i] <= score[opponent] ) wins[i]++;
		}
	}
}

void RankPopulation()
{
	int i, j, best, index;
	double temp;

	for(i=0; i< K ; i++) compare[i] = wins[i];
	for(i=0; i< maximum; i++){
		best = i;
		for(j= i+1; j< K; j++)
			if( compare[j] > compare[best] ) best = j;
		temp = compare[i];
		compare[i] = compare[best];
		compare[best] = temp;
	}

	for(i=0; i< K ; i++){
		for(j=0; j< NumParameter; j++) save_vector[i][j] = 0;
		save_score[i] = 0;
	}

	for(i=0; i< maximum; i++){
		flag = False;
		index = 0;
		do{
			if( compare[i] == wins[index] ){
				flag = True;
				for(j=0; j< NumParameter; j++)
					save_vector[i][j] = vector[index][j];
				save_score[i] = score[index];
				wins[index] = -1;
			}
			index++;
		} while( index <= K && flag == False );
	}

	for(i=0; i< K ; i++){
		for(j=0; j< NumParameter; j++) vector[i][j] = 0;
		score[i] = 0;
	}

	for(i=0; i< maximum; i++){
		for(j=0; j< NumParameter; j++) vector[i][j] = save_vector[i][j];
		score[i] = save_score[i];
	}

	K = maximum;
}

void Statistics()
{
	double sum, sumsq, mean, var;
	int i;
	
	smallest = score[0];
	BestSolution = 0;
	sum = sumsq = 0;

	for(i=0; i< K ; i++){
		sum+= score[i];
		sumsq+= score[i] * score[i];
		if( score[i] < smallest ){
			smallest = score[i];
			BestSolution = i;
		}
	}

	mean = sum /K;
	if(K==1) var = 0;
	else var = (sumsq - sum*sum/K)/(K-1);
	
	printf("count=%d, mean=%lf, var=%lf, K=%d\n", count, mean, var, K);
	printf("smallest=%lf\n", smallest);
	fprintf(outfile,"count=%d, mean=%lf, var=%lf, K=%d\n", count, mean, var, K);
	fprintf(outfile2,"%lf\n", smallest);
}

