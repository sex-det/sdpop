#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <fenv.h>

// void gsl_ran_multinomial (const gsl_rng * r, size_t K, unsigned int N, const double p[], unsigned int n[])
// unsigned int gsl_ran_binomial (const gsl_rng * r, double p, unsigned int n)

double gsl_ran_beta_improved (const gsl_rng * r, const double a, const double b) //a bug in the standard implementation of gsl_ran_beta gives "nan" for small a and b
 {
  if ( (a <= 1.0) && (b <= 1.0) )
  {
    double U, V, X, Y;
    while (1)
    {
      U = gsl_rng_uniform_pos(r);
      V = gsl_rng_uniform_pos(r);
      X = pow(U, 1.0/a);
      Y = pow(V, 1.0/b);
      if ((X + Y ) <= 1.0)
      {
          if (X + Y > 0)
          {
              return X/ (X + Y);
          }
          else
          {
              double logX = log(U)/a;
              double logY = log(V)/b;
              double logM = logX > logY ? logX: logY;
              logX -= logM;
              logY -= logM;
              return exp(logX - log(exp(logX) + exp(logY)));
          }
      }
    }
  }
  else
  {
    double x1 = gsl_ran_gamma(r, a, 1.0);
    double x2 = gsl_ran_gamma(r, b, 1.0);
    return x1 / (x1 + x2);
  }
 }
 
double* errormult(double ematrix[3][3], double x[3]) {
	static double vector[3];
	long double temp;
	int i,j;
	
	for (i=0;i<3;i++) {
		temp=0.;
		for (j=0;j<3;j++) {
			temp+=(long double)x[j]*(long double)ematrix[i][j];
		}
		vector[i]=temp;
	}
	return vector;
}

void calcQ(double m33[3][3], double e0, double e1, double e2) {
	m33[0][0]=1.-e0-e2;
	m33[0][1]=e1;
	m33[0][2]=e2;
	m33[1][0]=e0;
	m33[1][1]=1.-2.*e1;
	m33[1][2]=e0;
	m33[2][0]=e2;
	m33[2][1]=e1;
	m33[2][2]=1.-e0-e2;
}


int main (int argc, char *argv[]) {
	//usage : simul n_individuals_per_sex nsites
	const gsl_rng_type *T;
	gsl_rng *r;
	unsigned int nmmulti[3],nfmulti[3],nerror[3],t,type,f_fixed=0;
	double *emulti,Q[3][3],pmulti[3],perror[3],f=0.5;
	double e0=0.01,e1=0.001,e2=0.001;
	double pi[4],rtype,theta;
//	int errorcount[6],totalerrors;
	int poly;
//	double errorrate=0.01;
	
	int i,k,ntemp;
	int n; //nombre d'individus
	int nsites; //nombre de sites polymorphes (par contig)
	int ncontigs; //nombre de contigs
	
//	feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);
	
	if (argc < 8) {
		fprintf(stderr,"Usage: %s n_individuals_per_sex nsites ncontigs e0 e1 e2 pi\n",argv[0]);
		exit(1);
	}
	
	n=atoi(argv[1]);
	nsites=atoi(argv[2]);
	ncontigs=atoi(argv[3]);
	e0=atof(argv[4]);
	e1=atof(argv[5]);
	e2=atof(argv[6]);
//	if (argc == 8) {
//		f = atof(argv[7]);
//		f_fixed=1;
//	}
	calcQ(Q,e0,e1,e2);
	theta=atof(argv[7]);
	
	for (i=0;i<3;i++){
		nfmulti[i]=nmmulti[i]=0;
	}
	
	gsl_rng_env_setup();
	
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	//todo : test if allocation was successfull
	
	//pi[0]-pi[2] : genome fractions
	pi[0]=0.33;
	pi[1]=0.33;
	pi[2]=1.-pi[0]-pi[1];
	//pi[3] : only concerns XY polymorphisms
	pi[3]=0.5;
	
	for (k=0;k<ncontigs;k++) {
		//choisir le type de ségrégation
		rtype=(double)rand()/RAND_MAX;
		if (rtype <= pi[0]) {
			type=0;
		}
		else if (rtype <= pi[0]+pi[1]) {
			type=1;
		}
		else {
//		else if (rtype <= pi[0]+pi[1]+pi[2]) {
			type=2;
		}
//		else {
//			type=3;
//		}
		printf(">contig%d\t%d\n",k,type);
		for (t=0;t<nsites;t++) {  
			//si on ne simule qu'un contig, choisir le type de ségrégation par site
			if(ncontigs==1) {
				rtype=(double)rand()/RAND_MAX;
				if (rtype <= pi[0]) {
					type=0;
				}
				else if (rtype <= pi[0]+pi[1]) {
					type=1;
				}
				else if (rtype <= pi[0]+pi[1]+pi[2]*pi[3]) {
					type=2;
				}
				else {
					type=3;
				}
			}
			else if (type==2 || type==3) {
				//à l'intérieur du type XY, il faut choisir si le polymorphisme est sur X ou sur Y
				if ((double)rand()/RAND_MAX <= pi[3]) {
					type=2;
				}
				else {
					type=3;
				}
			}
			
			if(!f_fixed) {
//				f=(double)rand()/RAND_MAX;
//				f=1.;
//				while(f>0.5){
				f=gsl_ran_beta_improved(r,theta,theta);
//				}
			}
			
			if (type == 0) { //autosomal ; f is the frequency of allele 1
				pmulti[0]=f*f;
				pmulti[1]=2.*f*(1.-f);
				pmulti[2]=(1.-f)*(1.-f);
				emulti=errormult(Q, pmulti);
				gsl_ran_multinomial(r, 3, (unsigned)n, emulti, nfmulti);
				gsl_ran_multinomial(r, 3, (unsigned)n, emulti, nmmulti);
				//check for polymorphism
				if ( (nmmulti[0]==n && nfmulti[0]==n) || (nmmulti[2]==n && nfmulti[2]==n)) {
//				if ( (nfmulti[0]==n) || (nfmulti[2]==n)) {
//					t--;
					continue;
				}
				//randomize: invert alleles 1 and 2 (heterozygote counts don't change)
				if ((double)rand()/RAND_MAX <=0.5 ) {
					ntemp=nfmulti[0];
					nfmulti[0]=nfmulti[2];
					nfmulti[2]=ntemp;
					ntemp=nmmulti[0];
					nmmulti[0]=nmmulti[2];
					nmmulti[2]=ntemp;
				}
			}
			else if (type == 1) { //X-hémi ; f is the frequency of allele 1 on X
				//femelles
				pmulti[0]=f*f;
				pmulti[1]=2.*f*(1.-f);
				pmulti[2]=(1.-f)*(1.-f);
				emulti=errormult(Q, pmulti);
				gsl_ran_multinomial(r, 3, (unsigned)n, emulti, nfmulti);
				//mâles
				pmulti[0]=f;
				pmulti[1]=0.;
				pmulti[2]=1.-f;
				emulti=errormult(Q, pmulti);
				gsl_ran_multinomial(r, 3, (unsigned)n, emulti, nmmulti);
				
//				nmmulti[0]=gsl_ran_binomial(r,f,n);
//				nmmulti[1]=0;
//				nmmulti[2]=n-nmmulti[0];

				//check for polymorphism
				if ( (nmmulti[0]==n && nfmulti[0]==n) || (nmmulti[2]==n && nfmulti[2]==n)) {
//					t--;
					continue;
				}
				//randomize: invert alleles 1 and 2 (heterozygote counts don't change)
				if ((double)rand()/RAND_MAX <=0.5 ) {
					ntemp=nfmulti[0];
					nfmulti[0]=nfmulti[2];
					nfmulti[2]=ntemp;
					ntemp=nmmulti[0];
					nmmulti[0]=nmmulti[2];
					nmmulti[2]=ntemp;
				}
			}
			else if (type == 2) { //XY ; x-polymorphism
				// deux cas
				//simulate case 1: allele 1 is fixed on Y. f is the frequency of allele 2 on X
				//femelles
				pmulti[0]=(1.-f)*(1.-f);
				pmulti[1]=2.*f*(1.-f);
				pmulti[2]=f*f;
				emulti=errormult(Q, pmulti);
				gsl_ran_multinomial(r, 3, (unsigned)n, emulti, nfmulti);
//				gsl_ran_multinomial(r, 3, (unsigned)n, pmulti, nfmulti);
				//mâles
//				nmmulti[1]=gsl_ran_binomial(r,f,n);
//				nmmulti[0]=n-nmmulti[1];
//				nmmulti[2]=0;
				pmulti[0]=1.-f;
				pmulti[1]=f;
				pmulti[2]=0;				
				emulti=errormult(Q, pmulti);
				gsl_ran_multinomial(r, 3, (unsigned)n, emulti, nmmulti);

				//check for polymorphism
				if ( (nmmulti[0]==n && nfmulti[0]==n) || (nmmulti[2]==n && nfmulti[2]==n)) {
//					t--;
					continue;
				}
				//case 2: invert alleles 1 and 2 (heterozygote counts don't change)
				if ((double)rand()/RAND_MAX <=0.5 ) {
					ntemp=nfmulti[0];
					nfmulti[0]=nfmulti[2];
					nfmulti[2]=ntemp;
					ntemp=nmmulti[0];
					nmmulti[0]=nmmulti[2];
					nmmulti[2]=ntemp;
				}
			}
			else if (type == 3) { //XY ; y-polymorphism
				// deux cas
				//simulate case 1: allele 1 is fixed on X; f is the frequency of allele 2 on Y
				//femelles
//				nfmulti[0]=n; //N11f
//				nfmulti[1]=0; //N12f
//				nfmulti[2]=0; //N22f
				pmulti[0]=1.;
				pmulti[1]=0.;
				pmulti[2]=0.;
				emulti=errormult(Q, pmulti);
				gsl_ran_multinomial(r, 3, (unsigned)n, emulti, nfmulti);
				//mâles
//				nmmulti[1]=gsl_ran_binomial(r,f,n); //N12m
//				nmmulti[0]=n-nmmulti[1]; //N11m
//				nmmulti[2]=0; //N22m
				pmulti[0]=1.-f;
				pmulti[1]=f;
				pmulti[2]=0;				
				emulti=errormult(Q, pmulti);
				gsl_ran_multinomial(r, 3, (unsigned)n, emulti, nmmulti);
				//check for polymorphism
				if ( (nmmulti[0]==n && nfmulti[0]==n) || (nmmulti[2]==n && nfmulti[2]==n)) {
//					t--;
					continue;
				}
				//case 2: invert alleles 1 and 2 (heterozygote counts don't change)
				if ((double)rand()/RAND_MAX <=0.5 ) {
					ntemp=nfmulti[0];
					nfmulti[0]=nfmulti[2];
					nfmulti[2]=ntemp;
					ntemp=nmmulti[0];
					nmmulti[0]=nmmulti[2];
					nmmulti[2]=ntemp;
				}
			}
			printf("%d\t AT \t",t);
			for (i=0;i<3;i++) {
				printf("%d\t",nfmulti[i]);
			}
			for (i=0;i<3;i++) {
				printf("%d\t",nmmulti[i]);
			}
			printf("%d\t%f\n",type,f);
		}
	}
	gsl_rng_free (r);
	
	return 0;
}

