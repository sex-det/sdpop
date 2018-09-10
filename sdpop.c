#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
//#include <gsl/gsl_randist.h>
#include <gsl/gsl_machine.h>
#include <fenv.h>

#define FEMALE 0
#define MALE 1
#define SEXES 2
#define NTYPES 3 //number of biological segregation types
#define NRTYPES 4 //number of XY segregation types
#define NATYPES 6 //number of detailed segregation types
#define SITE 0
#define CONTIG 1
#define NONE 0
#define ERRORS 1
#define FIXED 2
#define N11 1
#define N12 2
#define N22 3
#define N11F 1
#define N12F 2
#define N22F 3
#define N11M 4
#define N12M 5
#define N22M 6
#define XY 1
#define ZW 2

double ***F; //vector containing allele frequencies per segregation type
double *****P; //vector containing P matrices (per sex and segregation type)
long double ***condsiteprob,***condsegprob; //conditional probabilities per site
double ***expS,***expA,**expR,******expTG;

double intpow(double x,int y){
	double result;
	int i;
	if(y==0){
		result=1.;
	}
	else if (y<0){
		result=1./x;
		for (i=-1;i>y;i--){
			result*=(1./x);
		}
	}
	else {
		result=x;
		for (i=1;i<y;i++){
			result*=x;
		}
	}
	return result;
}

void initEM(int ncontigs, int *npolysites, int ***polysite) {
	int n11f,n12f,n22f,n11m,n12m,n22m,nfem,nmal,ntot;
	int k,t,j,s,g;
	double f;
	
	//P : a npolysites x 2 x 4 x 3 matrix
	//P[kt][s][j][gp] : 
	//kt: sites
	//s: sexes, 0=female, 1=male
	//j: segregation types
	//gp: true genotypes, 0="11", 1="12", 2="22"
	
	if((P=(double *****)calloc((size_t)ncontigs,sizeof(double ****)))==NULL) { 
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}
	if((F=(double ***)calloc((size_t)ncontigs,sizeof(double **)))==NULL) { 
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}	
	if((expS=(double ***)calloc((size_t)ncontigs,sizeof(double **)))==NULL) { 
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}
	if((expA=(double ***)calloc((size_t)ncontigs,sizeof(double **)))==NULL) { 
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}
	if((condsiteprob=(long double ***)calloc((size_t)ncontigs,sizeof(long double **)))==NULL) { 
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}
	if((condsegprob=(long double ***)calloc((size_t)ncontigs,sizeof(long double **)))==NULL) { 
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}
	if((expR=(double **)calloc((size_t)ncontigs,sizeof(double *)))==NULL) { 
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}
	if((expTG=(double ******)calloc((size_t)ncontigs,sizeof(double *****)))==NULL) { 
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}
	
	
	for (k=0;k<ncontigs;k++) {
		if (npolysites[k] > 0) {
			if((P[k]=(double ****)calloc((size_t)npolysites[k],sizeof(double ***)))==NULL) { 
				fprintf(stderr,"error in memory allocation\n");
				exit(1);
			}
			if((F[k]=(double **)calloc((size_t)npolysites[k],sizeof(double *)))==NULL) { 
				fprintf(stderr,"error in memory allocation\n");
				exit(1);
			}
			if((expS[k]=(double **)calloc((size_t)npolysites[k],sizeof(double *)))==NULL) { 
				fprintf(stderr,"error in memory allocation\n");
				exit(1);
			}
			if((expA[k]=(double **)calloc((size_t)npolysites[k],sizeof(double *)))==NULL) { 
				fprintf(stderr,"error in memory allocation\n");
				exit(1);
			}
			if((condsiteprob[k]=(long double **)calloc((size_t)npolysites[k],sizeof(long double *)))==NULL) { 
				fprintf(stderr,"error in memory allocation\n");
				exit(1);
			}
			if((condsegprob[k]=(long double **)calloc((size_t)npolysites[k],sizeof(long double *)))==NULL) { 
				fprintf(stderr,"error in memory allocation\n");
				exit(1);
			}
			if((expR[k]=(double *)calloc((size_t)NTYPES,sizeof(double)))==NULL) { 
				fprintf(stderr,"error in memory allocation\n");
				exit(1);
			}
			if((expTG[k]=(double *****)calloc((size_t)npolysites[k],sizeof(double ****)))==NULL) { 
				fprintf(stderr,"error in memory allocation\n");
				exit(1);
			}
			for (t=0; t<npolysites[k]; t++){
				if((P[k][t]= (double ***)calloc((size_t)SEXES,sizeof(double **)))==NULL) { 
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}
				if((F[k][t]= (double *)calloc((size_t)NTYPES,sizeof(double)))==NULL) { 
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}
				if((expS[k][t]= (double *)calloc((size_t)NTYPES,sizeof(double)))==NULL) { 
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}
				if((expA[k][t]= (double *)calloc((size_t)NRTYPES,sizeof(double)))==NULL) { 
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}
				if((condsiteprob[k][t]= (long double *)calloc((size_t)NATYPES,sizeof(long double)))==NULL) { 
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}
				if((condsegprob[k][t]= (long double *)calloc((size_t)NTYPES,sizeof(long double)))==NULL) { 
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}
				if((expTG[k][t]= (double ****)calloc((size_t)SEXES,sizeof(double ***)))==NULL) { 
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}
				for(s=0;s<SEXES;s++){
					if((P[k][t][s]= (double **)calloc((size_t)NATYPES,sizeof(double *)))==NULL) { 
						fprintf(stderr,"error in memory allocation\n");
						exit(1);
					}
					if((expTG[k][t][s]= (double ***)calloc((size_t)NATYPES,sizeof(double **)))==NULL) { 
						fprintf(stderr,"error in memory allocation\n");
						exit(1);
					}
					for(j=0;j<NATYPES;j++){
						if((P[k][t][s][j]= (double *)calloc((size_t)3,sizeof(double)))==NULL) { 
							fprintf(stderr,"error in memory allocation\n");
							exit(1);
						}
						if((expTG[k][t][s][j]= (double **)calloc((size_t)3,sizeof(double *)))==NULL) { 
							fprintf(stderr,"error in memory allocation\n");
							exit(1);
						}
						for(g=0;g<3;g++){
							if((expTG[k][t][s][j][g]= (double *)calloc((size_t)3,sizeof(double)))==NULL) { 
								fprintf(stderr,"error in memory allocation\n");
								exit(1);
							}
						}
					}
				}
			}							
			
			for (t=0;t<npolysites[k];t++) {
				n11f=polysite[k][t][N11F]; //female counts
				n12f=polysite[k][t][N12F]; 
				n22f=polysite[k][t][N22F]; 			
				n11m=polysite[k][t][N11M]; //male counts
				n12m=polysite[k][t][N12M]; 
				n22m=polysite[k][t][N22M]; 			
				
				nfem=n11f+n12f+n22f;
				nmal=n11m+n12m+n22m;
				ntot=nfem+nmal;
				
				//autosomal snps. f0 is the frequency of allele 1 (symmetry)
				j=0;
				f=(double)(2*n11f+2*n11m+n12f+n12m)/(double)(2*ntot);
				F[k][t][0]=f;
				P[k][t][FEMALE][j][0]=f*f;
				P[k][t][FEMALE][j][1]=2.*f*(1.-f);
				P[k][t][FEMALE][j][2]=(1.-f)*(1.-f);
				P[k][t][MALE][j][0]=P[k][t][FEMALE][j][0];
				P[k][t][MALE][j][1]=P[k][t][FEMALE][j][1];
				P[k][t][MALE][j][2]=P[k][t][FEMALE][j][2];
				
				//x-hemizygous snps. f1 is the frequency of allele 1 (symmetry)
				j=1;
				f=(double)(2*n11f+n12f+n11m)/(double)(2*nfem+nmal);
				F[k][t][1]=f;
				P[k][t][FEMALE][j][0]=f*f;
				P[k][t][FEMALE][j][1]=2.*f*(1.-f);
				P[k][t][FEMALE][j][2]=(1.-f)*(1.-f);
				P[k][t][MALE][j][0]=f;
				P[k][t][MALE][j][1]=0.;
				P[k][t][MALE][j][2]=1-f;
				
				//x/y snsp ; x-polymorphism
				//allele 1 is fixed on Y; f21 is the frequency of allele 2 on X
				j=2;
				f=(double)(2*n22f+n12f+n12m)/(double)(2*nfem+n11m+n12m);
				P[k][t][FEMALE][j][0]=(1.-f)*(1.-f);
				P[k][t][FEMALE][j][1]=2.*f*(1.-f);
				P[k][t][FEMALE][j][2]=f*f;
				P[k][t][MALE][j][0]=1.-f;
				P[k][t][MALE][j][1]=f;
				P[k][t][MALE][j][2]=0.;
				//allele 2 is fixed on Y; f22 is the frequency of allele 1 on X
				j=3;
				f=(double)(2*n11f+n12f+n12m)/(double)(2*nfem+n12m+n22m);
				P[k][t][FEMALE][j][0]=f*f;
				P[k][t][FEMALE][j][1]=2.*f*(1.-f);
				P[k][t][FEMALE][j][2]=(1.-f)*(1.-f);
				P[k][t][MALE][j][0]=0.;
				P[k][t][MALE][j][1]=f;
				P[k][t][MALE][j][2]=1.-f;
				//x/y snsp ; y-polymorphism. f3 is the frequency, on Y, of the allele that is not fixed on X
				//case 1: allele 1 is fixed on X; f3 is the frequency of allele 2 on Y
				j=4;
				//f=(double)(n12m)/(double)(n11m+n12m);
				f=(double)(n12m)/(double)(nmal);
				P[k][t][FEMALE][j][0]=1.;
				P[k][t][FEMALE][j][1]=0.;
				P[k][t][FEMALE][j][2]=0.;
				P[k][t][MALE][j][0]=1.-f;
				P[k][t][MALE][j][1]=f;
				P[k][t][MALE][j][2]=0.;
				//case 2: allele 2 is fixed on X; f3 is the frequency of allele 1 on Y
				j=5;
				//f=(double)(n12m)/(double)(n22m+n12m);
				f=(double)(n12m)/(double)(nmal);
				P[k][t][FEMALE][j][0]=0.;
				P[k][t][FEMALE][j][1]=0.;
				P[k][t][FEMALE][j][2]=1.;
				P[k][t][MALE][j][0]=0.;
				P[k][t][MALE][j][1]=f;
				P[k][t][MALE][j][2]=1.-f;
			}		
		}
	}
}

void freeEM(int ncontigs, int *npolysites) {
	int k,t,s,j,g;
	
	for(k=0; k<ncontigs; k++){
		for (t=0; t<npolysites[k]; t++){
			for(s=0;s<SEXES;s++){
				for(j=0;j<NTYPES;j++){
					free(P[k][t][s][j]);
				}
				for(g=0;g<3;g++){
					free(expTG[k][t][s][g]);
				}
				free(P[k][t][s]);			
				free(expTG[k][t][s]);			
			}
			free(P[k][t]);
			free(expTG[k][t]);
			free(F[k][t]);
			free(expS[k][t]);
			free(expA[k][t]);
			free(condsiteprob[k][t]);
			free(condsegprob[k][t]);
		}
		free(expTG[k]);
		free(P[k]);
		free(F[k]);
		free(expS[k]);
		free(expA[k]);
		free(condsiteprob[k]);
		free(condsegprob[k]);
		free(expR[k]);
	}
	free(P);
	free(expTG);
	free(F);
	free(expS);
	free(expA);
	free(condsiteprob);
	free(condsegprob);
	free(expR);
}

void pisort(int n, long double *pi, int *piorder) {
	//create an array piorder with an order of j's such that pi[piorder[j]] are in descending order
	int i,ii,j,jmax,jtest;
	double pmax;
	for (i=0; i<n; i++) {
		jmax=-1;
		pmax=-1.0;
		for (j=0;j<n;j++) {
			jtest=1;
			for (ii=0;ii<i;ii++) {
				if(j==piorder[ii]) {
					jtest=0;
				}
			}
			if (jtest && pi[j]>pmax) {
				jmax=j;
				pmax=pi[j];
			}
		}
		piorder[i]=jmax;
	}
}

double* vectorize_d(double x0, double x1, double x2) {
	static double vector[3];
	vector[0]=x0;
	vector[1]=x1;
	vector[2]=x2;
	return vector;
}

unsigned int* vectorize_ui(unsigned int x0, unsigned int x1, unsigned int x2) {
	static unsigned int vector[3];
	vector[0]=x0;
	vector[1]=x1;
	vector[2]=x2;
	return vector;
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

void CondSiteProbs(int ncontigs, int *npolysites, int ***polysite, double Q[3][3])
{
	int k,t,jl,s,g,gp;
//	unsigned int *nmulti;
//	double *pmulti, *emulti;
	double tempgp;
	
	for (k=0;k<ncontigs;k++) {
		for(t=0;t<npolysites[k];t++){
			for(jl=0;jl<NATYPES;jl++) {
//				nmulti=vectorize_ui(polysite[k][t][N11F],polysite[k][t][N12F],polysite[k][t][N22F]);
//				pmulti=vectorize_d(P[k][t][FEMALE][jl][0],P[k][t][FEMALE][jl][1],P[k][t][FEMALE][jl][2]);
//				emulti=errormult(Q, pmulti);
//				condsiteprob[k][t][jl]=gsl_ran_multinomial_pdf(3,emulti,nmulti);
//				nmulti=vectorize_ui(polysite[k][t][N11M],polysite[k][t][N12M],polysite[k][t][N22M]);
//				pmulti=vectorize_d(P[k][t][MALE][jl][0],P[k][t][MALE][jl][1],P[k][t][MALE][jl][2]);
//				emulti=errormult(Q, pmulti);
//				condsiteprob[k][t][jl]*=gsl_ran_multinomial_pdf(3,emulti,nmulti);
//				if(isnan(condsiteprob[k][t][jl])) { //As gsl_ran_multinomial_pdf uses gsl_ran_multinomial_lnpdf (logarithmic version), 0 probabilities result in NaN
//					printf("Warning: NaN produced: contig %d, site %d, type %d \n",k,t,jl);
//					condsiteprob[k][t][jl]=0.;
//				}
				condsiteprob[k][t][jl]=1;
				for(s=0;s<2;s++){
					for(g=0;g<3;g++){
						tempgp=0.;
						for(gp=0;gp<3;gp++){
							tempgp+=P[k][t][s][jl][gp]*Q[g][gp];
						}
						condsiteprob[k][t][jl]*=intpow(tempgp,polysite[k][t][1+3*s+g]);
					}
				}
			}
		}
	}
}

void CondSegProbs(int ncontigs, int *npolysites, double *rho)
{
	int k,t,j,l;
	
	for (k=0;k<ncontigs;k++) {
		for(t=0;t<npolysites[k];t++){
			for(j=0;j<NTYPES;j++){
				if(j<2){
					condsegprob[k][t][j]=condsiteprob[k][t][j];
				}
				else {
					condsegprob[k][t][j]=0;
					for(l=0;l<NRTYPES;l++) {
						condsegprob[k][t][j]+=(long double)rho[l]*condsiteprob[k][t][j+l];
					}
				}
			}
		}
	}
}

long double totalcontigloglik(int ncontigs, int *npolysites, double *pi, double *rho) {
	long double loglik=0;
	long double suml,sumj,temp[NTYPES];
	int t,j,k,l;
	
	for (k=0;k<ncontigs;k++) {
		sumj=0;
		for(j=0;j<NTYPES;j++) {
			if(j<2){
				temp[j]=log(pi[j]);
				for(t=0;t<npolysites[k];t++){
					temp[j]+=logl(condsiteprob[k][t][j]);
				}
			}
			else {
				temp[j]=log(pi[2]); //XY probabilities
				for(t=0;t<npolysites[k];t++){
					suml=0;
					for(l=0;l<NRTYPES;l++){
						suml+=(long double)rho[l]*condsiteprob[k][t][j+l];
					}
					temp[j]+=logl(suml);
				}
			}
		}
		for(j=0;j<NTYPES;j++) {
			sumj+=expl(temp[j]);
		}
		loglik+=logl(sumj);
	}
	return loglik;
}

long double totalsiteloglik(int ncontigs, int *npolysites, double *pi, double *rho) {
	long double loglik=0;
	long double sumj,sum3;
	int t,j,k,l;
	
	for (k=0;k<ncontigs;k++) {
		for(t=0;t<npolysites[k];t++){
			sumj=0;
			for(j=0;j<NTYPES;j++) {
				if(j<2){
					sumj+=condsiteprob[k][t][j]*(long double)pi[j];
				}
				else {
					sum3=0;
					for(l=0;l<NRTYPES;l++){
						sum3+=condsiteprob[k][t][j+l]*(long double)rho[l];
					}
					sumj+=sum3*(long double)pi[j];
				}
			}
			loglik+=logl(sumj);
		}
	}
	return loglik;
}

long double condsiteexpect(int ncontigs, int *npolysites, int ***polysite, double *pi, double *rho, double Q[3][3])
{
	long double tempgs,tempjl,tempgp,logexp=0;
	int k,t,s,jl,j,l,g,gp;
	
	for (k=0;k<ncontigs;k++) {
		for(t=0;t<npolysites[k];t++){
			tempgs=0;
			for(s=0;s<2;s++){
				for(g=0;g<3;g++){
					tempjl=0;
					for(jl=0;jl<NATYPES;jl++){
						tempgp=0;
						for(gp=0;gp<3;gp++){
							if(P[k][t][s][jl][gp]>GSL_DBL_MIN){
								tempgp+=expTG[k][t][s][jl][g][gp]*(log(Q[g][gp])+log(P[k][t][s][jl][gp]));
							}
						}
						if(jl<2){
							tempjl+=tempgp*expS[k][t][jl];
						}
						else {
							tempjl+=tempgp*expS[k][t][2]*expA[k][t][jl-2];
						}
					}					
				tempgs+=polysite[k][t][1+g+3*s]*tempjl;
				}
			}
			logexp+=tempgs;
			for(l=0;l<NRTYPES;l++){
				logexp+=expS[k][t][2]*expA[k][t][l]*log(rho[l]);
			}
			for(j=0;j<NTYPES;j++){
				logexp+=expS[k][t][j]*log(pi[j]);
			}
		}
	}
	return logexp;
}

void horner(int length, int lmax, long double* array, double* param, double* result){
	long double tot;
	int i;
	tot=param[lmax];
	for(i=0;i<length;i++){
		if(i!=lmax){
			tot+=param[i]*array[i]/array[lmax];
		}
	}
	tot*=array[lmax];
	for(i=0;i<length;i++){
		result[i]=param[i]*array[i]/tot;
	}
}

void loghorner(int length, int lmax, long double* logarray, double* param, double* result){
	long double tot,ltot;
	int i;
//	long double temp[3];
	tot=(long double)param[lmax];
	for(i=0;i<length;i++){
		if(i!=lmax){
			tot+=(long double)param[i]*expl(logarray[i]-logarray[lmax]);
		}
	}
	ltot=logl(tot)+logarray[lmax];
	for(i=0;i<length;i++){
//		temp[0]=logarray[i];
//		temp[1]=logl((long double)param[i]);
//		temp[2]=ltot;
		if(log(param[i])>1){
			fprintf(stdout,"error: strange log behaviour: %f\t%f\n",param[i],log(param[i]));
		}
		result[i]=(double)expl(logarray[i]+(long double)log(param[i])-ltot);
	}
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


int main(int argc, char *argv[]) {
	
	FILE *fp,*outfile;
	int CUR_MAX = 4095;
	int NAME_LEN = 100;
	int NCONTIG_BATCH = 100;
	int ncontigs_allocated;
	char *line = calloc((size_t)CUR_MAX,sizeof(char));
	char *tmpline = calloc((size_t)CUR_MAX,sizeof(char));
	char *contig;
	int count = 0; 
	int length = 0;
	char ch;
	int i,j,k,l,jl,firstcontig,t,it,s,g,gp;
	int ***polysite;
	int *npolysites;
	int ncontigs,totsites;
	double sumET,weights;     
	double pi[NTYPES],rho[NRTYPES],oldpi[NTYPES],oldrho[NRTYPES],sumpi,sumrho;
	long double pidelta[NTYPES],rhodelta[NRTYPES],edelta[3],deltamax;
	int mode=CONTIG,errormodel=ERRORS,plateausteps,nnoncontigs,chromosomes=XY;
	double stop=0.01; //relative difference between parameter values used to signal convergence
	double minimumvalue=1e-10; //if a parameter value falls below this value, stop evaluating its optimisation
	double Q[3][3],e[3],olde[3];
	double maxj,maxl,sumgp;
	long double temp[NTYPES],templ[NRTYPES],maxgp;
//	long double logexp,sumj,lsumj;
	double pmax;
	int	jmax,lmax,gpmax;
	long double oldloglik,loglik;
	double U[3][3];
	int ni,nI=0;
	
//	feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);
	
	for(i=0;i<argc;i++) {
		fprintf(stdout,"%s ",argv[i]);
	}
	fprintf(stdout,"\n");
	
	if (argc != 6) {
		fprintf(stdout,"Usage: %s infile outfile mode errormodel heterogamety\n",argv[0]);
		exit(1);
	}
	
	if((fp=fopen(argv[1],"r"))==NULL){
		fprintf(stderr,"error opening input file %s\n",argv[1]);
		exit(1);
	}
	if((outfile=fopen(argv[2],"w"))==NULL){
		fprintf(stderr,"error opening output file %s\n",argv[2]);
		exit(1);
	}
	if(strcmp(argv[3],"c")==0 || strcmp(argv[3],"1")==0) {
		mode=CONTIG;
	}
	else if (strcmp(argv[3],"s")==0 || strcmp(argv[3],"0")==0) {
		mode=SITE;
	}
	else {
		fprintf(stderr,"Usage: %s infile outfile mode errormodel heterogamety\n",argv[0]);
		fprintf(stderr,"Mode should be either \"c\" or \"1\" for contig-mode, or \"s\" or \"0\" for site-wise optimisation\n");
		exit(1);	
	}
	if(strcmp(argv[4],"n")==0 || strcmp(argv[4],"0")==0) {
		errormodel=NONE;
	}
	else if (strcmp(argv[4],"e")==0 || strcmp(argv[4],"1")==0) {
		errormodel=ERRORS;
	}
	else if (strcmp(argv[4],"f")==0 || strcmp(argv[4],"2")==0) {
		errormodel=FIXED;
	}
	else {
		fprintf(stderr,"Usage: %s infile outfile mode errormodel heterogamety\n",argv[0]);
		fprintf(stderr,"Errormodel should be either \"e\" or \"1\" to estimate errors, \"f\" or \"2\" to use fixed error rates, or \"n\" or \"0\" for no errors\n");
		exit(1);	
	}
	if(strcmp(argv[5],"n")==0 || strcmp(argv[5],"0")==0) {
		chromosomes=NONE;
	}
	else if (strcmp(argv[5],"x")==0 || strcmp(argv[5],"1")==0) {
		chromosomes=XY;
	}
	else if (strcmp(argv[5],"z")==0 || strcmp(argv[5],"2")==0) {
		chromosomes=ZW;
	}
	else {
		fprintf(stderr,"Usage: %s infile outfile mode errormodel heterogamety\n",argv[0]);
		fprintf(stderr,"Heterogamety should be either be \"x\" or \"1\" for XY type, \"z\" or \"2\" for ZW type, or \"n\" or \"0\" for no sex chromosomes\n");
		exit(1);	
	}
		
	//reading counts from file 
	l=0; //line number
	ch='a';
	if((npolysites=(int *)malloc(sizeof(int)*NCONTIG_BATCH))==NULL) {
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}
	if((polysite=(int ***)malloc(sizeof(int **)*NCONTIG_BATCH))==NULL) {
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}
	if((contig=(char *)malloc(sizeof(char)*NCONTIG_BATCH*NAME_LEN))==NULL) { 
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}
	ncontigs_allocated=NCONTIG_BATCH;
	
	firstcontig=1;
	k=0;
	
	fprintf(stdout,"Reading data...\n");
	
	while (ch != EOF) { //loop through the file
		ch='a';
		count = 0;
		length = 0;
		while ( (ch != '\n') && (ch != EOF) ) { //loop through the line
			if(count == CUR_MAX) { // time to expand (for unexepectedly large line lengths) ?
				CUR_MAX *= 2; 
				count = 0;
				line = realloc(line, sizeof(char) * CUR_MAX); 
				tmpline = realloc(tmpline, sizeof(char) * CUR_MAX); 
			}
			ch = getc(fp); // read from stream.
			line[length] = ch;
			length++;
			count++;
		}
		line[length] = '\0';
		if (length <= 1) { //empty line : suppose it's the end of the file
			break;
		}
		//We've read one line :
		l++;
		
		if ( line[0] == '>' ){
			if (firstcontig==0) {	//Filtering polymorphisms for the last read contig
				npolysites[k]=t;
				t=0;
				k++;
			}
			
			if (firstcontig==1) {
				firstcontig=0;
				t=0;
			}
			//Start preparing to read a new contig
			if ( k >= ncontigs_allocated ) {
				ncontigs_allocated+=NCONTIG_BATCH;
				if((contig=(char *)realloc(contig,sizeof(char)*ncontigs_allocated*NAME_LEN))==NULL) { 
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}	
				if((polysite=(int ***)realloc(polysite,sizeof(int **)*ncontigs_allocated))==NULL) {
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}		
				if((npolysites=(int *)realloc(npolysites,sizeof(int **)*ncontigs_allocated))==NULL) {
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}		
			}
			if ( strlen(line) >= NAME_LEN ) {
				fprintf(stderr,"Error: a contig name is too long (more than %d characters) on line %d.\n",NAME_LEN-1,l);
				exit(1);
			}
			sscanf(line,">%s",&contig[k*NAME_LEN]);
			//			printf("%s\n",&contig[k*NAME_LEN]);
		}
		else { //line contains counts
				if (t==0) {
					if((polysite[k]=(int **)malloc(sizeof(int *)))==NULL){
						fprintf(stderr,"error in memory allocation\n");
						exit(1);
					}
					if((polysite[k][t]=(int *)malloc(sizeof(int)*7))==NULL){
						fprintf(stderr,"error in memory allocation\n");
						exit(1);
					}
				}
				else {
					if((polysite[k]=realloc(polysite[k],sizeof(int *)*(t+1)))==NULL){
						fprintf(stderr,"error in memory allocation\n");
						exit(1);
					}
					if((polysite[k][t]=(int *)malloc(sizeof(int)*7))==NULL){
						fprintf(stderr,"error in memory allocation\n");
						exit(1);
					}
				}
//				polysite[k][t][0]=t;
//				sscanf(line,"%*d\t%*f\t%d\t%d\t%d\t%d\t%d\t%d\t",&polysite[k][t][1],&polysite[k][t][2],&polysite[k][t][3],&polysite[k][t][4],&polysite[k][t][5],&polysite[k][t][6]);
//				if(chromosomes==XY){
//					sscanf(line,"%*d\t%*d\t%*f\t%d\t%d\t%d\t%d\t%d\t%d\t",&polysite[k][t][1],&polysite[k][t][2],&polysite[k][t][3],&polysite[k][t][4],&polysite[k][t][5],&polysite[k][t][6]);
//				}
//				else {
//					sscanf(line,"%*d\t%*d\t%*f\t%d\t%d\t%d\t%d\t%d\t%d\t",&polysite[k][t][4],&polysite[k][t][5],&polysite[k][t][6],&polysite[k][t][1],&polysite[k][t][2],&polysite[k][t][3]);
//				}
				if(chromosomes==XY){
					sscanf(line,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t",&polysite[k][t][0],&polysite[k][t][1],&polysite[k][t][2],&polysite[k][t][3],&polysite[k][t][4],&polysite[k][t][5],&polysite[k][t][6]);
				}
				else {
					sscanf(line,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t",&polysite[k][t][0],&polysite[k][t][4],&polysite[k][t][5],&polysite[k][t][6],&polysite[k][t][1],&polysite[k][t][2],&polysite[k][t][3]);
				}
				ni=0;
				for(i=1;i<=6;i++){
					ni+=polysite[k][t][i];
				}
				if(ni>nI){
					nI=ni;
				}
				t++;
		}
		memset(line, '\0', strlen(line)*(sizeof line));
	}
//	fprintf(stdout,"%d sites (individuals * positions), error rate (reads2snp) = %f\n",gensite,1.-totgenprob/gensite);
	
	npolysites[k]=t;
	ncontigs=k+1;
		
	//All reading has been done : clear memory.
	free(tmpline);
	free(line);
	fclose(fp);

	//end of data reading and filtering
	
	fprintf(stdout,"Found %d contigs\n",ncontigs);
	totsites=0;
	nnoncontigs=0;
	for (k=0;k<ncontigs;k++) {
		totsites+=npolysites[k];
		if(npolysites[k]==0) {
			nnoncontigs++;
		}
	}
	fprintf(stdout,"...and %d polymorphic sites\n",totsites);
	if(totsites==0){
		fprintf(stdout,"No polymorphic sites found; nothing to do.\n");
		exit(0);
	}
	
	//EM algorithm
	
	// Allocation
	initEM(ncontigs, npolysites, polysite);
	
	// Initialize parameters
	if(chromosomes==NONE){
		pi[0]=1.;
		pi[1]=0.;
		pi[2]=0.;
	}
	else {
	for(j=0;j<NTYPES;j++) {
		pi[j]=1./NTYPES;
	}
//	pi[0]=0.03;
//	pi[1]=0.9;
//	pi[2]=0.07;
	}
	for(l=0;l<NRTYPES;l++) {
		rho[l]=1./NRTYPES;
	}
//	rho[0]=0.03;
//	rho[1]=rho[0];
//	rho[2]=0.97;
//	rho[3]=rho[2];
	
	if(errormodel==NONE){
		for(i=0;i<3;i++){
			e[i]=0.;
		}
//	e[0]=0.3;
//	e[1]=0.2;
//	e[2]=0.5;
	}
	else {
		for(i=0;i<3;i++){
			e[i]=0.001;
		}
	}
	calcQ(Q,e[0],e[1],e[2]);		
	
	// Calculate conditional probabilities per site
	CondSiteProbs(ncontigs,npolysites,polysite,Q);
	CondSegProbs(ncontigs,npolysites,rho);
		
	//initial likelihood
	if (mode == SITE) {
		loglik=totalsiteloglik(ncontigs,npolysites,pi,rho);
	}
	else {
		loglik=totalcontigloglik(ncontigs,npolysites,pi,rho);
	}
	
	fprintf(stdout,"Initial log-likelihood: %Lf\n",loglik);
	
	it=0;
	plateausteps=0;
	while(plateausteps<10){
		it++;
		fprintf(stdout,"Iteration %d: ",it);

		//E-step
		for (k=0;k<ncontigs;k++) {
			if(npolysites[k]>0) {
			// XY detailed types
			for(t=0;t<npolysites[k];t++){
					maxl=-INFINITY;
					lmax=-1;
					j=2;
					for(l=0;l<NRTYPES;l++) {
					templ[l]=logl(condsiteprob[k][t][j+l]);
//					templ[l]=condsiteprob[k][t][j+l];
						if(templ[l]>maxl){
							maxl=templ[l];
							lmax=l;
						}
					}
					loghorner(NRTYPES,lmax,templ,rho,expA[k][t]);				
//					horner(NRTYPES,lmax,templ,rho,expA[k][t]);				
			}
			if (mode==SITE) {		//site-wise
				for(t=0;t<npolysites[k];t++){
					//Expectations of segregation types
					maxj=-INFINITY;
					jmax=-1;
					for(j=0;j<NTYPES;j++) {
						temp[j]=logl(condsegprob[k][t][j]);
//						temp[j]=condsegprob[k][t][j];
						if(temp[j]>maxj){
							maxj=temp[j];
							jmax=j;
						}
					}
					loghorner(NTYPES,jmax,temp,pi,expS[k][t]);
//					horner(NTYPES,jmax,temp,pi,expS[k][t]);
				}
			}
			else { // mode == CONTIG
				maxj=-INFINITY;
				jmax=-1;
				for(j=0;j<NTYPES;j++) {
					temp[j]=0; // for log calculations
//					temp[j]=1; // not log
					for(t=0;t<npolysites[k];t++){
						temp[j]+=logl(condsegprob[k][t][j]);
//						temp[j]*=condsegprob[k][t][j];
					}
					if(temp[j]>maxj){
						maxj=temp[j];
						jmax=j;
					}			
				}	
				loghorner(NTYPES,jmax,temp,pi,expR[k]);
//				horner(NTYPES,jmax,temp,pi,expR[k]);
//				maxj=-1.;
//				jmax=-1;
//				for(j=0;j<NTYPES;j++) {
//					temp[j]=1.;
//					for(t=0;t<npolysites[k];t++){
//						temp[j]*=condsegprob[k][t][j];
//					}
//					if(temp[j]>maxj){
//						maxj=temp[j];
//						jmax=j;
//					}			
//				}
//				horner(NTYPES,jmax,temp,pi,expR[k]);
			}
			//Expectations of true genotypes (independent of site- or contig-wise mode)
			if(errormodel==ERRORS) {
				for(t=0;t<npolysites[k];t++){
					for(jl=0;jl<NATYPES;jl++) {
						for(s=0;s<2;s++){
							for(g=0;g<3;g++){
								maxgp=-1;
								gpmax=-1;
								sumgp=0;
								for(gp=0;gp<3;gp++){
									temp[gp]=(long double)P[k][t][s][jl][gp];
									//								sumgp+=temp[gp];
									if(temp[gp]>maxgp){
										maxgp=temp[gp];
										gpmax=gp;
									}
								}
								horner(3,gpmax,temp,Q[g],expTG[k][t][s][jl][g]);	
							}
						}
					}
				}
			}
			}
		}
		
//		if(mode == SITE){
//			fprintf(stdout,"\nE-step: %Lf\n",condsiteexpect(ncontigs,npolysites,polysite,pi,rho,Q));
//		}
		
		//M-step
		if(chromosomes!=NONE){
		//Calculate new values for pi
		if (mode==SITE) {		//site-wise
			sumpi=0;
			for(j=0;j<NTYPES;j++) {
				oldpi[j]=pi[j];
				sumET=0;
				totsites=0;
				for (k=0;k<ncontigs;k++) {
					if(npolysites[k]>0){
					for(t=0;t<npolysites[k];t++){
						sumET+=expS[k][t][j];
					}
					totsites+=npolysites[k];
					}
				}
				pi[j]=(long double)sumET/(long double)totsites;
				sumpi+=pi[j];
			}
		}
		else { // mode == CONTIG
			for(j=0;j<NTYPES;j++) {
				oldpi[j]=pi[j];
			}
			sumpi=0;
			for(j=0;j<NTYPES;j++){ 
				sumET=0;
				for (k=0;k<ncontigs;k++) {
					if (npolysites[k]>0) {					
						sumET+=expR[k][j];
					}
				}
				pi[j]=sumET/(ncontigs-nnoncontigs);
				sumpi+=pi[j];
			}
		}
		
		//Calculate new values for rho
		for(l=0;l<NRTYPES;l++){
			oldrho[l]=rho[l];
		}
		sumET=0;
		weights=0;
		for (k=0;k<ncontigs;k++) {
			if (npolysites[k]>0) {					
				for(t=0;t<npolysites[k];t++){
					if(mode==SITE){
						sumET+=(expA[k][t][0]+expA[k][t][1])*expS[k][t][2];
						weights+=expS[k][t][2];
					}
					else {
						sumET+=(expA[k][t][0]+expA[k][t][1])*expR[k][2];
						weights+=expR[k][2];
					}
				}
			}
		}
		rho[0]=(long double)sumET/(long double)(2.*weights);
		rho[1]=rho[0];
		rho[2]=0.5-rho[0];
		rho[3]=rho[2];
		sumrho=0;
		for(l=0;l<NRTYPES;l++){
			sumrho+=rho[l];
		}
		}
		//Calculate new values for e
		if(errormodel==ERRORS){
			for(i=0;i<3;i++){
				olde[i]=e[i];
			}
			for(g=0;g<3;g++){
				//				printf("\n");
				for(gp=0;gp<3;gp++){
					U[g][gp]=0;
					for (k=0;k<ncontigs;k++) {
						if(npolysites[k]>0){
							for(t=0;t<npolysites[k];t++){
								for(s=0;s<2;s++) {
									temp[s]=0;
									for(j=0;j<NTYPES;j++) {
										if(j<2){
											if(mode==SITE){
												temp[s]+=(long double)expTG[k][t][s][j][g][gp]*(long double)expS[k][t][j];
											}
											else {
												temp[s]+=(long double)expTG[k][t][s][j][g][gp]*(long double)expR[k][j];
											}
										}
										else {
											for(l=0;l<NRTYPES;l++){
												if(mode==SITE){
													temp[s]+=(long double)expTG[k][t][s][j+l][g][gp]*(long double)expS[k][t][j]*(long double)expA[k][t][l];
												}
												else {
													temp[s]+=(long double)expTG[k][t][s][j+l][g][gp]*(long double)expR[k][j]*(long double)expA[k][t][l];
												}
											}
										}
									}
									temp[s]*=(long double)polysite[k][t][1+g+3*s];
								}
								U[g][gp]+=temp[0]+temp[1];
							}
						}
					}
					//					printf("%f ",U[g][gp]);
				}
				
			}
			//			printf("\n");
			
			e[0]=(U[1][0]+U[1][2])/(U[1][0]+U[1][2]+U[0][2]+U[2][0]+U[0][0]+U[2][2]);
			e[1]=(U[0][1]+U[2][1])/(2.*(U[0][1]+U[1][1]+U[2][1]));
			e[2]=(U[0][2]+U[2][0])/(U[1][0]+U[1][2]+U[0][2]+U[2][0]+U[0][0]+U[2][2]);
			//Recalculate Q
			calcQ(Q,e[0],e[1],e[2]);
		}
		
		//evaluate convergence
		deltamax=0.0;
		if(chromosomes!=NONE){
		fprintf(stdout,"pi");
		for(j=0;j<NTYPES;j++) {
			pidelta[j]= pi[j]> minimumvalue ? (pi[j]-oldpi[j])/oldpi[j] : 0.0;
			if (fabsl(pidelta[j]) > fabsl(deltamax)) {
				deltamax=pidelta[j];
			}
			fprintf(stdout," %f (%Le)",pi[j],pidelta[j]);
		}	
		fprintf(stdout,"; total: %f",sumpi);
		fprintf(stdout,"; rho");
		for(l=0;l<NRTYPES;l++) {
			rhodelta[l]= rho[l]> minimumvalue ? (rho[l]-oldrho[l])/oldrho[l] : 0.0;
			if (fabsl(rhodelta[l]) > fabsl(deltamax)) {
				deltamax=rhodelta[l];
			}
			if(l==0 || l==2){
				fprintf(stdout," %f (%Le)",rho[l],rhodelta[l]);
			}
		}	
		fprintf(stdout,"; total: %f",sumrho);
		}
		if(errormodel==ERRORS){
			for(i=0;i<3;i++){
				edelta[i]= e[i]> minimumvalue ? (e[i]-olde[i])/olde[i] : 0.0;
				if (fabsl(edelta[i]) > fabsl(deltamax)) {
					deltamax=edelta[i];
				}
			}
			fprintf(stdout,"; e0: %e; e1: %e; e2: %e",e[0],e[1],e[2]);
		}
		
		if (fabsl(deltamax) < stop) {
			plateausteps++;
		}
		
		// Calculate conditional probabilities per site
		CondSiteProbs(ncontigs,npolysites,polysite,Q);
		CondSegProbs(ncontigs,npolysites,rho);

		//new log-likelihood
		oldloglik=loglik;
		if(mode == SITE){
			loglik=totalsiteloglik(ncontigs,npolysites,pi,rho);
		}
		else {
			loglik=totalcontigloglik(ncontigs,npolysites,pi,rho);
		}
		fprintf(stdout,", log-likelihood: %Lf\n",loglik);

//		if(mode == SITE){
//		fprintf(stdout,"M-step: %Lf\n",condsiteexpect(ncontigs,npolysites,polysite,pi,rho,Q));
//		}
	}		
	//End of EM algorithm. Outputting
	
	if (mode == SITE ) {
		//calculate posterior probabilities per contig 
		for (k=0;k<ncontigs;k++) {
			if(npolysites[k]>0) {
				//autosomal and X-hemizygous: calculate average expectation
				for(j=0;j<2;j++) {
					sumET=0;
					for (t=0; t<npolysites[k]; t++){
						sumET+=expS[k][t][j];
					}
					expR[k][j]=sumET/npolysites[k];
				}
				//XY: calculate average of summed probabilities of all types
				sumET=0;
				for (t=0; t<npolysites[k]; t++){
					for(j=2;j<NTYPES;j++) {
						sumET+=expS[k][t][j];
					}
				}
				expR[k][2]=sumET/npolysites[k];
			}
		}
	}
	
	fprintf(outfile,"pi");
	for(j=0;j<NTYPES;j++) {
		fprintf(outfile," %f ",pi[j]);
	}	
	fprintf(outfile,"; rho");
	for(l=0;l<NRTYPES;l++) {
		if(l==0 || l==2){
			fprintf(outfile," %f ",rho[l]);
		}
	}
	fprintf(outfile,"; e0: %e; e1: %e; e2: %e",e[0],e[1],e[2]);
	fprintf(outfile,", log-likelihood: %Lf",loglik);
	if(chromosomes!=NONE){
		if(errormodel==ERRORS){
			fprintf(outfile,", BIC: %Lf\n",-2.*loglik+6.*log(nI));
		}
		else {
			fprintf(outfile,", BIC: %Lf\n",-2.*loglik+3.*log(nI));
		}
	}
	else {
		if(errormodel==ERRORS){
			fprintf(outfile,", BIC: %Lf\n",-2.*loglik+3.*log(nI));
		}
		else {
			fprintf(outfile,", BIC: %Lf\n",-2.*loglik+0.*log(nI));
		}
	}
	for (k=0;k<ncontigs;k++) {
		fprintf(outfile,">%s\t%d",&contig[k*NAME_LEN],npolysites[k]);
		ni=0;
		for (t=0; t<npolysites[k]; t++){
			for (i=1; i<7; i++) {
				ni+=polysite[k][t][i];
			}
		}
		fprintf(outfile,"\t%f",(double)ni/(double)npolysites[k]);
		if(npolysites[k]>0) {
			pmax=0.0;
			lmax=-1;
			for(j=0;j<NTYPES;j++) {
				if (expR[k][j]>pmax) {
					lmax=j;
					pmax=expR[k][j];
				}
				fprintf(outfile,"\tj : %d; pp : %f",j,expR[k][j]);
			}
			fprintf(outfile,"\t%d\n",lmax);
			//		if(mode==SITE){
			for (t=0; t<npolysites[k]; t++){
				for (i=0; i<7; i++) {
					fprintf(outfile,"%d\t",polysite[k][t][i]);
				}
					pmax=0.0;
					jmax=-1;
				if(mode==SITE){
					for(j=0;j<NTYPES;j++) {
						if (expS[k][t][j]>pmax) {
							jmax=j;
							pmax=expS[k][t][j];
						}
						fprintf(outfile,"%f\t",expS[k][t][j]);
					}
					for(l=0;l<NRTYPES;l++){
						fprintf(outfile,"%f\t",expA[k][t][l]);
					}
					fprintf(outfile,"%d\n",jmax);
				}
				else {
					for(jl=0;jl<NATYPES;jl++) {
						fprintf(outfile,"%Le\t",condsiteprob[k][t][jl]);
						if (condsiteprob[k][t][jl]>pmax) {
							jmax=jl;
							pmax=condsiteprob[k][t][jl];
						}
					}
					fprintf(outfile,"%d\n",jmax);
				}
			}
			//		}
		}
		else {
			fprintf(outfile,"\n");
		}
	}
	
	fclose(outfile);

	freeEM(ncontigs, npolysites);
	free(contig);
	for(k=0;k<ncontigs;k++) {
		if(npolysites[k]>0) {
			for(t=0; t<npolysites[k]; t++){
				free(polysite[k][t]);
			}
			free(polysite[k]);
		}
	}
	free(polysite);
	free(npolysites);

	return 0;
	
}
