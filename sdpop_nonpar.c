#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
//#include <gsl/gsl_randist.h>
#include <gsl/gsl_machine.h>
#include <fenv.h>
#include <time.h>

#define FEMALE 0
#define MALE 1
#define SEXES 2
#define JTYPES 7 //number of biological segregation types
#define J_AUTO 0
#define J_HAPLOID 1
#define J_PARA 2
#define J_XHEMI 3
#define J_XY 4
#define J_ZHEMI 5
#define J_ZW 6
#define JLTYPES 14 //number of detailed segregation types
#define JL_AUTO 0
#define JL_HAPLOID 1
#define JL_PARA1 2
#define JL_PARA2 3
#define JL_HEMIX 4
#define JL_XY1 5
#define JL_XY2 6
#define JL_XY3 7
#define JL_XY4 8
#define JL_HEMIZ 9
#define JL_ZW1 10
#define JL_ZW2 11
#define JL_ZW3 12
#define JL_ZW4 13
#define LTYPES 4 //number of XY segregation types
#define N11 1
#define N12 2
#define N22 3
#define N11F 0
#define N12F 1
#define N22F 2
#define N11M 3
#define N12M 4
#define N22M 5

double ***F; //vector containing allele frequencies per segregation type
double *****P; //vector containing P matrices (per sex and segregation type)
long double ***condsiteprob; //conditional probabilities per site
int nI=0,nFem=0,nMal=0; //maximum number of individuals, females and males

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

int lreturnmax(long double *pmax, long double *larray, int n) {
	int i;
	int imax=0;
	long double max=0;
	double r;
	
	max=larray[0];
	for (i=1; i<n; i++){
		if (larray[i]>max){
			max=larray[i];
			imax=i;
		}
		else if (larray[i]<max){
			continue;
		}
		else { //in the rare case that both are exactly equal, choose one at random
			for (;;) {
				if((r=(double)rand()/RAND_MAX)<0.5){
					max=larray[i];
					imax=i;
					break;
				}
				else if (r>0.5) {
					break;
				}
			}
		}
	}
	*pmax=max;
	return imax;
	
}

void initEM(int ncontigs, int *npolysites, int ***polysite) {
	int n11f,n12f,n22f,n11m,n12m,n22m,nfem,nmal,ntot;
	int k,t,jl,s;
	double f;
	
	//P : a npolysites x 2 x 4 x 3 matrix
	//P[kt][s][j][gp] : 
	//kt: sites
	//s: sexes, 0=homogametic, 1=heterogametic
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
	if((condsiteprob=(long double ***)calloc((size_t)ncontigs,sizeof(long double **)))==NULL) { 
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
			if((condsiteprob[k]=(long double **)calloc((size_t)npolysites[k],sizeof(long double *)))==NULL) { 
				fprintf(stderr,"error in memory allocation\n");
				exit(1);
			}
			for (t=0; t<npolysites[k]; t++){
				if((P[k][t]= (double ***)calloc((size_t)SEXES,sizeof(double **)))==NULL) { 
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}
				if((F[k][t]= (double *)calloc((size_t)JLTYPES,sizeof(double)))==NULL) { 
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}
				if((condsiteprob[k][t]= (long double *)calloc((size_t)JLTYPES,sizeof(long double)))==NULL) { 
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}
				for(s=0;s<SEXES;s++){
					if((P[k][t][s]= (double **)calloc((size_t)JLTYPES,sizeof(double *)))==NULL) { 
						fprintf(stderr,"error in memory allocation\n");
						exit(1);
					}
					for(jl=0;jl<JLTYPES;jl++){
						if((P[k][t][s][jl]= (double *)calloc((size_t)3,sizeof(double)))==NULL) { 
							fprintf(stderr,"error in memory allocation\n");
							exit(1);
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
				nFem=(nfem>nFem)?nfem:nFem;
				nMal=(nmal>nMal)?nmal:nMal;
				
				//autosomal snps. f0 is the frequency of allele 1 (symmetry)
				jl=JL_AUTO;
				f=(double)(2*n11f+2*n11m+n12f+n12m)/(double)(2*ntot);
				F[k][t][jl]=f;
				P[k][t][FEMALE][jl][0]=f*f;
				P[k][t][FEMALE][jl][1]=2.*f*(1.-f);
				P[k][t][FEMALE][jl][2]=(1.-f)*(1.-f);
				P[k][t][MALE][jl][0]=P[k][t][FEMALE][jl][0];
				P[k][t][MALE][jl][1]=P[k][t][FEMALE][jl][1];
				P[k][t][MALE][jl][2]=P[k][t][FEMALE][jl][2];
				
				//haploid snps (mitochondria, chloroplasts...)
				jl=JL_HAPLOID;
				f=0.;
				if(n11f+n11m>0) {
					f=(double)(n11f+n11m)/(double)(n11f+n11m+n22f+n22m);
				}
				F[k][t][jl]=f;
				P[k][t][FEMALE][jl][0]=f;
				P[k][t][FEMALE][jl][1]=0.;
				P[k][t][FEMALE][jl][2]=1.-f;
				P[k][t][MALE][jl][0]=P[k][t][FEMALE][jl][0];
				P[k][t][MALE][jl][1]=P[k][t][FEMALE][jl][1];
				P[k][t][MALE][jl][2]=P[k][t][FEMALE][jl][2];
				//paralogous snps
				//allele 1 is fixed in one of the paralogs
				jl=JL_PARA1;
				f=0.;
				if(n12f+n12m>0) {
					f=1.-sqrt(1.-(double)(n12f+n12m)/(double)(n11f+n11m+n12f+n12m));
				}
				P[k][t][FEMALE][jl][0]=1.-f;
				P[k][t][FEMALE][jl][1]=f;
				P[k][t][FEMALE][jl][2]=0.;
				P[k][t][MALE][jl][0]=P[k][t][FEMALE][jl][0];
				P[k][t][MALE][jl][1]=P[k][t][FEMALE][jl][1];
				P[k][t][MALE][jl][2]=P[k][t][FEMALE][jl][2];
				//allele 2 is fixed in one of the paralogs
				jl=JL_PARA2;
				f=0.;
				if(n12f+n12m>0) {
					f=1.-sqrt(1.-(double)(n12f+n12m)/(double)(n22f+n22m+n12f+n12m));
				}
				P[k][t][FEMALE][jl][0]=0.;
				P[k][t][FEMALE][jl][1]=f;
				P[k][t][FEMALE][jl][2]=1.-f;
				P[k][t][MALE][jl][0]=P[k][t][FEMALE][jl][0];
				P[k][t][MALE][jl][1]=P[k][t][FEMALE][jl][1];
				P[k][t][MALE][jl][2]=P[k][t][FEMALE][jl][2];
				
				//x-hemizygous snps. f1 is the frequency of allele 1 (symmetry)
				jl=JL_HEMIX;
				f=(double)(2*n11f+n12f+n11m)/(double)(2*nfem+nmal);
				F[k][t][jl]=f;
				P[k][t][FEMALE][jl][0]=f*f;
				P[k][t][FEMALE][jl][1]=2.*f*(1.-f);
				P[k][t][FEMALE][jl][2]=(1.-f)*(1.-f);
				P[k][t][MALE][jl][0]=f;
				P[k][t][MALE][jl][1]=0.;
				P[k][t][MALE][jl][2]=1-f;
				
				//x/y snsp ; x-polymorphism
				//allele 1 is fixed on Y; f21 is the frequency of allele 2 on X
				jl=JL_XY1;
				f=(double)(2*n22f+n12f+n12m)/(double)(2*nfem+n11m+n12m);
				P[k][t][FEMALE][jl][0]=(1.-f)*(1.-f);
				P[k][t][FEMALE][jl][1]=2.*f*(1.-f);
				P[k][t][FEMALE][jl][2]=f*f;
				P[k][t][MALE][jl][0]=1.-f;
				P[k][t][MALE][jl][1]=f;
				P[k][t][MALE][jl][2]=0.;
				//allele 2 is fixed on Y; f22 is the frequency of allele 1 on X
				jl=JL_XY2;
				f=(double)(2*n11f+n12f+n12m)/(double)(2*nfem+n12m+n22m);
				P[k][t][FEMALE][jl][0]=f*f;
				P[k][t][FEMALE][jl][1]=2.*f*(1.-f);
				P[k][t][FEMALE][jl][2]=(1.-f)*(1.-f);
				P[k][t][MALE][jl][0]=0.;
				P[k][t][MALE][jl][1]=f;
				P[k][t][MALE][jl][2]=1.-f;
				//x/y snsp ; y-polymorphism. f3 is the frequency, on Y, of the allele that is not fixed on X
				//case 1: allele 1 is fixed on X; f3 is the frequency of allele 2 on Y
				jl=JL_XY3;
				//f=(double)(n12m)/(double)(n11m+n12m);
				f=(double)(n12m)/(double)(nmal);
				P[k][t][FEMALE][jl][0]=1.;
				P[k][t][FEMALE][jl][1]=0.;
				P[k][t][FEMALE][jl][2]=0.;
				P[k][t][MALE][jl][0]=1.-f;
				P[k][t][MALE][jl][1]=f;
				P[k][t][MALE][jl][2]=0.;
				//case 2: allele 2 is fixed on X; f3 is the frequency of allele 1 on Y
				jl=JL_XY4;
				//f=(double)(n12m)/(double)(n22m+n12m);
				f=(double)(n12m)/(double)(nmal);
				P[k][t][FEMALE][jl][0]=0.;
				P[k][t][FEMALE][jl][1]=0.;
				P[k][t][FEMALE][jl][2]=1.;
				P[k][t][MALE][jl][0]=0.;
				P[k][t][MALE][jl][1]=f;
				P[k][t][MALE][jl][2]=1.-f;

				//z-hemizygous snps. f1 is the frequency of allele 1 (symmetry)
				jl=JL_HEMIZ;
				f=(double)(2*n11m+n12m+n11f)/(double)(2*nmal+nfem);
				F[k][t][jl]=f;
				P[k][t][MALE][jl][0]=f*f;
				P[k][t][MALE][jl][1]=2.*f*(1.-f);
				P[k][t][MALE][jl][2]=(1.-f)*(1.-f);
				P[k][t][FEMALE][jl][0]=f;
				P[k][t][FEMALE][jl][1]=0.;
				P[k][t][FEMALE][jl][2]=1-f;
				
				//x/y snsp ; x-polymorphism
				//allele 1 is fixed on Y; f21 is the frequency of allele 2 on X
				jl=JL_ZW1;
				f=(double)(2*n22m+n12m+n12f)/(double)(2*nmal+n11f+n12f);
				P[k][t][MALE][jl][0]=(1.-f)*(1.-f);
				P[k][t][MALE][jl][1]=2.*f*(1.-f);
				P[k][t][MALE][jl][2]=f*f;
				P[k][t][FEMALE][jl][0]=1.-f;
				P[k][t][FEMALE][jl][1]=f;
				P[k][t][FEMALE][jl][2]=0.;
				//allele 2 is fixed on Y; f22 is the frequency of allele 1 on X
				jl=JL_ZW2;
				f=(double)(2*n11m+n12m+n12f)/(double)(2*nmal+n12f+n22f);
				P[k][t][MALE][jl][0]=f*f;
				P[k][t][MALE][jl][1]=2.*f*(1.-f);
				P[k][t][MALE][jl][2]=(1.-f)*(1.-f);
				P[k][t][FEMALE][jl][0]=0.;
				P[k][t][FEMALE][jl][1]=f;
				P[k][t][FEMALE][jl][2]=1.-f;
				//x/y snsp ; y-polymorphism. f3 is the frequency, on Y, of the allele that is not fixed on X
				//case 1: allele 1 is fixed on X; f3 is the frequency of allele 2 on Y
				jl=JL_ZW3;
				//f=(double)(n12m)/(double)(n11m+n12m);
				f=(double)(n12f)/(double)(nfem);
				P[k][t][MALE][jl][0]=1.;
				P[k][t][MALE][jl][1]=0.;
				P[k][t][MALE][jl][2]=0.;
				P[k][t][FEMALE][jl][0]=1.-f;
				P[k][t][FEMALE][jl][1]=f;
				P[k][t][FEMALE][jl][2]=0.;
				//case 2: allele 2 is fixed on X; f3 is the frequency of allele 1 on Y
				jl=JL_ZW4;
				//f=(double)(n12m)/(double)(n22m+n12m);
				f=(double)(n12f)/(double)(nfem);
				P[k][t][MALE][jl][0]=0.;
				P[k][t][MALE][jl][1]=0.;
				P[k][t][MALE][jl][2]=1.;
				P[k][t][FEMALE][jl][0]=0.;
				P[k][t][FEMALE][jl][1]=f;
				P[k][t][FEMALE][jl][2]=1.-f;
			}
		}		
	}
}

void freeEM(int ncontigs, int *npolysites) {
	int k,t,s,j;
	
	for(k=0; k<ncontigs; k++){
		for (t=0; t<npolysites[k]; t++){
			for(s=0;s<SEXES;s++){
				for(j=0;j<JLTYPES;j++){
					free(P[k][t][s][j]);
				}
				free(P[k][t][s]);			
			}
			free(P[k][t]);
			free(F[k][t]);
			free(condsiteprob[k][t]);
		}
		free(P[k]);
		free(F[k]);
		free(condsiteprob[k]);
	}
	free(P);
	free(F);
	free(condsiteprob);
}

void CondSiteProbs(int ncontigs, int *npolysites, int ***polysite)
{
	int k,t,jl,s,g;
	
	for (k=0;k<ncontigs;k++) {
		for(t=0;t<npolysites[k];t++){
			for(jl=0;jl<JLTYPES;jl++) {
				condsiteprob[k][t][jl]=1;
				for(s=0;s<2;s++){
					for(g=0;g<3;g++){
//						condsiteprob[k][t][jl]*=intpow(P[k][t][s][jl][g],polysite[k][t][1+3*s+g]);
						condsiteprob[k][t][jl]*=intpow(P[k][t][s][jl][g],polysite[k][t][3*s+g]);
					}
				}
			}
		}
	}
}


int read_cnt(FILE *fp, 	int NAME_LEN, int chromosomes, char **contig_p, int **npolysites_p, char ****site_p, int ****polysite_p) 
{
	char *contig;
	int *npolysites,***polysite;
	char ***site;
	int CUR_MAX = 4095;
	int NCONTIG_BATCH = 100;
	int ncontigs_allocated,ncontigs;
	char *line = calloc((size_t)CUR_MAX,sizeof(char));
	char *tmpline = calloc((size_t)CUR_MAX,sizeof(char));
	int count = 0; 
	int length = 0;
	char ch;
	int i,k,l,firstcontig,t;
	int ni;

	//reading counts from file 
	l=0; //line number
	ch='a';
	if((npolysites=(int *)malloc(sizeof(int)*NCONTIG_BATCH))==NULL) {
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}
	if((site=(char ***)malloc(sizeof(char **)*NCONTIG_BATCH))==NULL) {
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
				if((site=(char ***)realloc(site,sizeof(char **)*ncontigs_allocated))==NULL) {
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
					if((site[k]=(char **)malloc(sizeof(int *)))==NULL){
						fprintf(stderr,"error in memory allocation\n");
						exit(1);
					}
					if((site[k][t]=(char *)malloc(sizeof(int)*NAME_LEN))==NULL){
						fprintf(stderr,"error in memory allocation\n");
						exit(1);
					}
					if((polysite[k]=(int **)malloc(sizeof(int *)))==NULL){
						fprintf(stderr,"error in memory allocation\n");
						exit(1);
					}
					if((polysite[k][t]=(int *)malloc(sizeof(int)*6))==NULL){
						fprintf(stderr,"error in memory allocation\n");
						exit(1);
					}
				}
				else {
					if((site[k]=realloc(site[k],sizeof(char *)*(t+1)))==NULL){
						fprintf(stderr,"error in memory allocation\n");
						exit(1);
					}
					if((site[k][t]=(char *)malloc(sizeof(char)*NAME_LEN))==NULL){
						fprintf(stderr,"error in memory allocation\n");
						exit(1);
					}
					if((polysite[k]=realloc(polysite[k],sizeof(int *)*(t+1)))==NULL){
						fprintf(stderr,"error in memory allocation\n");
						exit(1);
					}
					if((polysite[k][t]=(int *)malloc(sizeof(int)*6))==NULL){
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
					sscanf(line,"%s\t%d\t%d\t%d\t%d\t%d\t%d\t",site[k][t],&polysite[k][t][0],&polysite[k][t][1],&polysite[k][t][2],&polysite[k][t][3],&polysite[k][t][4],&polysite[k][t][5]);
				ni=0;
				for(i=0;i<6;i++){
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
	
	*npolysites_p=npolysites;
	*polysite_p=polysite;
	*site_p=site;
	*contig_p=contig;
	
	return ncontigs;
}

int main(int argc, char *argv[]) {
	
	FILE *fp,*outfile;
	char *contig;
	int NAME_LEN=100;
	int i,k,jl,t,j;
	int ***polysite;
	char ***site;
	int *npolysites;
	int ncontigs,totsites;
	int nnoncontigs;
	long double pmax;
	int	jmax;
	int ni,nXY,nZW;
	long double contigprob[JTYPES];
	long double segprob[JTYPES];
	int maxprob[JTYPES];

	srand(time(NULL));
//	feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);
	
	for(i=0;i<argc;i++) {
		fprintf(stdout,"%s ",argv[i]);
	}
	fprintf(stdout,"\n");
	
	if (argc != 3) {
		fprintf(stdout,"Usage: %s infile outfile mode errormodel heterogamety ploidy paralogs\n",argv[0]);
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
	
	ncontigs=read_cnt(fp,NAME_LEN,1,&contig,&npolysites,&site,&polysite);
	fclose(fp);
	
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
	
	// Allocation
	initEM(ncontigs, npolysites, polysite);	
	CondSiteProbs(ncontigs,npolysites,polysite);
	
	for (k=0;k<ncontigs;k++) {
		fprintf(outfile,">%s\t%d",&contig[k*NAME_LEN],npolysites[k]);
		ni=0;
		nXY=0;
		nZW=0;
		for (t=0; t<npolysites[k]; t++){
			for (i=0; i<6; i++) {
				ni+=polysite[k][t][i];
			}
		}
		fprintf(outfile,"\t%f",(double)ni/(double)npolysites[k]);
		if(npolysites[k]>0) {
			for(j=0;j<JTYPES;j++){
				contigprob[j]=0;
				maxprob[j]=0;
			}
			for (t=0; t<npolysites[k]; t++){
				segprob[J_AUTO]=condsiteprob[k][t][JL_AUTO];
				segprob[J_HAPLOID]=condsiteprob[k][t][JL_HAPLOID];
				segprob[J_PARA]=(condsiteprob[k][t][JL_PARA1]>condsiteprob[k][t][JL_PARA2]) ? condsiteprob[k][t][JL_PARA1] : condsiteprob[k][t][JL_PARA2] ;
				segprob[J_XHEMI]=condsiteprob[k][t][JL_HEMIX];
				pmax=condsiteprob[k][t][JL_XY1];
				for(jl=JL_XY2;jl<=JL_XY4;jl++){
					if(condsiteprob[k][t][jl]>pmax){
						pmax=condsiteprob[k][t][jl];
					}
				}
				segprob[J_XY]=pmax;
				segprob[J_ZHEMI]=condsiteprob[k][t][JL_HEMIZ];
				pmax=condsiteprob[k][t][JL_ZW1];
				for(jl=JL_ZW2;jl<=JL_ZW4;jl++){
					if(condsiteprob[k][t][jl]>pmax){
						pmax=condsiteprob[k][t][jl];
					}
				}
				segprob[J_ZW]=pmax;
				for(j=0;j<JTYPES;j++){
					contigprob[j]+=segprob[j];
				}
				jmax=lreturnmax(&pmax,segprob,JTYPES);
//				pmax=segprob[0];
//				jmax=0;
//				for (j=1; j<JTYPES; j++){
//					if (segprob[j]>pmax){
//						pmax=segprob[j];
//						jmax=j;
//					}
//					else if (segprob[j]==pmax){
//						r=(double)rand()/RAND_MAX;
//						while(r==0.5) {
//							r=(double)rand()/RAND_MAX;
//						}
//						if(r<0.5){
//							pmax=segprob[j];
//							jmax=j;
//						}
//					}
//				}
				maxprob[jmax]++;
				//count fixed XY or ZW differences
				if (polysite[k][t][N12F]==nFem && (polysite[k][t][N11M]==nMal || polysite[k][t][N22M]==nMal)) {
					nZW++;
				}
				else if (polysite[k][t][N12M]==nMal && (polysite[k][t][N11F]==nFem || polysite[k][t][N22F]==nFem)) {
					nXY++;
				}
				

			}
//			pmax=0;
//			jmax=0;
			for(j=0;j<JTYPES;j++){
				contigprob[j]/=(long double)npolysites[k];
				fprintf(outfile,"\t%Lf\t%d",contigprob[j],maxprob[j]);
//				if(contigprob[j]>=pmax){
//					pmax=contigprob[j];
//					jmax=j;
//				}
			}
			fprintf(outfile,"\t%d\t%d",nXY,nZW);
			jmax=lreturnmax(&pmax,contigprob,JTYPES);
			fprintf(outfile,"\t%d\t%Lf",jmax,logl(pmax)-logl(contigprob[J_AUTO]));
		}
		fprintf(outfile,"\n");
		
		//permutation
		
		
		if(npolysites[k]>0) {
			//		if(mode==SITE){
			for (t=0; t<npolysites[k]; t++){
					fprintf(outfile,"%s\t",site[k][t]);
				for (i=0; i<6; i++) {
					fprintf(outfile,"%d\t",polysite[k][t][i]);
				}
//					pmax=0.0;
//					jmax=-1;

					for(jl=0;jl<JLTYPES;jl++) {
//						if(condsiteprob[k][t][jl]>GSL_DBL_MIN){
//						fprintf(outfile,"%Lf\t",logl(condsiteprob[k][t][jl]));
//						}
//						else {
							fprintf(outfile,"%Lf\t",condsiteprob[k][t][jl]);
//						}
//						if (condsiteprob[k][t][jl]>pmax) {
//							jmax=jl;
//							pmax=condsiteprob[k][t][jl];
//						}
					}
					jmax=lreturnmax(&pmax,condsiteprob[k][t],JLTYPES);
					fprintf(outfile,"%d\t%Lf\n",jmax,logl(pmax)-logl(condsiteprob[k][t][0]));
				}
			}
	}
	
	fclose(outfile);

	freeEM(ncontigs, npolysites);
	free(contig);
	for(k=0;k<ncontigs;k++) {
		if(npolysites[k]>0) {
			for(t=0; t<npolysites[k]; t++){
				free(polysite[k][t]);
				free(site[k][t]);
			}
			free(polysite[k]);
			free(site[k]);
		}
	}
	free(polysite);
	free(site);
	free(npolysites);

	return 0;
	
}