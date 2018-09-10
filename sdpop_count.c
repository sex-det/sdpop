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
#define N11 0
#define N12 1
#define N22 2
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

int **countgenotypes(int *npos, int nind, int nmin, char *genotypes, int *sex) {
	char nucvec[4] = {0};
	int i,j,randbit;
	int div,n11f,n12f,n22f,n11m,n12m,n22m,nA,nT,nG,nC,nN;
	int **polysite=NULL;
	char pos1,pos2,nuc1,nuc2;
	int n;
		
	n=0;
	for (j=0; j<*npos; j++){
		//count the number of alleles : A, T, G, C, N
		nA=nT=nG=nC=nN=0;
		for (i=0; i<2*nind; i++){ //loop through all chromosomes
			switch (genotypes[j*nind*2+i]) {
			case 'A' :
				nA++;
				break;
			case 'T' :
				nT++;
				break;
			case 'G' :
				nG++;
				break;
			case 'C' :
				nC++;
				break;
			case 'n' :
			case 'N' :
				nN++;
				break;
			default :
				fprintf(stderr,"Error in parsing genotypes at position %d: allowed are only A, T, G, C, and N (or n).\n",j+1);
				exit(1);
			}
		}
		div=0;
		memset(nucvec, 0, strlen(nucvec)*(sizeof nucvec));
		n11f=n22f=n12f=0;
		n11m=n22m=n12m=0;
		if (nA>0) {
			nucvec[div]='A';
			div++;
		}
		if (nT>0) {
			nucvec[div]='T';
			div++;
		}
		if (nG>0) {
			nucvec[div]='G';
			div++;
		}
		if (nC>0) {
			nucvec[div]='C';
			div++;
		}
		if (div == 1) {
			for (i=0; i<nind; i++){
				if(genotypes[j*nind*2+2*i] == nucvec[0] ) {
					if (sex[i]==FEMALE){
						n11f++;
					}
					else if (sex[i]==MALE){
						n11m++;
					}
				}
			}
		}
		else if (div == 2) { //simple polymorphism
			//randomize 
			nuc1=nucvec[0];
			nuc2=nucvec[1];
			randbit=rand()%2;
			nucvec[randbit]=nuc1;
			nucvec[!randbit]=nuc2;
			for (i=0; i<nind; i++){ //loop through all chromosomes
				pos1=genotypes[j*nind*2+2*i];
				pos2=genotypes[j*nind*2+2*i+1];
				if (sex[i]==FEMALE) {
					if (pos1 == pos2 && pos1 == nucvec[0]) { //homozygous
						n11f++;
					}
					else if (pos1 == pos2 && pos1 == nucvec[1]) { //homozygous
						n22f++;
					}
					else if (pos1 != pos2) { //heterozygous
						n12f++;
					}
				}
				else if (sex[i]==MALE) {
					if (pos1 == pos2 && pos1 == nucvec[0]) { //homozygous
						n11m++;
					}
					else if (pos1 == pos2 && pos1 == nucvec[1]) { //homozygous
						n22m++;
					}
					else if (pos1 != pos2) { //heterozygous
						n12m++;
					}
				}
			}
		}
		if (n11f+n12f+n22f>=nmin && n11m+n12m+n22m>=nmin){
			if (n==0) {
				if((polysite=(int **)malloc(sizeof(int *)))==NULL){
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}
				if((polysite[n]=(int *)malloc(sizeof(int)*6))==NULL){
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}
			}
			else {
				if((polysite=realloc(polysite,sizeof(int *)*(n+1)))==NULL){
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}
				if((polysite[n]=(int *)malloc(sizeof(int)*6))==NULL){
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}
			}
			polysite[n][N11F]=n11f; //female counts
			polysite[n][N12F]=n12f; 
			polysite[n][N22F]=n22f; 			
			polysite[n][N11M]=n11m; //male counts
			polysite[n][N12M]=n12m; 
			polysite[n][N22M]=n22m; 	
			n++;
		}
	}

	if(n==0) {
		*npos=0;
		return NULL;
	}
	else {
		*npos=n;
		return polysite;
	}
}

long double **condsiteprobs(int npolysites, int **polysite) {
	int n11f,n12f,n22f,n11m,n12m,n22m,nfem,nmal,ntot;
	int t,jl,s,g;
	double f;
	double ****P;
	long double **condsiteprob;
	
	//P : a npolysites x 2 x 4 x 3 matrix
	//P[t][s][j][gp] : 
	//t: sites
	//s: sexes, 0=homogametic, 1=heterogametic
	//j: segregation types
	//gp: true genotypes, 0="11", 1="12", 2="22"
		
	
	if (npolysites > 0) {
		if((P=(double ****)calloc((size_t)npolysites,sizeof(double ***)))==NULL) { 
			fprintf(stderr,"error in memory allocation\n");
			exit(1);
		}
		if((condsiteprob=(long double **)calloc((size_t)npolysites,sizeof(long double *)))==NULL) { 
			fprintf(stderr,"error in memory allocation\n");
			exit(1);
		}
		for (t=0; t<npolysites; t++){
			if((P[t]= (double ***)calloc((size_t)SEXES,sizeof(double **)))==NULL) { 
				fprintf(stderr,"error in memory allocation\n");
				exit(1);
			}
			if((condsiteprob[t]= (long double *)calloc((size_t)JLTYPES,sizeof(long double)))==NULL) { 
				fprintf(stderr,"error in memory allocation\n");
				exit(1);
			}
			for(s=0;s<SEXES;s++){
				if((P[t][s]= (double **)calloc((size_t)JLTYPES,sizeof(double *)))==NULL) { 
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}
				for(jl=0;jl<JLTYPES;jl++){
					if((P[t][s][jl]= (double *)calloc((size_t)3,sizeof(double)))==NULL) { 
						fprintf(stderr,"error in memory allocation\n");
						exit(1);
					}
				}
			}
		}							
		
		for (t=0;t<npolysites;t++) {
			n11f=polysite[t][N11F]; //female counts
			n12f=polysite[t][N12F]; 
			n22f=polysite[t][N22F]; 			
			n11m=polysite[t][N11M]; //male counts
			n12m=polysite[t][N12M]; 
			n22m=polysite[t][N22M]; 			
			
			nfem=n11f+n12f+n22f;
			nmal=n11m+n12m+n22m;
			ntot=nfem+nmal;
			//				nFem=(nfem>nFem)?nfem:nFem;
			//				nMal=(nmal>nMal)?nmal:nMal;
			
			//autosomal snps. f0 is the frequency of allele 1 (symmetry)
			jl=JL_AUTO;
			f=(double)(2*n11f+2*n11m+n12f+n12m)/(double)(2*ntot);
			P[t][FEMALE][jl][0]=f*f;
			P[t][FEMALE][jl][1]=2.*f*(1.-f);
			P[t][FEMALE][jl][2]=(1.-f)*(1.-f);
			P[t][MALE][jl][0]=P[t][FEMALE][jl][0];
			P[t][MALE][jl][1]=P[t][FEMALE][jl][1];
			P[t][MALE][jl][2]=P[t][FEMALE][jl][2];
			
			//haploid snps (mitochondria, chloroplasts...)
			jl=JL_HAPLOID;
			f=0.;
			if(n11f+n11m>0) {
				f=(double)(n11f+n11m)/(double)(n11f+n11m+n22f+n22m);
			}
			P[t][FEMALE][jl][0]=f;
			P[t][FEMALE][jl][1]=0.;
			P[t][FEMALE][jl][2]=1.-f;
			P[t][MALE][jl][0]=P[t][FEMALE][jl][0];
			P[t][MALE][jl][1]=P[t][FEMALE][jl][1];
			P[t][MALE][jl][2]=P[t][FEMALE][jl][2];
			//paralogous snps
			//allele 1 is fixed in one of the paralogs
			jl=JL_PARA1;
			f=0.;
			if(n12f+n12m>0) {
				f=1.-sqrt(1.-(double)(n12f+n12m)/(double)(n11f+n11m+n12f+n12m));
			}
			P[t][FEMALE][jl][0]=1.-f;
			P[t][FEMALE][jl][1]=f;
			P[t][FEMALE][jl][2]=0.;
			P[t][MALE][jl][0]=P[t][FEMALE][jl][0];
			P[t][MALE][jl][1]=P[t][FEMALE][jl][1];
			P[t][MALE][jl][2]=P[t][FEMALE][jl][2];
			//allele 2 is fixed in one of the paralogs
			jl=JL_PARA2;
			f=0.;
			if(n12f+n12m>0) {
				f=1.-sqrt(1.-(double)(n12f+n12m)/(double)(n22f+n22m+n12f+n12m));
			}
			P[t][FEMALE][jl][0]=0.;
			P[t][FEMALE][jl][1]=f;
			P[t][FEMALE][jl][2]=1.-f;
			P[t][MALE][jl][0]=P[t][FEMALE][jl][0];
			P[t][MALE][jl][1]=P[t][FEMALE][jl][1];
			P[t][MALE][jl][2]=P[t][FEMALE][jl][2];
			
			//x-hemizygous snps. f1 is the frequency of allele 1 (symmetry)
			jl=JL_HEMIX;
			f=(double)(2*n11f+n12f+n11m)/(double)(2*nfem+nmal);
			P[t][FEMALE][jl][0]=f*f;
			P[t][FEMALE][jl][1]=2.*f*(1.-f);
			P[t][FEMALE][jl][2]=(1.-f)*(1.-f);
			P[t][MALE][jl][0]=f;
			P[t][MALE][jl][1]=0.;
			P[t][MALE][jl][2]=1-f;
			
			//x/y snsp ; x-polymorphism
			//allele 1 is fixed on Y; f21 is the frequency of allele 2 on X
			jl=JL_XY1;
			f=(double)(2*n22f+n12f+n12m)/(double)(2*nfem+n11m+n12m);
			P[t][FEMALE][jl][0]=(1.-f)*(1.-f);
			P[t][FEMALE][jl][1]=2.*f*(1.-f);
			P[t][FEMALE][jl][2]=f*f;
			P[t][MALE][jl][0]=1.-f;
			P[t][MALE][jl][1]=f;
			P[t][MALE][jl][2]=0.;
			//allele 2 is fixed on Y; f22 is the frequency of allele 1 on X
			jl=JL_XY2;
			f=(double)(2*n11f+n12f+n12m)/(double)(2*nfem+n12m+n22m);
			P[t][FEMALE][jl][0]=f*f;
			P[t][FEMALE][jl][1]=2.*f*(1.-f);
			P[t][FEMALE][jl][2]=(1.-f)*(1.-f);
			P[t][MALE][jl][0]=0.;
			P[t][MALE][jl][1]=f;
			P[t][MALE][jl][2]=1.-f;
			//x/y snsp ; y-polymorphism. f3 is the frequency, on Y, of the allele that is not fixed on X
			//case 1: allele 1 is fixed on X; f3 is the frequency of allele 2 on Y
			jl=JL_XY3;
			//f=(double)(n12m)/(double)(n11m+n12m);
			f=(double)(n12m)/(double)(nmal);
			P[t][FEMALE][jl][0]=1.;
			P[t][FEMALE][jl][1]=0.;
			P[t][FEMALE][jl][2]=0.;
			P[t][MALE][jl][0]=1.-f;
			P[t][MALE][jl][1]=f;
			P[t][MALE][jl][2]=0.;
			//case 2: allele 2 is fixed on X; f3 is the frequency of allele 1 on Y
			jl=JL_XY4;
			//f=(double)(n12m)/(double)(n22m+n12m);
			f=(double)(n12m)/(double)(nmal);
			P[t][FEMALE][jl][0]=0.;
			P[t][FEMALE][jl][1]=0.;
			P[t][FEMALE][jl][2]=1.;
			P[t][MALE][jl][0]=0.;
			P[t][MALE][jl][1]=f;
			P[t][MALE][jl][2]=1.-f;
			
			//z-hemizygous snps. f1 is the frequency of allele 1 (symmetry)
			jl=JL_HEMIZ;
			f=(double)(2*n11m+n12m+n11f)/(double)(2*nmal+nfem);
			P[t][MALE][jl][0]=f*f;
			P[t][MALE][jl][1]=2.*f*(1.-f);
			P[t][MALE][jl][2]=(1.-f)*(1.-f);
			P[t][FEMALE][jl][0]=f;
			P[t][FEMALE][jl][1]=0.;
			P[t][FEMALE][jl][2]=1-f;
			
			//x/y snsp ; x-polymorphism
			//allele 1 is fixed on Y; f21 is the frequency of allele 2 on X
			jl=JL_ZW1;
			f=(double)(2*n22m+n12m+n12f)/(double)(2*nmal+n11f+n12f);
			P[t][MALE][jl][0]=(1.-f)*(1.-f);
			P[t][MALE][jl][1]=2.*f*(1.-f);
			P[t][MALE][jl][2]=f*f;
			P[t][FEMALE][jl][0]=1.-f;
			P[t][FEMALE][jl][1]=f;
			P[t][FEMALE][jl][2]=0.;
			//allele 2 is fixed on Y; f22 is the frequency of allele 1 on X
			jl=JL_ZW2;
			f=(double)(2*n11m+n12m+n12f)/(double)(2*nmal+n12f+n22f);
			P[t][MALE][jl][0]=f*f;
			P[t][MALE][jl][1]=2.*f*(1.-f);
			P[t][MALE][jl][2]=(1.-f)*(1.-f);
			P[t][FEMALE][jl][0]=0.;
			P[t][FEMALE][jl][1]=f;
			P[t][FEMALE][jl][2]=1.-f;
			//x/y snsp ; y-polymorphism. f3 is the frequency, on Y, of the allele that is not fixed on X
			//case 1: allele 1 is fixed on X; f3 is the frequency of allele 2 on Y
			jl=JL_ZW3;
			//f=(double)(n12m)/(double)(n11m+n12m);
			f=(double)(n12f)/(double)(nfem);
			P[t][MALE][jl][0]=1.;
			P[t][MALE][jl][1]=0.;
			P[t][MALE][jl][2]=0.;
			P[t][FEMALE][jl][0]=1.-f;
			P[t][FEMALE][jl][1]=f;
			P[t][FEMALE][jl][2]=0.;
			//case 2: allele 2 is fixed on X; f3 is the frequency of allele 1 on Y
			jl=JL_ZW4;
			//f=(double)(n12m)/(double)(n22m+n12m);
			f=(double)(n12f)/(double)(nfem);
			P[t][MALE][jl][0]=0.;
			P[t][MALE][jl][1]=0.;
			P[t][MALE][jl][2]=1.;
			P[t][FEMALE][jl][0]=0.;
			P[t][FEMALE][jl][1]=f;
			P[t][FEMALE][jl][2]=1.-f;
			
			for(jl=0;jl<JLTYPES;jl++) {
				condsiteprob[t][jl]=1;
				for(s=0;s<2;s++){
					for(g=0;g<3;g++){
						condsiteprob[t][jl]*=intpow(P[t][s][jl][g],polysite[t][3*s+g]);
					}
				}
			}
		}
	}
	for (t=0; t<npolysites; t++){
		for(s=0;s<SEXES;s++){
			for(jl=0;jl<JLTYPES;jl++){
				free(P[t][s][jl]);
			}
			free(P[t][s]);
		}
		free(P[t]);
	}
	free(P);
	
	return condsiteprob;
}

char **parsenames(char *string,int *nind,int maxlength) {
	int n,i,j,ichar;
	char **names;
	n=1;
	for (ichar=0;ichar<strlen(string);ichar++) {
		if ( string[ichar] == ',' ) {
			n++;
		}
	}
	if((names=(char **)calloc((size_t)n,sizeof(char *)))==NULL) { 
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}
	for (i=0; i<n; i++){
		if((names[i]= (char *)malloc(sizeof(char) * maxlength))==NULL) { 
			fprintf(stderr,"error in memory allocation\n");
			exit(1);
		}
	}
	i=0;
	j=0;
	for (ichar=0;ichar<strlen(string);ichar++) {
		if (string[ichar] != ',') {
			names[i][j]=string[ichar];
			j++;
		}
		else {
			names[i][j]='\0';
			i++;
			j=0;
		}
		if(i>=maxlength) {
			fprintf(stderr,"In parsenames(): individual names seem too long (>%d). Check argument format or reduce length.\n",maxlength);
			exit(1);
		}
	}
	names[i][j]='\0';
	*nind=n;
	return names;
}

void countsites(int *fixedxy, int *fixedzw, int *polyxy, int *polyzw, int *homozygous, int *other, int nsites, int **site){
	int t;			
	int nZW=0,nXY=0,nhom=0,nzw=0,nxy=0,nother=0;
	int n11f,n12f,n22f,n11m,n12m,n22m,nfem,nmal,ntot,n11,n12,n22;
	
	if(nsites>0) {
		for (t=0; t<nsites; t++){
			n11f=site[t][N11F]; //female counts
			n12f=site[t][N12F]; 
			n22f=site[t][N22F]; 			
			n11m=site[t][N11M]; //male counts
			n12m=site[t][N12M]; 
			n22m=site[t][N22M]; 	
			nfem=n11f+n12f+n22f;
			nmal=n11m+n12m+n22m;
			ntot=nfem+nmal;
			n11=n11f+n11m;
			n12=n12f+n12m;
			n22=n22f+n22m;
			
			if(n11==ntot) { //site homozygous
				nhom++;
			}
			//fixed differences
			else if (n12f==nfem && (n11m==nmal || n22m==nmal)) {
				nZW++;
			}
			else if (n12m==nmal && (n11f==nfem || n22f==nfem)) {
				nXY++;
			}
			//non-fixed xy or zw polymorphisms
			else if (n12f>0 && ((n11f==0 && n22m==nmal) || (n22f==0 && n11m==nmal))) {
				nzw++;
			}
			else if (n12m>0  && ((n22m==0 && n11f==nfem) || (n11m==0 && n22f==nfem))) {
				nxy++;
			}
			else {
				nother++;
			}
		}
	}
	*fixedxy=nXY;
	*fixedzw=nZW;
	*polyxy=nxy;
	*polyzw=nzw;
	*homozygous=nhom;
	*other=nother;
}

long double *segprobs(int npolysites, long double **condsiteprob){	
	int j,jl,t;	
	long double segprob[JTYPES];
	long double *contigprob;
	long double pmax;
	
	if((contigprob=(long double *)calloc((size_t)JTYPES,sizeof(long double)))==NULL) { 
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}
	
	if(npolysites>0) {
		for(j=0;j<JTYPES;j++){
			contigprob[j]=0;
		}
		for (t=0; t<npolysites; t++){
			segprob[J_AUTO]=condsiteprob[t][JL_AUTO];
			segprob[J_HAPLOID]=condsiteprob[t][JL_HAPLOID];
			segprob[J_PARA]=(condsiteprob[t][JL_PARA1]>condsiteprob[t][JL_PARA2]) ? condsiteprob[t][JL_PARA1] : condsiteprob[t][JL_PARA2] ;
			segprob[J_XHEMI]=condsiteprob[t][JL_HEMIX];
			pmax=condsiteprob[t][JL_XY1];
			for(jl=JL_XY2;jl<=JL_XY4;jl++){
				if(condsiteprob[t][jl]>pmax){
					pmax=condsiteprob[t][jl];
				}
			}
			segprob[J_XY]=pmax;
			segprob[J_ZHEMI]=condsiteprob[t][JL_HEMIZ];
			pmax=condsiteprob[t][JL_ZW1];
			for(jl=JL_ZW2;jl<=JL_ZW4;jl++){
				if(condsiteprob[t][jl]>pmax){
					pmax=condsiteprob[t][jl];
				}
			}
			segprob[J_ZW]=pmax;
			for(j=0;j<JTYPES;j++){
				contigprob[j]+=segprob[j];
			}
		}
	}
	return contigprob;
}

int main(int argc, char *argv[]) {
	
	FILE *fp,*outfile;
	char *contig;
	int NAME_LEN=100;
	int CUR_MAX = 4095;
	int i,k,t,j,l,ip,n1,n3,n4;
	int **counts;
//	int **permcounts,nperm;
//	long double **siteprobs,**permsiteprobs;
	int nsites,nf,nm,pos,ta,nA,nT,nG,nC,nN,fi,div;
	int ncontigs,totsites;
	int nnoncontigs;
	int ni,nXY,nZW,nxy,nzw,nother,nhom,nmin;
//	long double *contigprobs,*permcontigprobs;
	char *line = calloc((size_t)CUR_MAX,sizeof(char));
	char *tmpline = calloc((size_t)CUR_MAX,sizeof(char));
	char c;
	int count = 0; 
	int length = 0;
	char ch;
	int firstcontig;
	int pipepos,ichar,nfem,nmal,ifem,imal,nplus=0;
	int *sex,*foundsex;
//	int testp[JTYPES],testnxy,testnzw;
//	double ptestp[JTYPES],ptestnxy,ptestnzw;
	char **name,**femname,**malname;
	int *ffound,*mfound,nfound;
	char *genotypes;
//	int ri,tempsex;

	srand(time(NULL));
//	feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);
	
	for(i=0;i<argc;i++) {
		fprintf(stdout,"%s ",argv[i]);
	}
	fprintf(stdout,"\n");
	
	if (argc != 6) {
		fprintf(stdout,"Usage: %s infile outfile namelist1 namelist2 nmin\n",argv[0]);
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

	femname=parsenames(argv[3],&nfem,NAME_LEN);
	malname=parsenames(argv[4],&nmal,NAME_LEN);	

	nmin=atoi(argv[5]);
	
	if((ffound=(int *)calloc((size_t)(nfem),sizeof(int)))==NULL) {
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}
	if((mfound=(int *)calloc((size_t)(nmal),sizeof(int)))==NULL) {
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}

	if((contig=(char *)malloc(sizeof(char)*NAME_LEN))==NULL) { 
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}			
		
	l=0; //line number
	ch='a';	
	firstcontig=1;
	k=0;
	fprintf(stdout,"Reading and outputting...\n");
		
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
			if (firstcontig==0) {	//treat last read contig
				nsites=t;
				if(nsites>0){
				//do the actual work

				//calculate statistic on data
				//count genotypes
				counts=countgenotypes(&nsites,nfound,nmin,genotypes,foundsex);
//				siteprobs=condsiteprobs(nsites,counts);
				//count characteristic patterns
				countsites(&nXY,&nZW,&nxy,&nzw,&nhom,&nother,nsites,counts);
//				contigprobs=segprobs(npolysites,siteprobs);

/*				//permute
				for(j=0;j<JTYPES;j++){
					testp[j]=0;
				}
				testnxy=0;
				testnzw=0;
				for(ip=0;ip<nperm;ip++){
					for (i=0;i<nfound-1;i++) {
						ri=i+rand()/(RAND_MAX/(nfound-i)+1);
						tempsex=foundsex[i];
						foundsex[i]=foundsex[ri];
						foundsex[ri]=tempsex;
					}
					//calculate permuted statistic
					permcounts=countgenotypes(npolysites,nfound,genotypes,foundsex);
					permsiteprobs=condsiteprobs(npolysites,permcounts);
					countfixed(&p_nXY,&p_nZW,npolysites,nfem,nmal,permcounts);
					permcontigprobs=segprobs(npolysites,permsiteprobs);
					
					//if real statistic > permuted statistic
					for(j=0;j<JTYPES;j++){
						if(permcontigprobs[j]>=contigprobs[j]){
							testp[j]++;
						}
					}
					if(p_nXY>=nXY){
						testnxy++;
					}
					if(p_nZW>=nZW){
						testnzw++;
					}
					for(t=0;t<npolysites;t++){
						free(permsiteprobs[t]);
						free(permcounts[t]);
					}
					free(permsiteprobs);
					free(permcounts);
					free(permcontigprobs);
				}
				//p-values
				for(j=0;j<JTYPES;j++){
					ptestp[j]=(double)testp[j]/(double)nperm;
				}
				ptestnxy=(double)testnxy/(double)nperm;
				ptestnzw=(double)testnzw/(double)nperm;	
				
				fprintf(outfile,">%s\t%d\t%d\t%f\t%d\t%f",contig,npolysites,nXY,ptestnxy,nZW,ptestnzw);
				for(j=0;j<JTYPES;j++){
					fprintf(outfile,"\t%Lf\t%f",contigprobs[j],ptestp[j]);
				}
				fprintf(outfile,"\n");
				*/
				fprintf(outfile,">%s\t%d\t%d\t%d\t%d\t%d\t%d\n",contig,nsites,nXY,nxy,nZW,nzw,nother);
			}
				
				totsites+=nsites;
				k++;
				
				if(k % 5000 == 0){
					fprintf(stdout,"%d contigs, %d polymorphic sites, and still working...\n",k,totsites);
				}
				
				t=0;

				for (i=0; i<ni; i++){
					free(name[i]);
				}
				free(name);	
				free(genotypes);
				free(sex);
				free(foundsex);
			}
			
			if (firstcontig==1) {
				firstcontig=0;
				t=0;
			}
			//Start preparing to read a new contig
			if ( strlen(line) >= NAME_LEN ) {
				fprintf(stderr,"Error: a contig name is too long (more than %d characters) on line %d.\n",NAME_LEN-1,l);
				exit(1);
			}
			sscanf(line,">%s",contig);
		}
		else if ( line[0] == 'p' ){
			for (ifem=0;ifem<nfem;ifem++) {
				ffound[ifem]=0;
			}
			for (imal=0;imal<nmal;imal++) {
				mfound[imal]=0;
			}
			//Important : individuals can be lacking in some contigs !
			//we'll have to count the number of individuals first : the number of words separated by tabs, minus one
			ni=0; //we'll have "position\tname1\tname2\t...\tnameN\0 : the number of tabs is the number of individuals.
			i=0;
			while ((c=line[i]) != '\0'){
				if (c == '\t') 
					ni++;
				i++;
			}
			
			if((name=(char **)calloc((size_t)(ni),sizeof(char *)))==NULL) { 
				fprintf(stderr,"error in memory allocation\n");
				exit(1);
			}
			sscanf(line,"position%[^\n]",tmpline);
			strcpy(line,tmpline);
			//find names iteratively ; shorten the line for each name found :
			for (i=0; i<ni; i++){
				if((name[i]= (char *)malloc(sizeof(char) * NAME_LEN))==NULL) { 
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}
				sscanf(line,"\t%s%[^\n]",name[i],tmpline);
				strcpy(line,tmpline);
				//strip the names unnecessary prefixes (cut after last pipe)
				ichar=0;
				pipepos=0;
				while (name[i][ichar]!='\0') {
					if (name[i][ichar]=='|') {
						pipepos=ichar;
					}
					ichar++;
					if(ichar == NAME_LEN) {
						fprintf(stderr,"Error: one of the individual names seems to have reached the maximum length of %d characters (including prefixes).\n",NAME_LEN);
						exit(1);
					}
				}
				if(pipepos>0) {
					for(ichar=pipepos+1;ichar<NAME_LEN;ichar++) {
						name[i][ichar-pipepos-1]=name[i][ichar];
					}
				}
				
			}
			//Find out the sex of the individuals.
			//sex has the sex of all individuals in the *gen file (sex can be male, female, or absent)
			//foundsex has the sex of the individuals we want to study, in the order of the *gen file
			//genotypes only has the genotypes of the individuals we want to study, in the order of the *gen file
			if((sex=(int *)calloc((size_t)(ni),sizeof(int)))==NULL) { 
				fprintf(stderr,"error in memory allocation\n");
				exit(1);
			}
			if((foundsex=(int *)calloc((size_t)(ni),sizeof(int)))==NULL) { 
				fprintf(stderr,"error in memory allocation\n");
				exit(1);
			}
			nf=0;
			nm=0;
			for (i=0; i<ni; i++){			
				sex[i]=-1;
				for (ifem=0;ifem<nfem;ifem++) {
					if (strcmp(name[i],femname[ifem])==0) {
						sex[i]=FEMALE;
						foundsex[nm+nf]=FEMALE;
						ffound[ifem]=1;
						nf++;
					}
				}
				for (imal=0;imal<nmal;imal++) {
					if (strcmp(name[i],malname[imal])==0) {
						sex[i]=MALE;
						foundsex[nm+nf]=MALE;
						mfound[imal]=1;
						nm++;
					}
				}
			}
			nfound=nf+nm;
			if(ni>nfem+nmal){
				if(nplus==0) {
					fprintf(stdout,"Contig %s seems to have more observations than names given:\n",&contig[k*NAME_LEN]);
					fprintf(stdout,"found %d individuals in contig, and %d names in command line\n",ni,nfem+nmal);
					fprintf(stdout,"No sex found for individual(s): ");
					for (i=0; i<ni; i++){
						if ( sex[i]==-1 ) {
							fprintf(stdout,"%s ",name[i]);
						}		
					}
					fprintf(stdout,"\n");
					fprintf(stdout,"Suppressing warnings about following contigs\n");
					nplus=1;
				}
			}
			firstcontig=0;
			t=0;
		}
		else if (nfound>0) { //line contains genotypes of individuals we're interested in
			sscanf(line,"%d%[^\n]",&pos,tmpline);
			strcpy(line,tmpline);
//			if(pos != t+1){
//				fprintf(stderr,"Error in input file line number %d: expecting to read \"%d ...\"\n",l,t+1);
//				fprintf(stderr,"read \"%d\" instead\n",pos);
//				exit(1);
//			}
			if (t==0) {
				if((genotypes=(char *)malloc(sizeof(char)*(nfound)*2))==NULL){
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}
				ta=1;
			}
			else if (t+1>ta) { //if the previous position was not polymorphic, don't add a new site
				if((genotypes=realloc(genotypes,sizeof(char)*(nfound)*2*(ta+1)))==NULL){
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}
				ta++;
			}
			fi=0;
			for (i=0; i<ni; i++){
				if (sex[i]>=0) {
						sscanf(line,"\t%c%c|%*f%[^\n]",&genotypes[t*nfound*2+2*fi],&genotypes[t*nfound*2+2*fi+1],tmpline);
//					sscanf(line,"\t%c|%c%[^\n]",&genotypes[t*nfound*2+2*fi],&genotypes[t*nfound*2+2*fi+1],tmpline);
					fi++;
				}
				else {
						sscanf(line,"\t%*c%*c|%*f%[^\n]",tmpline);
//					sscanf(line,"\t%*c|%*c%[^\n]",tmpline);
				}
				strcpy(line,tmpline);
			}
			//test for polymorphisms here
			nA=nT=nG=nC=nN=0;
			for (i=0; i<2*ni; i++){ //loop through all chromosomes
				switch (genotypes[t*ni*2+i]) {
				case 'A' :
					nA++;
					break;
				case 'T' :
					nT++;
					break;
				case 'G' :
					nG++;
					break;
				case 'C' :
					nC++;
					break;
				case 'n' :
				case 'N' :
					nN++;
					break;
				default :
					fprintf(stderr,"Error in parsing genotypes at position %d: allowed are only A, T, G, C, and N (or n).\n",j+1);
					exit(1);
				}
			}
			div=0;
			if (nA>0) {
				div++;
			}
			if (nT>0) {
				div++;
			}
			if (nG>0) {
				div++;
			}
			if (nC>0) {
				div++;
			}
			if (div <= 2) { //homozygotes or simple polymorphism; retain
				t++;
			}
			else if (div ==3) { //three alleles observed; discard
				n3++;
			}
			else if (div ==4) { //four alleles observed; discard
				n4++;
			}			
		}
		memset(line, '\0', strlen(line)*(sizeof line));
	}
	
	nsites=t;
	ncontigs=k+1;
		
	//All reading has been done : clear memory.
	free(tmpline);
	free(line);
	
	fclose(fp);
	
	fprintf(stdout,"Found %d contigs\n",ncontigs);
	totsites=0;
	nnoncontigs=0;
	for (k=0;k<ncontigs;k++) {
		totsites+=nsites;
		if(nsites==0) {
			nnoncontigs++;
		}
	}
	fprintf(stdout,"...and %d polymorphic sites\n",totsites);
	if(totsites==0){
		fprintf(stdout,"No polymorphic sites found; nothing to do.\n");
		exit(0);
	}
	
	
	fclose(outfile);

//	for(k=0; k<ncontigs; k++){
//		for (t=0; t<npolysites[k]; t++){
//			free(condsiteprob[k][t]);
//		}
//		free(condsiteprob[k]);
//	}
//	free(condsiteprob);

	free(contig);

	return 0;
	
}