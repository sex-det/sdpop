#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
//#include <gsl/gsl_randist.h>
#include <gsl/gsl_machine.h>
#include <fenv.h>
#include <time.h> 
#include "reading.h"
#include "calc.h"

#include <string>
#include <vector>


double ***F; //vector containing allele frequencies per segregation type
double *****P; //vector containing P matrices (per sex and segregation type)
long double ***condsiteprob,***condsegprob; //conditional probabilities per site
long double **contigllik; //conditional log likelihood per contig
double ***expS,***expA,**expR,******expTG;

extern int NAME_LEN;

void initEM(int ncontigs, int *npolysites, int ***polysite) {
	int n11f,n12f,n22f,n11m,n12m,n22m,nfem,nmal,ntot;
	int k,t,jl,s,g;
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
	if((contigllik=(long double **)calloc((size_t)ncontigs,sizeof(long double *)))==NULL) { 
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
			if((expR[k]=(double *)calloc((size_t)JTYPES,sizeof(double)))==NULL) { 
				fprintf(stderr,"error in memory allocation\n");
				exit(1);
			}
			if((contigllik[k]=(long double *)calloc((size_t)JTYPES,sizeof(long double)))==NULL) { 
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
				if((F[k][t]= (double *)calloc((size_t)JLTYPES,sizeof(double)))==NULL) { 
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}
				if((expS[k][t]= (double *)calloc((size_t)JTYPES,sizeof(double)))==NULL) { 
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}
				if((expA[k][t]= (double *)calloc((size_t)LTYPES,sizeof(double)))==NULL) { 
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}
				if((condsiteprob[k][t]= (long double *)calloc((size_t)JLTYPES,sizeof(long double)))==NULL) { 
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}
				if((condsegprob[k][t]= (long double *)calloc((size_t)JTYPES,sizeof(long double)))==NULL) { 
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}
				if((expTG[k][t]= (double ****)calloc((size_t)SEXES,sizeof(double ***)))==NULL) { 
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}
				for(s=0;s<SEXES;s++){
					if((P[k][t][s]= (double **)calloc((size_t)JLTYPES,sizeof(double *)))==NULL) { 
						fprintf(stderr,"error in memory allocation\n");
						exit(1);
					}
					if((expTG[k][t][s]= (double ***)calloc((size_t)JLTYPES,sizeof(double **)))==NULL) { 
						fprintf(stderr,"error in memory allocation\n");
						exit(1);
					}
					for(jl=0;jl<JLTYPES;jl++){
						if((P[k][t][s][jl]= (double *)calloc((size_t)3,sizeof(double)))==NULL) { 
							fprintf(stderr,"error in memory allocation\n");
							exit(1);
						}
						if((expTG[k][t][s][jl]= (double **)calloc((size_t)3,sizeof(double *)))==NULL) { 
							fprintf(stderr,"error in memory allocation\n");
							exit(1);
						}
						for(g=0;g<3;g++){
							if((expTG[k][t][s][jl][g]= (double *)calloc((size_t)3,sizeof(double)))==NULL) { 
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
				F[k][t][jl]=f;
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
				F[k][t][jl]=f;
				P[k][t][FEMALE][jl][0]=0.;
				P[k][t][FEMALE][jl][1]=f;
				P[k][t][FEMALE][jl][2]=1.-f;
				P[k][t][MALE][jl][0]=P[k][t][FEMALE][jl][0];
				P[k][t][MALE][jl][1]=P[k][t][FEMALE][jl][1];
				P[k][t][MALE][jl][2]=P[k][t][FEMALE][jl][2];
				
				//x-hemizygous snps. f is the frequency of allele 1 (symmetry)
				jl=JL_HEMI;
//				f=(double)(2*n11f+n12f+n11m)/(double)(2*nfem+nmal);
				f=(double)(2*n11f+n12f+n11m)/(double)(2*nfem+n11m+n22m);
				F[k][t][jl]=f;
				P[k][t][FEMALE][jl][0]=f*f;
				P[k][t][FEMALE][jl][1]=2.*f*(1.-f);
				P[k][t][FEMALE][jl][2]=(1.-f)*(1.-f);
				P[k][t][MALE][jl][0]=f;
				P[k][t][MALE][jl][1]=0.;
				P[k][t][MALE][jl][2]=1-f;
				
				//x/y snsp ; x-polymorphism
				//allele 1 is fixed on Y; f is the frequency of allele 2 on X
				jl=JL_SEX1;
				f=(double)(2*n22f+n12f+n12m)/(double)(2*nfem+n11m+n12m);
				F[k][t][jl]=f;
				P[k][t][FEMALE][jl][0]=(1.-f)*(1.-f);
				P[k][t][FEMALE][jl][1]=2.*f*(1.-f);
				P[k][t][FEMALE][jl][2]=f*f;
				P[k][t][MALE][jl][0]=1.-f;
				P[k][t][MALE][jl][1]=f;
				P[k][t][MALE][jl][2]=0.;
				//allele 2 is fixed on Y; f is the frequency of allele 1 on X
				jl=JL_SEX2;
				f=(double)(2*n11f+n12f+n12m)/(double)(2*nfem+n12m+n22m);
				F[k][t][jl]=f;
				P[k][t][FEMALE][jl][0]=f*f;
				P[k][t][FEMALE][jl][1]=2.*f*(1.-f);
				P[k][t][FEMALE][jl][2]=(1.-f)*(1.-f);
				P[k][t][MALE][jl][0]=0.;
				P[k][t][MALE][jl][1]=f;
				P[k][t][MALE][jl][2]=1.-f;
				//x/y snsp ; y-polymorphism. f is the frequency, on Y, of the allele that is not fixed on X
				//case 1: allele 1 is fixed on X; f3 is the frequency of allele 2 on Y
				jl=JL_SEX3;
//				f= (n11m+n12m>0) ? (double)(n12m)/(double)(n11m+n12m) : 0;
				f=(double)(n12m+n22m)/(double)(nmal);
				F[k][t][jl]=f;
				P[k][t][FEMALE][jl][0]=1.;
				P[k][t][FEMALE][jl][1]=0.;
				P[k][t][FEMALE][jl][2]=0.;
				P[k][t][MALE][jl][0]=1.-f;
				P[k][t][MALE][jl][1]=f;
				P[k][t][MALE][jl][2]=0.;
				//case 2: allele 2 is fixed on X; f is the frequency of allele 1 on Y
				jl=JL_SEX4;
//				f= (n22m+n12m>0) ? (double)(n12m)/(double)(n22m+n12m) : 0;
				f=(double)(n12m+n11m)/(double)(nmal);
				F[k][t][jl]=f;
				P[k][t][FEMALE][jl][0]=0.;
				P[k][t][FEMALE][jl][1]=0.;
				P[k][t][FEMALE][jl][2]=1.;
				P[k][t][MALE][jl][0]=0.;
				P[k][t][MALE][jl][1]=f;
				P[k][t][MALE][jl][2]=1.-f;
			}
		}		
	}
}

void freeEM(int ncontigs, int *npolysites) {
	int k,t,s,j,g;
	
	for(k=0; k<ncontigs; k++){
		for (t=0; t<npolysites[k]; t++){
			for(s=0;s<SEXES;s++){
				for(j=0;j<JTYPES;j++){
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


long double totalcontigloglik(int ncontigs, int *npolysites, double *pi, double *rho) {
	long double loglik=0;
	long double suml,sumj,temp[JTYPES];
	int t,j,k,l;
	
	for (k=0;k<ncontigs;k++) {
		sumj=0;
		for(j=0;j<JTYPES;j++) {
			if(pi[j]>GSL_DBL_MIN){
				temp[j]=log(pi[j]);
				switch (j) {
				case J_AUTO : case J_HAPLOID :
					for(t=0;t<npolysites[k];t++){
						temp[j]+=logl(condsiteprob[k][t][j]);
					}
					break;
				case J_PARA :
					for(t=0;t<npolysites[k];t++){
						suml=0;
						for(l=0;l<2;l++){
							suml+=(long double)0.5*condsiteprob[k][t][JL_PARA1+l];
						}
						temp[j]+=logl(suml);
					}
					break;
				case J_HEMI :
					for(t=0;t<npolysites[k];t++){
						temp[j]+=logl(condsiteprob[k][t][JL_HEMI]);
					}
					break;
				case J_SEX :
					for(t=0;t<npolysites[k];t++){
						suml=0;
						for(l=0;l<LTYPES;l++){
							suml+=(long double)rho[l]*condsiteprob[k][t][JL_SEX1+l];
						}
						temp[j]+=logl(suml);
					}
					break;
				}
				sumj+=expl(temp[j]);
			}
		}
/*		if(isinf(logl(sumj))){
			fprintf(stdout,"Infinity produced: contig %d, %d sites, %Le\n",k,npolysites[k],sumj);
			for(j=0;j<JTYPES;j++) {
				fprintf(stdout,"%d %Le %Le\n",j,temp[j],expl(temp[j]));
			}
		}
		*/
		if(sumj>LDBL_MIN){
		loglik+=logl(sumj);
		}
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
			for(j=0;j<JTYPES;j++) {
				switch (j) {
				case J_AUTO : case J_HAPLOID :
					sumj+=condsiteprob[k][t][j]*(long double)pi[j];
					break;
				case J_PARA :
					sumj+=(condsiteprob[k][t][JL_PARA1]+condsiteprob[k][t][JL_PARA2])*0.5*(long double)pi[j];
					break;
				case J_HEMI :
					sumj+=condsiteprob[k][t][JL_HEMI]*(long double)pi[j];
					break;
				case J_SEX :
					sum3=0;
					for(l=0;l<LTYPES;l++){
						sum3+=condsiteprob[k][t][JL_SEX1+l]*(long double)rho[l];
					}
					sumj+=sum3*(long double)pi[j];
					break;
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
					for(jl=0;jl<JLTYPES;jl++){
						tempgp=0;
						for(gp=0;gp<3;gp++){
							if(P[k][t][s][jl][gp]>GSL_DBL_MIN){
								tempgp+=expTG[k][t][s][jl][g][gp]*(log(Q[g][gp])+log(P[k][t][s][jl][gp]));
							}
						}
						if(jl<JL_SEX1){
							tempjl+=tempgp*expS[k][t][jl];
						}
						else {
							tempjl+=tempgp*expS[k][t][jl]*expA[k][t][jl-JL_SEX1];
						}
					}					
				tempgs+=polysite[k][t][1+g+3*s]*tempjl;
				}
			}
			logexp+=tempgs;
			for(l=0;l<LTYPES;l++){
				logexp+=expS[k][t][J_SEX]*expA[k][t][l]*log(rho[l]);
			}
			for(j=0;j<JTYPES;j++){
				logexp+=expS[k][t][j]*log(pi[j]);
			}
		}
	}
	return logexp;
}

int contigCONTIGprobs(int nj,int nt,long double **csp,long double *temp)
{
	int j,jmax,t;
	long double maxj;
				maxj=-INFINITY;
				jmax=-1;
				for(j=0;j<nj;j++) {
					temp[j]=0; // for log calculations
//					temp[j]=1; // not log
					for(t=0;t<nt;t++){
						temp[j]+=logl(csp[t][j]);
//						temp[j]*=csp[t][j];
					}
					if(temp[j]>maxj){
						maxj=temp[j];
						jmax=j;
					}			
				}
				return jmax;
}

int SITEprobs(int nj,long double *csp, long double *temp)
{
	int j,jmax;
	double maxj;
	
	maxj=-INFINITY;
	jmax=-1;
	for(j=0;j<nj;j++) {
		temp[j]=logl(csp[j]);
		//						temp[j]=csp[j];
		if(temp[j]>maxj){
			maxj=temp[j];
			jmax=j;
		}
	}
	return jmax;
}

int main(int argc, char *argv[]) {
	
	FILE *fp,*outfile;
	std::vector<std::string> contig;
//	char *contig;
	int i,j,k,l,jl,t,it,s,g,gp;
	int ***polysite;
	int *npolysites;
	int ncontigs,totsites;
	double sumET,weights;     
	double pi[JTYPES],rho[LTYPES],oldpi[JTYPES],oldrho[LTYPES],sumpi,sumrho;
	long double pidelta[JTYPES],rhodelta[LTYPES],edelta[3],deltamax;
	int mode=CONTIG,errormodel=ERRORS,plateausteps,nnoncontigs,chromosomes=XY,ploidy=DIP,paralogs=NOPARA;
	double stop=0.0001; //relative difference between parameter values used to signal convergence
	double minimumvalue=1e-10; //if a parameter value falls below this value, stop evaluating its optimisation
	double Q[3][3],e,olde;
	double maxj,maxl,sumgp;
	long double temp[JTYPES],templ[LTYPES],maxgp,sumL;
	int **bmin;
//	long double logexp,sumj,lsumj;
	double pmin,pmax,gmax,fx,fy;
	int	jmin,jmax,lmax,gpmax;
	long double oldloglik,loglik;
	double U[3][3],Usim,Udis;
	int ni,nI=0,npi=1,sites_individuals,npar;
	double *geom_score;
	
//	feenableexcept(FE_INVALID & FE_OVERFLOW);
//	feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);
	
	for(i=0;i<argc;i++) {
		fprintf(stdout,"%s ",argv[i]);
	}
	fprintf(stdout,"\n");
	
	if (argc != 8) {
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
	if(strcmp(argv[3],"c")==0 || strcmp(argv[3],"1")==0) {
		mode=CONTIG;
	}
	else if (strcmp(argv[3],"s")==0 || strcmp(argv[3],"0")==0) {
		mode=SITE;
	}
	else {
		fprintf(stderr,"Usage: %s infile outfile mode errormodel heterogamety ploidy paralogs\n",argv[0]);
		fprintf(stderr,"Mode should be either \"c\" or \"1\" for contig-mode, or \"s\" or \"0\" for site-wise optimisation\n");
		exit(1);	
	}

	e=0.0001;
	
	if(strcmp(argv[4],"n")==0 || strcmp(argv[4],"0")==0) {
		errormodel=NONE;
		e=0.;
	}
	else if (strcmp(argv[4],"e")==0 || strcmp(argv[4],"1")==0) {
		errormodel=ERRORS;
	}
	else if (strcmp(argv[4],"f")==0 || strcmp(argv[4],"2")==0) {
		errormodel=FIXED;
	}
	else {
		fprintf(stderr,"Usage: %s infile outfile mode errormodel heterogamety ploidy paralogs\n",argv[0]);
		fprintf(stderr,"Errormodel should be either \"e\" or \"1\" to estimate errors, \"f\" or \"2\" to use fixed error rates, or \"n\" or \"0\" for no errors\n");
		exit(1);	
	}
	
	pi[0]=1.;
	for(j=1;j<JTYPES;j++){
		pi[j]=0.;
	}	
	if(strcmp(argv[5],"n")==0 || strcmp(argv[5],"0")==0) {
		chromosomes=NONE;
	}
	else if (strcmp(argv[5],"x")==0 || strcmp(argv[5],"1")==0) {
		chromosomes=XY;
		pi[J_HEMI]=1.;
		pi[J_SEX]=1.;
		npi+=2;
	}
	else if (strcmp(argv[5],"z")==0 || strcmp(argv[5],"2")==0) {
		chromosomes=ZW;
		pi[J_HEMI]=1.;
		pi[J_SEX]=1.;
		npi+=2;
	}
	else {
		fprintf(stderr,"Usage: %s infile outfile mode errormodel heterogamety ploidy paralogs\n",argv[0]);
		fprintf(stderr,"Heterogamety should be either be \"x\" or \"1\" for XY type, \"z\" or \"2\" for ZW type, or \"n\" or \"0\" for no sex chromosomes\n");
		exit(1);	
	}

	if(strcmp(argv[6],"d")==0 || strcmp(argv[6],"0")==0) {
		ploidy=DIP;
	}
	else if (strcmp(argv[6],"h")==0 || strcmp(argv[6],"1")==0) {
		ploidy=HAPDIP;
		pi[J_HAPLOID]=1.;
		npi++;
	}
	else {
		fprintf(stderr,"Usage: %s infile outfile mode errormodel heterogamety ploidy paralogs\n",argv[0]);
		fprintf(stderr,"Ploidy should be either be \"d\" or \"0\" for diploid only, \"h\" or \"1\" for including haploid genes\n");
		exit(1);	
	}
	if(strcmp(argv[7],"o")==0 || strcmp(argv[7],"0")==0) {
		paralogs=NOPARA;
	}
	else if (strcmp(argv[7],"p")==0 || strcmp(argv[7],"1")==0) {
		paralogs=PARALOGS;
		pi[J_PARA]=1.;
		npi++;
	}
	else {
		fprintf(stderr,"Usage: %s infile outfile mode errormodel heterogamety ploidy paralogs\n",argv[0]);
		fprintf(stderr,"Paralogs should be either be \"o\" or \"0\" for orthologs only, \"p\" or \"1\" for including paralogy\n");
		exit(1);	
	}
	
	fprintf(stdout,"Some hardcoded or system-dependent values:\n");
	fprintf(stdout,"Minimum positive value for double: %e\n",DBL_MIN);
	fprintf(stdout,"Minimum positive value for double (GSL): %e\n",GSL_DBL_MIN);
	fprintf(stdout,"Minimum positive value for long double: %Le\n",LDBL_MIN);
	fprintf(stdout,"Criterium for convergence: delta < %e\n",stop);
	fprintf(stdout,"\n");
	fprintf(stdout,"Reading...\n");
	ncontigs=read_cnt2(fp,NAME_LEN,chromosomes,contig,&npolysites,&polysite);
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
	
	//EM algorithm
	
	// Allocation
	initEM(ncontigs, npolysites, polysite);
	
	// Initialize parameters
/*	for(j=0;j<JTYPES;j++) {
		if(pi[j]>0.1){
			pi[j]/=(double)npi;
		}
	}
*/
    // Random pi initialization
    int seed = time(NULL);
    srand(seed);
    sumpi=0;
    for(j=0;j<JTYPES;j++) {
		if(pi[j]>0.1){
			pi[j] = rand()+1;
			sumpi+=pi[j];
		}
    }
    for(j=0;j<JTYPES;j++) {
		if(pi[j]>0.1){
			pi[j]/=sumpi;
    	}
    }

	for(l=0;l<LTYPES;l++) {
		rho[l]=1./LTYPES;
//		rho[l]=2./LTYPES;
	}		
	calcQ(Q,e,e,e);		
	
	// Calculate conditional probabilities per site
	CondSiteProbs(ncontigs,npolysites,polysite,Q,P,condsiteprob);
	CondSegProbs(ncontigs,npolysites,rho,condsiteprob,condsegprob);
	k=ncontigs-1;
//	for(t=0;t<npolysites[k];t++){
//		for (i=0; i<7; i++) {
//			fprintf(stdout,"%d\t",polysite[k][t][i]);
//		}
//		for(j=0;j<JTYPES;j++){
//			fprintf(stdout,"%d: %Lf;\t",j,condsegprob[k][t][j]);
//		}
//		fprintf(stdout,"\n");
//	}
	
	
	//initial likelihood
	if (mode == SITE) {
		loglik=totalsiteloglik(ncontigs,npolysites,pi,rho);
	}
	else {
		loglik=totalcontigloglik(ncontigs,npolysites,pi,rho);
	}
	fprintf(stdout,"Initial values:\npi:");
	for(j=0;j<JTYPES;j++) {
		fprintf(stdout," %f",pi[j]);
	}
	if(chromosomes!=NONE){
		fprintf(stdout,"; rho");
		for(l=0;l<LTYPES;l++) {
			if(l==0 || l==2){
				fprintf(stdout," %f (%Le)",rho[l],rhodelta[l]);
			}
		}
	}
	fprintf(stdout,"; e: %e",e);				
	fprintf(stdout,", log-likelihood: %Lf\n",loglik);
	
	it=0;
	plateausteps=0;
	while(plateausteps<10){
		it++;
		fprintf(stdout,"Iteration %d: ",it);

		//E-step
		for (k=0;k<ncontigs;k++) {
			if(npolysites[k]>0) {
			// XY detailed types
			if(pi[J_SEX]>GSL_DBL_MIN){
				for(t=0;t<npolysites[k];t++){
					maxl=-INFINITY;
					lmax=-1;
					for(l=0;l<LTYPES;l++) {
						templ[l]=logl(condsiteprob[k][t][JL_SEX1+l]);
						//templ[l]=condsiteprob[k][t][j+l];
						if(templ[l]>maxl){
							maxl=templ[l];
							lmax=l;
						}
					}
					loghorner(LTYPES,lmax,templ,rho,expA[k][t]);				
					//horner(LTYPES,lmax,templ,rho,expA[k][t]);				
				}
			}
			//segregation types
			if (mode==SITE) {		//site-wise
				for(t=0;t<npolysites[k];t++){
					//Expectations of segregation types
					jmax=SITEprobs(JTYPES,condsegprob[k][t],temp);
					loghorner(JTYPES,jmax,temp,pi,expS[k][t]);			
					for(j=0;j<JTYPES;j++){
						if(isnan(expS[k][t][j])){
								fprintf(stdout,"NaN produced (expS): contig %d, site %d, type %d: %f\n",k,polysite[k][t][0],j,expS[k][t][j]);
								exit(1);
						}
						else if(isinf(expS[k][t][j])){
							fprintf(stdout,"Inf produced (expS): contig %d, site %d, type %d: %f\n",k,polysite[k][t][0],j,expS[k][t][j]);
						}
					}
					//					horner(JTYPES,jmax,temp,pi,expS[k][t]);
				}
			}
			else { // mode == CONTIG
				jmax=contigCONTIGprobs(JTYPES,npolysites[k],condsegprob[k],contigllik[k]);
				loghorner(JTYPES,jmax,contigllik[k],pi,expR[k]);
//				horner(JTYPES,jmax,temp,pi,expR[k]);
//				maxj=-1.;
//				jmax=-1;
//				for(j=0;j<JTYPES;j++) {
//					temp[j]=1.;
//					for(t=0;t<npolysites[k];t++){
//						temp[j]*=condsegprob[k][t][j];
//					}
//					if(temp[j]>maxj){
//						maxj=temp[j];
//						jmax=j;
//					}			
//				}
//				horner(JTYPES,jmax,temp,pi,expR[k]);
					for(j=0;j<JTYPES;j++){
						if(isnan(expR[k][j])){
								fprintf(stdout,"NaN produced (expR): contig %d, type %d: %f\n",k,j,expR[k][j]);
								exit(1);
						}
						else if(isinf(expR[k][j])){
							fprintf(stdout,"Inf produced (expR): contig %d, type %d: %f\n",k,j,expR[k][j]);
						}
					}
			}
			//Expectations of true genotypes (independent of site- or contig-wise mode)
			if(errormodel==ERRORS) {
				for(t=0;t<npolysites[k];t++){
					for(jl=0;jl<JLTYPES;jl++) {
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
								for(gp=0;gp<3;gp++){
									if(isnan(expTG[k][t][s][jl][g][gp])){
											fprintf(stdout,"NaN produced (expTG): k=%d t=%d s=%d jl=%d g=%d gp=%d: %f\n",k,t,s,jl,g,gp,expR[k][j]);
											exit(1);
									}
									else if(isinf(expTG[k][t][s][jl][g][gp])){
										fprintf(stdout,"Inf produced (expTG): k=%d t=%d s=%d jl=%d g=%d gp=%d: %f\n",k,t,s,jl,g,gp,expR[k][j]);
									}
								}
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
		//Calculate new values for pi
		if (mode==SITE) {		//site-wise
			sumpi=0;
			for(j=0;j<JTYPES;j++) {
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
			for(j=0;j<JTYPES;j++) {
				oldpi[j]=pi[j];
			}
			sumpi=0;
			for(j=0;j<JTYPES;j++){ 
				sumET=0;
				for (k=0;k<ncontigs;k++) {
					if (npolysites[k]>0) {					
						sumET+=expR[k][j];
					}
				}
				pi[j]=(long double)sumET/(long double)(ncontigs-nnoncontigs);
				sumpi+=pi[j];
			}
		}
		if (chromosomes!=NONE) {
		//Calculate new values for rho
		for(l=0;l<LTYPES;l++){
			oldrho[l]=rho[l];
		}
		sumET=0;
		weights=0;
		for (k=0;k<ncontigs;k++) {
			if (npolysites[k]>0) {					
				for(t=0;t<npolysites[k];t++){
					if(mode==SITE){
						sumET+=(expA[k][t][0]+expA[k][t][1])*expS[k][t][J_SEX];
						weights+=expS[k][t][J_SEX];
					}
					else {
						sumET+=(expA[k][t][0]+expA[k][t][1])*expR[k][J_SEX];
						weights+=expR[k][J_SEX];
					}
				}
			}
		}
//		rho[0]=(long double)sumET/(long double)(weights);
//		rho[1]=rho[0];
//		rho[2]=1.-rho[0];
//		rho[3]=rho[2];
		rho[0]=(long double)sumET/(long double)(2.*weights);
		rho[1]=rho[0];
		rho[2]=0.5-rho[0];
		rho[3]=rho[2];
		sumrho=0;
		for(l=0;l<LTYPES;l++){
			sumrho+=rho[l];
		}
		}
		//Calculate new values for e
		if(errormodel==ERRORS){
			olde=e;
			for(g=0;g<3;g++){
				//				printf("\n");
				for(gp=0;gp<3;gp++){
					U[g][gp]=0;
					for (k=0;k<ncontigs;k++) {
						if(npolysites[k]>0){
							for(t=0;t<npolysites[k];t++){
								for(s=0;s<2;s++) {
									temp[s]=0;
									for(j=0;j<JTYPES;j++) {
										switch (j) {
										case J_AUTO : case J_HAPLOID :
											if(mode==SITE){
												temp[s]+=(long double)expTG[k][t][s][j][g][gp]*(long double)expS[k][t][j];
											}
											else {
												temp[s]+=(long double)expTG[k][t][s][j][g][gp]*(long double)expR[k][j];
											}
											break;
										case J_PARA :
											if(mode==SITE){
												temp[s]+=(long double)(expTG[k][t][s][JL_PARA1][g][gp]+expTG[k][t][s][JL_PARA2][g][gp])*0.5*(long double)expS[k][t][j];
											}
											else {
												temp[s]+=(long double)(expTG[k][t][s][JL_PARA1][g][gp]+expTG[k][t][s][JL_PARA2][g][gp])*0.5*(long double)expTG[k][t][s][j][g][gp]*(long double)expR[k][j];
											}
											break;
										case J_HEMI :
											if(mode==SITE){
												temp[s]+=(long double)expTG[k][t][s][JL_HEMI][g][gp]*(long double)expS[k][t][j];
											}
											else {
												temp[s]+=(long double)expTG[k][t][s][JL_HEMI][g][gp]*(long double)expR[k][j];
											}
											break;
										case J_SEX :
											for(l=0;l<LTYPES;l++){
												if(mode==SITE){
													temp[s]+=(long double)expTG[k][t][s][JL_SEX1+l][g][gp]*(long double)expS[k][t][j]*(long double)expA[k][t][l];
												}
												else {
													temp[s]+=(long double)expTG[k][t][s][JL_SEX1+l][g][gp]*(long double)expR[k][j]*(long double)expA[k][t][l];
												}
											}
											break;
										}
									}
									temp[s]*=(long double)polysite[k][t][1+g+3*s];
								}
								U[g][gp]+=temp[0]+temp[1];
							}
						}
					}
					//				printf("%d %d %f ",g,gp,U[g][gp]);
				}
				
			}
			//	printf("\n");
			Usim=0;
			Udis=0;
			for(g=0;g<3;g++){
				for(gp=0;gp<3;gp++){
					if(g==gp){
						Usim+=U[g][gp];
					}
					else {
						Udis+=U[g][gp];
					}
				}
			}
			//printf("%f %f\n",Udis,Usim);
			e=Udis/(2*(Udis+Usim));
			if(e<GSL_DBL_MIN){
				e=GSL_DBL_MIN;
			}
			//Recalculate Q
			calcQ(Q,e,e,e);
		}
		
		//evaluate convergence
		deltamax=0.0;
		if(chromosomes!=NONE || ploidy!=DIP || paralogs!=NOPARA){
			fprintf(stdout,"pi");
			for(j=0;j<JTYPES;j++) {
				pidelta[j]= pi[j]> minimumvalue ? (pi[j]-oldpi[j])/oldpi[j] : 0.0;
				if (fabsl(pidelta[j]) > fabsl(deltamax)) {
					deltamax=pidelta[j];
				}
				fprintf(stdout," %f (%Le)",pi[j],pidelta[j]);
			}	
			fprintf(stdout,"; total: %f",sumpi);
			if(chromosomes!=NONE){
				fprintf(stdout,"; rho");
				for(l=0;l<LTYPES;l++) {
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
		}
		if(errormodel==ERRORS){
//			for(i=0;i<3;i++){
//				edelta[i]= e[i]> minimumvalue ? (e[i]-olde[i])/olde[i] : 0.0;
//				if (fabsl(edelta[i]) > fabsl(deltamax)) {
//					deltamax=edelta[i];
//				}
//			}
			fprintf(stdout,"; e: %e",e);
		}
		
		if (fabsl(deltamax) < stop) {
			plateausteps++;
		}
		
		// Calculate conditional probabilities per site
		CondSiteProbs(ncontigs,npolysites,polysite,Q,P,condsiteprob);
		CondSegProbs(ncontigs,npolysites,rho,condsiteprob,condsegprob);

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
	
	if((geom_score=(double *)malloc(sizeof(double)*JTYPES))==NULL) {
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}
	
	fprintf(outfile,"pi");
	for(j=0;j<JTYPES;j++) {
		fprintf(outfile," %f ",pi[j]);
	}	
	fprintf(outfile,"; rho");                                      
	for(l=0;l<LTYPES;l++) {
		if(l==0 || l==2){
			fprintf(outfile," %f ",rho[l]);
		}
	}
	fprintf(outfile,"; e0: %e",e);
	fprintf(outfile,", log-likelihood: %Lf",loglik);
	npar=ploidy+paralogs+(chromosomes!=NONE)*3;
	if(errormodel==ERRORS){
		npar+=1;
	}
	fprintf(outfile,", BIC (sites): %Lf",-2.*loglik+npar*log(totsites));
	fprintf(outfile,", BIC (contigs): %Lf",-2.*loglik+npar*log(ncontigs-nnoncontigs));
	sites_individuals=0;
	for (k=0;k<ncontigs;k++) {
		for (t=0; t<npolysites[k]; t++){
			for (i=1; i<7; i++) {
				sites_individuals+=polysite[k][t][i];
			}			
		}
	}
	fprintf(outfile,", BIC (sites*individuals): %Lf\n",-2.*loglik+npar*log(sites_individuals));

	for (k=0;k<ncontigs;k++) {
		//calculate likelihoods and posterior probabilities per contig 
			if(npolysites[k]>0) {
				if (mode == SITE ) { //calculate mean posterior
					for(j=0;j<JTYPES;j++) {
						sumET=0;
						for (t=0; t<npolysites[k]; t++){
							sumET+=expS[k][t][j];
							//						sumET+=log(condsegprob[k][t][j]/(condsegprob[k][t][bmin[k][t]]));
						}
						expR[k][j]=sumET/npolysites[k];
					}
				}
				//calculate geometric mean scores
				pmax=0.0;
				jmax=-1;
				for(j=0;j<JTYPES;j++) {
					sumL=0;
					for (t=0; t<npolysites[k]; t++){
						sumL+=log(condsegprob[k][t][j]);
					}
					temp[j]=sumL/npolysites[k];
					if(expl(temp[j])>pmax){
							pmax=expl(temp[j]);
							jmax=j;
					}
				}
				loghorner(JTYPES,jmax,temp,pi,geom_score);

//			fprintf(outfile,">%s\t%d",&contig[k*NAME_LEN],npolysites[k]);
			fprintf(outfile,">%s\t%d",contig[k].data(),npolysites[k]);
			ni=0;
			for (t=0; t<npolysites[k]; t++){
				for (i=1; i<7; i++) {
					ni+=polysite[k][t][i];
				}
			}
			fprintf(outfile,"\t%f",(double)ni/(double)npolysites[k]);
				pmax=0.0;
				gmax=0.0;
				jmax=-1;
				lmax=-1;
			for(j=0;j<JTYPES;j++) {
				if(geom_score[j] > gmax){
					gmax=geom_score[j];
					lmax=j;
				}				
				if(expR[k][j] > pmax){
					pmax=expR[k][j];
					jmax=j;
				}				
			}
			for(j=0;j<JTYPES;j++) {
				fprintf(outfile,"\t%d\t%e\t%e",j,expR[k][j],geom_score[j]);
			}
//			fprintf(outfile,"\t%f\n",zscore(k,npolysites[k],polysite[k],rho,Q,lmax,1000,P[k],condsegprob[k]));			
			fprintf(outfile,"\t%d/%d\n",jmax,lmax);
			
			for (t=0; t<npolysites[k]; t++){
				fprintf(outfile,"%d\t",polysite[k][t][0]);
				fprintf(outfile,"%c%c\t",int2DNA(polysite[k][t][7]),int2DNA(polysite[k][t][8]));			
				for (i=1; i<7; i++) {
					fprintf(outfile,"%d\t",polysite[k][t][i]);
				}
				//calculate estimated frequencies of allele 1 on X and Y
				//first method: fx and fy correspond to those matching with the most probable XY segregation subtype
				pmax=0.0;
				lmax=-1;
				for(l=0;l<LTYPES;l++){
						if (expA[k][t][l]>pmax) {
							lmax=l;
							pmax=expA[k][t][l];
						}
				}
				if(lmax==0) {
					fx=1-F[k][t][JL_SEX1+lmax];
					fy=1;
				}
				else if(lmax==1) {
					fx=F[k][t][JL_SEX1+lmax];
					fy=0;
				}
				else if(lmax==2) {
					fx=1;
					fy=1-F[k][t][JL_SEX1+lmax];
				}
				else if(lmax==3) {
					fx=0;
					fy=F[k][t][JL_SEX1+lmax];
				}
				fprintf(outfile,"%f %f %d\t",fx,fy,lmax);
				//second method: fx and fy correspond to those matching with the most probable XY segregation subtype
				fx=expA[k][t][0]*(1-F[k][t][JL_SEX1])+expA[k][t][1]*F[k][t][JL_SEX2]+expA[k][t][2];				
				fy=expA[k][t][0]+expA[k][t][2]*(1-F[k][t][JL_SEX3])+expA[k][t][3]*F[k][t][JL_SEX4];				
				fprintf(outfile,"%f %f\t",fx,fy);
					for(l=0;l<LTYPES;l++){
						fprintf(outfile,"%f\t",expA[k][t][l]);
					}
				
				pmax=0.0;
				jmax=-1;
				if(mode==SITE){
					for(j=0;j<JTYPES;j++) {
						if (expS[k][t][j]>pmax) {
							jmax=j;
							pmax=expS[k][t][j];
						}
						fprintf(outfile,"%Le %f\t",condsegprob[k][t][j],expS[k][t][j]);
					}
					fprintf(outfile,"%d\n",jmax);
				}
				else {
					for(jl=0;jl<JLTYPES;jl++) {
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
	}
	
	fclose(outfile);

	freeEM(ncontigs, npolysites);
//	free(contig);
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
	free(geom_score);

	return 0;
	
}