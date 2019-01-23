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
double ***expS,****expA,**expR,******expTG;

Model model;

extern int NAME_LEN;

void initEM(std::vector<Contig>& contigs) {
	int n11f,n12f,n22f,n11m,n12m,n22m,nfem,nmal,ntot;
	int k,t,jl,s,g,j;
	int ncontigs,npolysites;
	double f;
	
	ncontigs=contigs.size();
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
	if((expA=(double ****)calloc((size_t)ncontigs,sizeof(double ***)))==NULL) { 
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
	
	
	for (k=0;k<contigs.size();k++) {
		Contig & current_contig = contigs[k];
		if((npolysites=current_contig.snps.size())>0) {
			if((P[k]=(double ****)calloc((size_t)npolysites,sizeof(double ***)))==NULL) { 
				fprintf(stderr,"error in memory allocation\n");
				exit(1);
			}
			if((F[k]=(double **)calloc((size_t)npolysites,sizeof(double *)))==NULL) { 
				fprintf(stderr,"error in memory allocation\n");
				exit(1);
			}
			if((expS[k]=(double **)calloc((size_t)npolysites,sizeof(double *)))==NULL) { 
				fprintf(stderr,"error in memory allocation\n");
				exit(1);
			}
			if((expA[k]=(double ***)calloc((size_t)npolysites,sizeof(double **)))==NULL) { 
				fprintf(stderr,"error in memory allocation\n");
				exit(1);
			}
			if((condsiteprob[k]=(long double **)calloc((size_t)npolysites,sizeof(long double *)))==NULL) { 
				fprintf(stderr,"error in memory allocation\n");
				exit(1);
			}
			if((condsegprob[k]=(long double **)calloc((size_t)npolysites,sizeof(long double *)))==NULL) { 
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
			if((expTG[k]=(double *****)calloc((size_t)npolysites,sizeof(double ****)))==NULL) { 
				fprintf(stderr,"error in memory allocation\n");
				exit(1);
			}
			for (t=0; t<npolysites; t++){
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
				if((expA[k][t]= (double **)calloc((size_t)JTYPES,sizeof(double *)))==NULL) { 
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
				for(j=0;j<JTYPES;j++){
					if((expA[k][t][j]= (double *)calloc((size_t)4,sizeof(double)))==NULL) { 
						fprintf(stderr,"error in memory allocation\n");
						exit(1);
					}
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
			
			for (t=0;t<npolysites;t++) {
				n11f=current_contig.snps[t].genotypes_by_sex[N11F]; //female counts
				n12f=current_contig.snps[t].genotypes_by_sex[N12F]; 
				n22f=current_contig.snps[t].genotypes_by_sex[N22F]; 			
				n11m=current_contig.snps[t].genotypes_by_sex[N11M]; //male counts
				n12m=current_contig.snps[t].genotypes_by_sex[N12M]; 
				n22m=current_contig.snps[t].genotypes_by_sex[N22M]; 			
				
				nfem=n11f+n12f+n22f;
				nmal=n11m+n12m+n22m;
				ntot=nfem+nmal;
				
				foreach_jl(model,[&](const auto jl){
						//autosomal snps. f0 is the frequency of allele 1 (symmetry)
						if(jl==JL_AUTO){
							f=(double)(2*n11f+2*n11m+n12f+n12m)/(double)(2*ntot);
							F[k][t][jl]=f;
							P[k][t][FEMALE][jl][0]=f*f;
							P[k][t][FEMALE][jl][1]=2.*f*(1.-f);
							P[k][t][FEMALE][jl][2]=(1.-f)*(1.-f);
							P[k][t][MALE][jl][0]=P[k][t][FEMALE][jl][0];
							P[k][t][MALE][jl][1]=P[k][t][FEMALE][jl][1];
							P[k][t][MALE][jl][2]=P[k][t][FEMALE][jl][2];
						}
						//haploid snps (mitochondria, chloroplasts...)
						if(jl==JL_HAPLOID){
//							f=0.;
//							if(n11f+n11m>0) {
//								f=(double)(n11f+n11m)/(double)(n11f+n11m+n22f+n22m);
//							}
							f=(double)(n11f+n12f+n11m+n12m)/(double)(ntot);
							F[k][t][jl]=f;
							P[k][t][FEMALE][jl][0]=f;
							P[k][t][FEMALE][jl][1]=0.;
							P[k][t][FEMALE][jl][2]=1.-f;
							P[k][t][MALE][jl][0]=P[k][t][FEMALE][jl][0];
							P[k][t][MALE][jl][1]=P[k][t][FEMALE][jl][1];
							P[k][t][MALE][jl][2]=P[k][t][FEMALE][jl][2];
						}
						//paralogous snps
						//allele 1 is fixed in one of the paralogs
						if(jl==JL_PARA1){
//							f=0.;
//							if(n12f+n12m>0) {
//								f=1.-sqrt(1.-(double)(n12f+n12m)/(double)(n11f+n11m+n12f+n12m));
//							}
							f=1.-sqrt(1.-(double)(n12f+n12m+n22f+n22m)/(double)(ntot));
							F[k][t][jl]=f;
							P[k][t][FEMALE][jl][0]=1.-f;
							P[k][t][FEMALE][jl][1]=f;
							P[k][t][FEMALE][jl][2]=0.;
							P[k][t][MALE][jl][0]=P[k][t][FEMALE][jl][0];
							P[k][t][MALE][jl][1]=P[k][t][FEMALE][jl][1];
							P[k][t][MALE][jl][2]=P[k][t][FEMALE][jl][2];
						}
						//allele 2 is fixed in one of the paralogs
						if(jl==JL_PARA2){
//							f=0.;
//							if(n12f+n12m>0) {
//								f=1.-sqrt(1.-(double)(n12f+n12m)/(double)(n22f+n22m+n12f+n12m));
//							}
							f=1.-sqrt(1.-(double)(n12f+n12m+n11f+n11m)/(double)(ntot));
							F[k][t][jl]=f;
							P[k][t][FEMALE][jl][0]=0.;
							P[k][t][FEMALE][jl][1]=f;
							P[k][t][FEMALE][jl][2]=1.-f;
							P[k][t][MALE][jl][0]=P[k][t][FEMALE][jl][0];
							P[k][t][MALE][jl][1]=P[k][t][FEMALE][jl][1];
							P[k][t][MALE][jl][2]=P[k][t][FEMALE][jl][2];
						}
						//x-hemizygous snps. f is the frequency of allele 1 (symmetry)
						if(jl==JL_HEMI){
//							f=(double)(2*n11f+n12f+n11m)/(double)(2*nfem+n11m+n22m);
							f=(double)(2*n11f+n12f+n11m+n12m)/(double)(2*nfem+nmal);
							F[k][t][jl]=f;
							P[k][t][FEMALE][jl][0]=f*f;
							P[k][t][FEMALE][jl][1]=2.*f*(1.-f);
							P[k][t][FEMALE][jl][2]=(1.-f)*(1.-f);
							P[k][t][MALE][jl][0]=f;
							P[k][t][MALE][jl][1]=0.;
							P[k][t][MALE][jl][2]=1-f;
						}
						//x/y snsp ; x-polymorphism
						//allele 1 is fixed on Y; f is the frequency of allele 2 on X
						if(jl==JL_SEX1){
//							f=(double)(2*n22f+n12f+n12m)/(double)(2*nfem+n11m+n12m);
							f=(double)(2*n22f+n12f+n12m+n22m)/(double)(2*nfem+nmal);
							F[k][t][jl]=f;
							P[k][t][FEMALE][jl][0]=(1.-f)*(1.-f);
							P[k][t][FEMALE][jl][1]=2.*f*(1.-f);
							P[k][t][FEMALE][jl][2]=f*f;
							P[k][t][MALE][jl][0]=1.-f;
							P[k][t][MALE][jl][1]=f;
							P[k][t][MALE][jl][2]=0.;
						}
						//allele 2 is fixed on Y; f is the frequency of allele 1 on X
						if(jl==JL_SEX2){
//							f=(double)(2*n11f+n12f+n12m)/(double)(2*nfem+n12m+n22m);
							f=(double)(2*n11f+n12f+n12m+n11m)/(double)(2*nfem+nmal);
							F[k][t][jl]=f;
							P[k][t][FEMALE][jl][0]=f*f;
							P[k][t][FEMALE][jl][1]=2.*f*(1.-f);
							P[k][t][FEMALE][jl][2]=(1.-f)*(1.-f);
							P[k][t][MALE][jl][0]=0.;
							P[k][t][MALE][jl][1]=f;
							P[k][t][MALE][jl][2]=1.-f;
						}
						//x/y snsp ; y-polymorphism. f is the frequency, on Y, of the allele that is not fixed on X
						//case 1: allele 1 is fixed on X; f3 is the frequency of allele 2 on Y
						if(jl==JL_SEX3){
							//				f= (n11m+n12m>0) ? (double)(n12m)/(double)(n11m+n12m) : 0;
							f=(double)(n12m+n22m)/(double)(nmal);
							F[k][t][jl]=f;
							P[k][t][FEMALE][jl][0]=1.;
							P[k][t][FEMALE][jl][1]=0.;
							P[k][t][FEMALE][jl][2]=0.;
							P[k][t][MALE][jl][0]=1.-f;
							P[k][t][MALE][jl][1]=f;
							P[k][t][MALE][jl][2]=0.;
						}
						//case 2: allele 2 is fixed on X; f is the frequency of allele 1 on Y
						if(jl==JL_SEX4){
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
						//z-hemizygous snps. f is the frequency of allele 1 (symmetry)
						if(jl==JL_ZHEMI){
//							f=(double)(2*n11m+n12m+n11f)/(double)(2*nmal+n11f+n22f);
							f=(double)(2*n11m+n12m+n11f+n12f)/(double)(2*nmal+nfem);
							F[k][t][jl]=f;
							P[k][t][MALE][jl][0]=f*f;
							P[k][t][MALE][jl][1]=2.*f*(1.-f);
							P[k][t][MALE][jl][2]=(1.-f)*(1.-f);
							P[k][t][FEMALE][jl][0]=f;
							P[k][t][FEMALE][jl][1]=0.;
							P[k][t][FEMALE][jl][2]=1-f;
						}
						//z/w snsp ; z-polymorphism
						//allele 1 is fixed on W; f is the frequency of allele 2 on Z
						if(jl==JL_ZW1){
//							f=(double)(2*n22m+n12m+n12f)/(double)(2*nmal+n11f+n12f);
							f=(double)(2*n22m+n12m+n12f+n22f)/(double)(2*nmal+nfem);
							F[k][t][jl]=f;
							P[k][t][MALE][jl][0]=(1.-f)*(1.-f);
							P[k][t][MALE][jl][1]=2.*f*(1.-f);
							P[k][t][MALE][jl][2]=f*f;
							P[k][t][FEMALE][jl][0]=1.-f;
							P[k][t][FEMALE][jl][1]=f;
							P[k][t][FEMALE][jl][2]=0.;
						}
						//allele 2 is fixed on W; f is the frequency of allele 1 on Z
						if(jl==JL_ZW2){
//							f=(double)(2*n11m+n12m+n12f)/(double)(2*nmal+n12f+n22f);
							f=(double)(2*n11m+n12m+n12f+n11f)/(double)(2*nmal+nfem);
							F[k][t][jl]=f;
							P[k][t][MALE][jl][0]=f*f;
							P[k][t][MALE][jl][1]=2.*f*(1.-f);
							P[k][t][MALE][jl][2]=(1.-f)*(1.-f);
							P[k][t][FEMALE][jl][0]=0.;
							P[k][t][FEMALE][jl][1]=f;
							P[k][t][FEMALE][jl][2]=1.-f;
						}
						//z/w snsp ; w-polymorphism. f is the frequency, on W, of the allele that is not fixed on Z
						//case 1: allele 1 is fixed on Z; f is the frequency of allele 2 on W
						if(jl==JL_ZW3){
							f=(double)(n12f+n22f)/(double)(nfem);
							F[k][t][jl]=f;
							P[k][t][MALE][jl][0]=1.;
							P[k][t][MALE][jl][1]=0.;
							P[k][t][MALE][jl][2]=0.;
							P[k][t][FEMALE][jl][0]=1.-f;
							P[k][t][FEMALE][jl][1]=f;
							P[k][t][FEMALE][jl][2]=0.;
						}
						//case 2: allele 2 is fixed on Z; f is the frequency of allele 1 on W
						if(jl==JL_ZW4){
							f=(double)(n12f+n11f)/(double)(nfem);
							F[k][t][jl]=f;
							P[k][t][MALE][jl][0]=0.;
							P[k][t][MALE][jl][1]=0.;
							P[k][t][MALE][jl][2]=1.;
							P[k][t][FEMALE][jl][0]=0.;
							P[k][t][FEMALE][jl][1]=f;
							P[k][t][FEMALE][jl][2]=1.-f;
						}
				});
			}
		}		
	}
}

void freeEM(std::vector<Contig>& contigs) {
	int k,t,s,jl,j,g;
	int npolysites;
	
	for(k=0; k<contigs.size(); k++){
		Contig & current_contig = contigs[k];
		npolysites=current_contig.snps.size();
		for (t=0; t<npolysites; t++){
			for(s=0;s<SEXES;s++){
				for(jl=0;jl<JLTYPES;jl++){
					free(P[k][t][s][jl]);
				}
				for(g=0;g<3;g++){
					free(expTG[k][t][s][g]);
				}
				free(P[k][t][s]);			
				free(expTG[k][t][s]);			
			}
			for(j=0;j<JTYPES;j++){
				free(expA[k][t][j]);
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

long double totalcontigloglik(std::vector<Contig>& contigs, double *pi, double **rho) {
	long double loglik=0;
	long double suml,sumj,temp[JTYPES];
	int t,k,l;
	
	for (k=0;k<contigs.size();k++) {
		Contig & current_contig = contigs[k];
		int npolysites=current_contig.snps.size();
		sumj=0;
		foreach_j(model,[&](const auto j){
				//for(j=0;j<JTYPES;j++) {
				if(pi[j]>GSL_DBL_MIN){
					temp[j]=log(pi[j]);
					if( j==J_AUTO || j== J_HAPLOID){
						for(t=0;t<npolysites;t++){
							temp[j]+=logl(condsiteprob[k][t][j]);
						}
					}
					if(j== J_PARA){
						for(t=0;t<npolysites;t++){
							suml=0;
							for(l=0;l<2;l++){
								suml+=(long double)rho[j][l]*condsiteprob[k][t][JL_PARA1+l];
							}
							temp[j]+=logl(suml);
						}
					}
					if(j== J_HEMI){
						for(t=0;t<npolysites;t++){
							temp[j]+=logl(condsiteprob[k][t][JL_HEMI]);
						}
					}
					if( j== J_SEX){
						for(t=0;t<npolysites;t++){
							suml=0;
							foreach_l_xy(model,[&](const auto l,const auto jl,const auto j){
									//for(l=0;l<LTYPES;l++){
									suml+=(long double)rho[j][l]*condsiteprob[k][t][jl];
							});
							temp[j]+=logl(suml);
						}
					}
					if(j== J_ZHEMI){
						for(t=0;t<npolysites;t++){
							temp[j]+=logl(condsiteprob[k][t][JL_ZHEMI]);
						}
					}
					if( j== J_ZW){
						for(t=0;t<npolysites;t++){
							suml=0;
							foreach_l_zw(model,[&](const auto l,const auto jl,const auto j){
									//for(l=0;l<LTYPES;l++){
									suml+=(long double)rho[j][l]*condsiteprob[k][t][jl];
							});
							temp[j]+=logl(suml);
						}
					}
					sumj+=expl(temp[j]);
				}
		});
		if(isinf(logl(sumj))){
			fprintf(stdout,"Infinity produced: contig %d, %d sites, %Le\n",k,npolysites,sumj);
			foreach_j(model,[&](const auto j){
				fprintf(stdout,"%d %Le %Le\n",j,temp[j],expl(temp[j]));
			});
		}
		
		if(sumj>LDBL_MIN){
		loglik+=logl(sumj);
		}
	}
	return loglik;
}

long double totalsiteloglik(std::vector<Contig>& contigs, double *pi, double **rho) {
	long double loglik=0;
	long double sumj,sum3;
	int t,k;
	
	for (k=0;k<contigs.size();k++) {
		Contig & current_contig = contigs[k];
		for(t=0;t<current_contig.snps.size();t++){
			sumj=0;
			foreach_j(model,[&](const auto j){
					//for(j=0;j<JTYPES;j++) {
					if(j== J_AUTO || j== J_HAPLOID){
						sumj+=condsiteprob[k][t][j]*(long double)pi[j];
					}
					if(j== J_PARA){
						foreach_l_para(model,[&](const auto l,const auto jl,const auto j){
								sumj+=condsiteprob[k][t][jl]*rho[j][l]*(long double)pi[j];
						});
					}
					if(j== J_HEMI){
						sumj+=condsiteprob[k][t][JL_HEMI]*(long double)pi[j];
					}
					if(j== J_SEX){
						sum3=0;
						foreach_l_xy(model,[&](const auto l,const auto jl,const auto j){
								//for(l=0;l<LTYPES;l++){
								sum3+=condsiteprob[k][t][jl]*(long double)rho[j][l];
						});
						sumj+=sum3*(long double)pi[j];
					}
					if(j== J_ZHEMI){
						sumj+=condsiteprob[k][t][JL_ZHEMI]*(long double)pi[j];
					}
					if(j== J_ZW){
						sum3=0;
						foreach_l_zw(model,[&](const auto l,const auto jl,const auto j){
								//for(l=0;l<LTYPES;l++){
								sum3+=condsiteprob[k][t][jl]*(long double)rho[j][l];
						});
						sumj+=sum3*(long double)pi[j];
					}
			});
			loglik+=logl(sumj);
		}
	}
	return loglik;
}

long double condsiteexpect(std::vector<Contig>& contigs, double *pi, double **rho, double Q[3][3])
{
	long double tempgs,tempjl,tempgp,logexp=0;
	int k,t,s,g,gp;
	
	for (k=0;k<contigs.size();k++) {
		Contig & current_contig = contigs[k];
		for(t=0;t<current_contig.snps.size();t++){
			tempgs=0;
			for(s=0;s<2;s++){
				for(g=0;g<3;g++){
					tempjl=0;
					foreach_jl(model,[&](const auto jl){
						tempgp=0;
						for(gp=0;gp<3;gp++){
							if(P[k][t][s][jl][gp]>GSL_DBL_MIN){
								tempgp+=expTG[k][t][s][jl][g][gp]*(log(Q[g][gp])+log(P[k][t][s][jl][gp]));
							}
						}
						if(jl==JL_AUTO || jl==JL_HAPLOID){
							tempjl+=tempgp*expS[k][t][jl];
						}
						if(jl==JL_HEMI){
							tempjl+=tempgp*expS[k][t][J_HEMI];
						}
						if(jl==JL_ZHEMI){
							tempjl+=tempgp*expS[k][t][J_ZHEMI];
						}
						if(jl == JL_PARA1 || jl == JL_PARA2 ) {
							tempjl+=tempgp*expS[k][t][J_PARA]*expA[k][t][J_PARA][jl-JL_PARA1];
						}						
						if(jl >= JL_SEX1 && jl <= JL_SEX4 ) {
							tempjl+=tempgp*expS[k][t][J_SEX]*expA[k][t][J_SEX][jl-JL_SEX1];
						}
						if(jl >= JL_ZW1 && jl <= JL_ZW4 ) {
							tempjl+=tempgp*expS[k][t][J_ZW]*expA[k][t][J_ZW][jl-JL_ZW1];
						}
					});					
				tempgs+=current_contig.snps[t].genotypes_by_sex[g+3*s]*tempjl;
				}
			}
			logexp+=tempgs;
			foreach_l_xy(model,[&](const auto l,const auto jl, const auto j){
					logexp+=expS[k][t][j]*expA[k][t][j][l]*log(rho[j][l]);
			});
			foreach_l_zw(model,[&](const auto l,const auto jl, const auto j){
					logexp+=expS[k][t][j]*expA[k][t][j][l]*log(rho[j][l]);
			});
			foreach_l_para(model,[&](const auto l,const auto jl, const auto j){
					logexp+=expS[k][t][j]*expA[k][t][j][l]*log(rho[j][l]);
			});
			foreach_j(model,[&](const auto j){
				logexp+=expS[k][t][j]*log(pi[j]);
			});
		}
	}
	return logexp;
}
long double condcontigexpect(std::vector<Contig>& contigs, double *pi, double **rho, double Q[3][3])
{
	long double tempgs,tempjl,tempgp,logexp=0;
	int k,t,s,g,gp;
	
	for (k=0;k<contigs.size();k++) {
		Contig & current_contig = contigs[k];
		for(t=0;t<current_contig.snps.size();t++){
			tempgs=0;
			for(s=0;s<2;s++){
				for(g=0;g<3;g++){
					tempjl=0;
					foreach_jl(model,[&](const auto jl){
						tempgp=0;
						for(gp=0;gp<3;gp++){
							if(P[k][t][s][jl][gp]>GSL_DBL_MIN){
								tempgp+=expTG[k][t][s][jl][g][gp]*(log(Q[g][gp])+log(P[k][t][s][jl][gp]));
							}
						}
						if(jl==JL_AUTO || jl==JL_HAPLOID){
							tempjl+=tempgp*expR[k][jl];
						}
						if(jl==JL_HEMI){
							tempjl+=tempgp*expR[k][J_HEMI];
						}
						if(jl==JL_ZHEMI){
							tempjl+=tempgp*expR[k][J_ZHEMI];
						}
						if(jl == JL_PARA1 || jl == JL_PARA2 ) {
							tempjl+=tempgp*expR[k][J_PARA]*expA[k][t][J_PARA][jl-JL_PARA1];
						}						
						if(jl >= JL_SEX1 && jl <= JL_SEX4 ) {
							tempjl+=tempgp*expR[k][J_SEX]*expA[k][t][J_SEX][jl-JL_SEX1];
						}
						if(jl >= JL_ZW1 && jl <= JL_ZW4 ) {
							tempjl+=tempgp*expR[k][J_ZW]*expA[k][t][J_ZW][jl-JL_ZW1];
						}
					});					
				tempgs+=current_contig.snps[t].genotypes_by_sex[g+3*s]*tempjl;
				}
			}
			logexp+=tempgs;
			foreach_l_xy(model,[&](const auto l,const auto jl, const auto j){
					logexp+=expR[k][j]*expA[k][t][j][l]*log(rho[j][l]);
			});
			foreach_l_zw(model,[&](const auto l,const auto jl, const auto j){
					logexp+=expR[k][j]*expA[k][t][j][l]*log(rho[j][l]);
			});
			foreach_l_para(model,[&](const auto l,const auto jl, const auto j){
					logexp+=expR[k][j]*expA[k][t][j][l]*log(rho[j][l]);
			});
			foreach_j(model,[&](const auto j){
					logexp+=expR[k][j]*log(pi[j]);
			});
		}
//		if(current_contig.snps.size()>0){
//			foreach_j(model,[&](const auto j){
//					logexp+=expR[k][j]*log(pi[j]);
//			});
//		}
	}
	return logexp;
}

int contigCONTIGprobs(int nj,int nt,long double **csp,long double *temp)
{
	int jmax,t;
	long double maxj;
				maxj=-INFINITY;
				jmax=-1;
				foreach_j(model,[&](const auto j){
				//for(j=0;j<nj;j++) {
					temp[j]=0; // for log calculations
//					temp[j]=1; // not log
					for(t=0;t<nt;t++){
						temp[j]+=logl(csp[t][j]);
//						temp[j]*=csp[t][j];
					}
					if(temp[j]>=maxj){
						maxj=temp[j];
						jmax=j;
					}			
				});
				return jmax;
}

int SITEprobs(int nj,long double *csp, long double *temp)
{
	int jmax;
	double maxj;
	
	maxj=-INFINITY;
	jmax=-1;
	foreach_j(model,[&](const auto j){
	//for(j=0;j<nj;j++) {
		temp[j]=logl(csp[j]);
		//						temp[j]=csp[j];
		if(temp[j]>=maxj){
			maxj=temp[j];
			jmax=j;
		}
	});
	return jmax;
}

int main(int argc, char *argv[]) {
	
	FILE *fp,*outfile;
	std::vector<std::string> contig;
//	char *contig;
	int i,j,k,l,t,it,s,g,gp;
	int npolysites;
	int ncontigs,totsites;
	double sumET,weights;     
	double dumpi[JTYPES],pi[JTYPES],**rho,oldpi[JTYPES],**oldrho,sumpi,sumrho;
	long double pidelta[JTYPES],**rhodelta,deltamax;
	int mode=CONTIG,errormodel=ERRORS,plateausteps,nnoncontigs;
	double stop=0.0001; //relative difference between parameter values used to signal convergence
	double minimumvalue=1e-10; //if a parameter value falls below this value, stop evaluating its optimisation
	double Q[3][3],e;
	double maxl,sumgp;
	long double temp[JTYPES],templ[LTYPES],maxgp,sumL;
	double pmax,gmax,rmax,fx,fy;
	int	jmax,lmax,nmax,gpmax;
	long double oldloglik,loglik;
	double U[3][3],Usim,Udis;
	int nj,ni,npi=1,sites_individuals,npar;
	double *geom_score,*geom_nopi_score;
	int warning=0;

	std::vector<Contig> contigs;

	for(j=0;j<JTYPES;j++){
		pi[j]=0;
		oldpi[j]=0;
		pidelta[j]=0;
		temp[j]=0;
	}
	for(l=0;l<4;l++){
		templ[l]=0;
	}
	
	
//	feenableexcept(FE_INVALID & FE_OVERFLOW);
//	feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);
	
	for(i=0;i<argc;i++) {
		fprintf(stdout,"%s ",argv[i]);
	}
	fprintf(stdout,"\n");
	
	if (argc != 8) {
		fprintf(stderr,"Usage: %s infile outfile mode errormodel heterogamety ploidy paralogs\n",argv[0]);
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
	
	if(strcmp(argv[5],"n")==0 || strcmp(argv[5],"0")==0) {
		model.xy=0;
		model.zw=0;
	}
	else if (strcmp(argv[5],"x")==0 || strcmp(argv[5],"1")==0) {
		model.xy=1;
		model.zw=0;
		npi+=2;
	}
	else if (strcmp(argv[5],"z")==0 || strcmp(argv[5],"2")==0) {
		model.xy=0;
		model.zw=1;
		npi+=2;
	}
	else if (strcmp(argv[5],"b")==0 || strcmp(argv[5],"3")==0) {
		model.xy=1;
		model.zw=1;
		npi+=2;
	}
	else {
		fprintf(stderr,"Usage: %s infile outfile mode errormodel heterogamety ploidy paralogs\n",argv[0]);
		fprintf(stderr,"Heterogamety should be either be \"x\" or \"1\" for XY type, \"z\" or \"2\" for ZW type, \"n\" or \"0\" for no sex chromosomes, or \"b\" or \"3\" for both.\n");
		exit(1);	
	}

	if(strcmp(argv[6],"d")==0 || strcmp(argv[6],"0")==0) {
		model.haploid=0;
	}
	else if (strcmp(argv[6],"h")==0 || strcmp(argv[6],"1")==0) {
		model.haploid=1;
		npi++;
	}
	else {
		fprintf(stderr,"Usage: %s infile outfile mode errormodel heterogamety ploidy paralogs\n",argv[0]);
		fprintf(stderr,"Ploidy should be either be \"d\" or \"0\" for diploid only, \"h\" or \"1\" for including haploid genes\n");
		exit(1);	
	}
	if(strcmp(argv[7],"o")==0 || strcmp(argv[7],"0")==0) {
		model.paralogs=0;
	}
	else if (strcmp(argv[7],"p")==0 || strcmp(argv[7],"1")==0) {
		model.paralogs=1;
		npi++;
	}
	else {
		fprintf(stderr,"Usage: %s infile outfile mode errormodel heterogamety ploidy paralogs\n",argv[0]);
		fprintf(stderr,"Paralogs should be either be \"o\" or \"0\" for orthologs only, \"p\" or \"1\" for including paralogy\n");
		exit(1);	
	}
	
	foreach_j(model,[&](const auto j){
			pi[j]=1./(double)npi;
	});
	
	fprintf(stdout,"Some hardcoded or system-dependent values:\n");
	fprintf(stdout,"Minimum positive value for double: %e\n",DBL_MIN);
	fprintf(stdout,"Minimum positive value for double (GSL): %e\n",GSL_DBL_MIN);
	fprintf(stdout,"Minimum positive value for long double: %Le\n",LDBL_MIN);
	fprintf(stdout,"Criterium for convergence: delta < %e\n",stop);
	fprintf(stdout,"\n");
	fprintf(stdout,"Reading...\n");
	ncontigs=read_cnt_model(fp,NAME_LEN,model,contigs);
	fclose(fp);
	
	fprintf(stdout,"Found %d contigs\n",ncontigs);
	totsites=0;
	nnoncontigs=0;
	for (k=0;k<contigs.size();k++) {
		Contig & current_contig = contigs[k];
		npolysites=current_contig.snps.size();
		totsites+=npolysites;
		if(npolysites==0) {
			nnoncontigs++;
		}
	}
	fprintf(stdout,"...and %d polymorphic sites\n",totsites);
	if(totsites==0){
		fprintf(stderr,"No polymorphic sites found; nothing to do.\n");
		exit(0);
	}
	
	//EM algorithm
	
	// Allocation
	initEM(contigs);
	
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
 	foreach_j(model,[&](const auto j){
			pi[j] = rand();
//			pi[j] = 1;
			sumpi+=pi[j];
    });
 	foreach_j(model,[&](const auto j){
			pi[j]/=sumpi;
    });
    
    if((rho=(double **)calloc((size_t)JTYPES,sizeof(double *)))==NULL) { 
    	fprintf(stderr,"error in memory allocation\n");
    	exit(1);
    }
    if((oldrho=(double **)calloc((size_t)JTYPES,sizeof(double *)))==NULL) { 
    	fprintf(stderr,"error in memory allocation\n");
    	exit(1);
    }
    if((rhodelta=(long double **)calloc((size_t)JTYPES,sizeof(long double *)))==NULL) { 
    	fprintf(stderr,"error in memory allocation\n");
    	exit(1);
    }
    for(j=0;j<JTYPES;j++){
    	if((rho[j]=(double *)calloc((size_t)4,sizeof(double)))==NULL) { 
    		fprintf(stderr,"error in memory allocation\n");
    		exit(1);
    	}
    	if((oldrho[j]=(double *)calloc((size_t)4,sizeof(double)))==NULL) { 
    		fprintf(stderr,"error in memory allocation\n");
    		exit(1);
    	}
    	if((rhodelta[j]=(long double *)calloc((size_t)4,sizeof(long double)))==NULL) { 
    		fprintf(stderr,"error in memory allocation\n");
    		exit(1);
    	}
    }
    
	foreach_l_xy(model,[&](const auto l,const auto jl,const auto j){
	//for(l=0;l<LTYPES;l++) {
		rho[J_SEX][l]=1./4;
	});		
	foreach_l_zw(model,[&](const auto l,const auto jl,const auto j){
	//for(l=0;l<LTYPES;l++) {
		rho[J_ZW][l]=1./4;
	});		
	foreach_l_para(model,[&](const auto l,const auto jl,const auto j){
	//for(l=0;l<LTYPES;l++) {
		rho[J_PARA][l]=1./2;
	});		
	calcQ(Q,e,e,e);		
	
	// Calculate conditional probabilities per site
	CondSiteProbs(contigs,model,Q,P,condsiteprob);
	CondSegProbs(contigs,model,rho,condsiteprob,condsegprob);
	//k=1;
	//Contig & current_contig = contigs[k];
	//npolysites=current_contig.snps.size();
	//for(t=0;t<npolysites;t++){
	//			for (i=0; i<6; i++) {
	//				fprintf(stdout,"%d\t",current_contig.snps[t].genotypes_by_sex[i]);
	//			}
	//	for(j=0;j<JTYPES;j++){
	//		fprintf(stdout,"%d: %Lf;\t",j,condsegprob[k][t][j]);
	//	}
	//	fprintf(stdout,"\n");
	//}
	
	
	//initial likelihood
	if (mode == SITE) {
		loglik=totalsiteloglik(contigs,pi,rho);
	}
	else {
		loglik=totalcontigloglik(contigs,pi,rho);
	}
	fprintf(stdout,"Initial values:\npi:");
	for(j=0;j<JTYPES;j++) {
		fprintf(stdout," %f",pi[j]);
	}
	if(model.xy || model.zw){
		fprintf(stdout,"; rho");
		foreach_l_xy(model,[&](const auto l,const auto jl,const auto j){
		//for(l=0;l<LTYPES;l++) {
			if(l==0 || l==2) {
				fprintf(stdout," %f (%Le)",rho[j][l],rhodelta[j][l]);
			}
		});
		foreach_l_zw(model,[&](const auto l,const auto jl,const auto j){
		//for(l=0;l<LTYPES;l++) {
			if(l==0 || l==2) {
				fprintf(stdout," %f (%Le)",rho[j][l],rhodelta[j][l]);
			}
		});
	}
	fprintf(stdout,"; e: %e",e);				
	fprintf(stdout,", log-likelihood: %Lf\n",loglik);
	
	it=0;
	plateausteps=0;
	while(plateausteps<10 && warning!=1){
		it++;
		fprintf(stdout,"Iteration %d: ",it);

		//E-step
		for (k=0;k<contigs.size();k++) {
			Contig & current_contig = contigs[k];
			if((npolysites=current_contig.snps.size())>0) {
				// para, XY and ZW detailed types
					for(t=0;t<npolysites;t++){
						if(model.xy && pi[J_SEX]>minimumvalue){
							maxl=-INFINITY;
							lmax=-1;
							foreach_l_xy(model,[&](const auto l,const auto jl,const auto j){
									//for(l=0;l<LTYPES;l++) {
									templ[l]=logl(condsiteprob[k][t][jl]);
									//templ[l]=condsiteprob[k][t][j+l];
									if(templ[l]>=maxl){
										maxl=templ[l];
										lmax=l;
									}
							});
							loghorner(4,lmax,templ,rho[J_SEX],expA[k][t][J_SEX]);
							foreach_l_xy(model,[&](const auto l,const auto jl,const auto j){
									if(isnan(expA[k][t][j][l])){
										fprintf(stderr,"NaN produced (E-step, new value for A): contig %d, site %d, type %d, subtype %d (%d): %e\n",k,current_contig.snps[t].position,j,jl,l,expA[k][t][j][l]);
										for (i=0; i<6; i++) {
											fprintf(stderr,"%d\t",current_contig.snps[0].genotypes_by_sex[i]);
										}
										fprintf(stderr,"\n");
										foreach_jl(model,[&](const auto jl){
										fprintf(stderr,"%Le\t",condsiteprob[k][t][jl]);
										});
										fprintf(stderr,"\n");
										foreach_j(model,[&](const auto j){
										fprintf(stderr,"%d %f %Le\t",j,pi[j],condsegprob[k][t][j]);
										});
										fprintf(stderr,"\n");
										warning=1;
									}
									else if(isinf(expA[k][t][j][l])){
										fprintf(stderr,"Inf produced (expA): contig %d, site %d, type %d, subtype %d (%d): %e\n",k,current_contig.snps[t].position,j,jl,l,expA[k][t][j][l]);
									}
							});
						}
						if(model.zw && pi[J_ZW]>minimumvalue){
							maxl=-INFINITY;
							lmax=-1;
							foreach_l_zw(model,[&](const auto l,const auto jl,const auto j){
									//for(l=0;l<LTYPES;l++) {
									templ[l]=logl(condsiteprob[k][t][jl]);
									//templ[l]=condsiteprob[k][t][j+l];
									if(templ[l]>=maxl){
										maxl=templ[l];
										lmax=l;
									}
							});
							loghorner(4,lmax,templ,rho[J_ZW],expA[k][t][J_ZW]);
						}
						if(model.paralogs && pi[J_PARA]>minimumvalue){
							maxl=-INFINITY;
							lmax=-1;
							foreach_l_para(model,[&](const auto l,const auto jl,const auto j){
									//for(l=0;l<LTYPES;l++) {
									templ[l]=logl(condsiteprob[k][t][jl]);
									//templ[l]=condsiteprob[k][t][j+l];
									if(templ[l]>=maxl){
										maxl=templ[l];
										lmax=l;
									}
							});
							loghorner(2,lmax,templ,rho[J_PARA],expA[k][t][J_PARA]);
						}
						//horner(LTYPES,lmax,templ,rho,expA[k][t]);				
					}
					//segregation types
					if (mode==SITE) {		//site-wise
						for(t=0;t<npolysites;t++){
							//Expectations of segregation types
							jmax=SITEprobs(JTYPES,condsegprob[k][t],temp);
							loghorner(JTYPES,jmax,temp,pi,expS[k][t]);			
							foreach_j(model,[&](const auto j){
									//for(j=0;j<JTYPES;j++){
									if(isnan(expS[k][t][j])){
										fprintf(stderr,"NaN produced (expS): contig %d, site %d, type %d: %f\n",k,current_contig.snps[t].position,j,expS[k][t][j]);
										for (i=0; i<6; i++) {
											fprintf(stderr,"%d\t",current_contig.snps[0].genotypes_by_sex[i]);
										}
										fprintf(stderr,"\n");
										foreach_jl(model,[&](const auto jl){
										fprintf(stderr,"%Le\t",condsiteprob[k][t][jl]);
										});
										fprintf(stderr,"\n");
										foreach_j(model,[&](const auto j){
										fprintf(stderr,"%d %f %Le\t",j,pi[j],condsegprob[k][t][j]);
										});
										fprintf(stderr,"\n");
										warning=1;
									}
									else if(isinf(expS[k][t][j])){
										fprintf(stderr,"Inf produced (expS): contig %d, site %d, type %d: %f\n",k,current_contig.snps[t].position,j,expS[k][t][j]);
									}
							});
							//					horner(JTYPES,jmax,temp,pi,expS[k][t]);
						}
					}
					else { // mode == CONTIG
						jmax=contigCONTIGprobs(JTYPES,npolysites,condsegprob[k],contigllik[k]);
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
				//auto do_stuff = [&exprR, k](const j) {
					foreach_j(model,[&](const auto j){
						if(isnan(expR[k][j])){
							fprintf(stderr,"NaN produced (expR): contig %d, type %d: %f\n",k,j,expR[k][j]);
							warning=1;
						}
						else if(isinf(expR[k][j])){
							fprintf(stderr,"Inf produced (expR): contig %d, type %d: %f\n",k,j,expR[k][j]);
						}
				});
				//foreach_j (mode, do_stuff);
				
/*				for(j=0;j<JTYPES;j++){
					if(isnan(expR[k][j])){
						fprintf(stderr,"NaN produced (expR): contig %d, type %d: %f\n",k,j,expR[k][j]);
						exit(1);
					}
						else if(isinf(expR[k][j])){
							fprintf(stderr,"Inf produced (expR): contig %d, type %d: %f\n",k,j,expR[k][j]);
						}
					}*/
			}
			//Expectations of true genotypes (independent of site- or contig-wise mode)
			//if(k==0){
			//	printf("\n");
			//	for (i=0; i<6; i++) {
			//		printf("%d\t",current_contig.snps[0].genotypes_by_sex[i]);
			//	}
			//	printf("\n");
			//}
			if(errormodel==ERRORS) {
				for(t=0;t<npolysites;t++){
					foreach_jl(model,[&](const auto jl){
					//for(jl=0;jl<JLTYPES;jl++) {
						for(s=0;s<2;s++){
							//if(k==0 && t==0){
							//printf("k=%d t=%d s=%d jl=%d:",k,t,s,jl);
							//	for(gp=0;gp<3;gp++){
							//		printf("\t%f",P[k][t][s][jl][gp]);
							//	}
							//	printf("\n");
							//}
							for(g=0;g<3;g++){
								maxgp=-INFINITY;
								gpmax=-1;
								sumgp=0;
								for(gp=0;gp<3;gp++){
									temp[gp]=logl((long double)P[k][t][s][jl][gp]);
									//								sumgp+=temp[gp];
									if(temp[gp]>=maxgp){
										maxgp=temp[gp];
										gpmax=gp;
									}
								}
								loghorner(3,gpmax,temp,Q[g],expTG[k][t][s][jl][g]);						
							//if(k==0 && t==0){
							//	printf("g=%d:",g);
							//	for(gp=0;gp<3;gp++){
							//		printf("\t%f",expTG[k][t][s][jl][g][gp]);
							//	}
							//	printf("\n");
							//}
									for(gp=0;gp<3;gp++){
									if(isnan(expTG[k][t][s][jl][g][gp])){
										fprintf(stderr,"NaN produced (expTG): k=%d t=%d s=%d jl=%d g=%d gp=%d: %f\n",k,t,s,jl,g,gp,expR[k][j]);
										warning=1;
									}
									else if(isinf(expTG[k][t][s][jl][g][gp])){
										fprintf(stderr,"Inf produced (expTG): k=%d t=%d s=%d jl=%d g=%d gp=%d: %f\n",k,t,s,jl,g,gp,expR[k][j]);
									}
								}
							}
						}
					});
				}
			}
			}
		}
		
//		if(mode == SITE){
//			fprintf(stdout,"\nE-step: %Lf\n",condsiteexpect(contigs,pi,rho,Q));
//		}
//		else{
//			fprintf(stdout,"\nE-step: %Lf\n",condcontigexpect(contigs,pi,rho,Q));
//		}
		
		//M-step
		//Calculate new values for pi
		if (mode==SITE) {		//site-wise
			sumpi=0;
			foreach_j(model,[&](const auto j){
			//for(j=0;j<JTYPES;j++) {
				oldpi[j]=pi[j];
				sumET=0;
				totsites=0;
				for (k=0;k<contigs.size();k++) {
					Contig & current_contig = contigs[k];
					if((npolysites=current_contig.snps.size())>0) {
						for (t=0; t<npolysites; t++){
							sumET+=expS[k][t][j];
						}
						totsites+=npolysites;
					}
				}
				pi[j]=(long double)sumET/(long double)totsites;
				sumpi+=pi[j];
			});
		}
		else { // mode == CONTIG
			foreach_j(model,[&](const auto j){
			//for(j=0;j<JTYPES;j++) {
				oldpi[j]=pi[j];
			});
			sumpi=0;
			foreach_j(model,[&](const auto j){
			//for(j=0;j<JTYPES;j++){ 
				sumET=0;
				for (k=0;k<contigs.size();k++) {
					Contig & current_contig = contigs[k];
					if((npolysites=current_contig.snps.size())>0) {
						sumET+=expR[k][j];
					}
				}
				pi[j]=(long double)sumET/(long double)(ncontigs-nnoncontigs);
				sumpi+=pi[j];
			});
		}
		if (model.xy || model.zw) {
		//Calculate new values for rho
		foreach_l_xy(model,[&](const auto l,const auto jl, const auto j){
			//for(l=0;l<LTYPES;l++){
			oldrho[j][l]=rho[j][l];
		});
		foreach_l_zw(model,[&](const auto l,const auto jl, const auto j){
			//for(l=0;l<LTYPES;l++){
			oldrho[j][l]=rho[j][l];
		});

		if(model.xy){
			sumET=0;
			weights=0;
			for (k=0;k<contigs.size();k++) {
				Contig & current_contig = contigs[k];
				if((npolysites=current_contig.snps.size())>0) {
					for (t=0; t<npolysites; t++){
						if(mode==SITE){
							sumET+=(expA[k][t][J_SEX][0]+expA[k][t][J_SEX][1])*expS[k][t][J_SEX];
							weights+=expS[k][t][J_SEX];
						}
						else {
							sumET+=(expA[k][t][J_SEX][0]+expA[k][t][J_SEX][1])*expR[k][J_SEX];
							weights+=expR[k][J_SEX];
						}
					}
				}
			}
			rho[J_SEX][0]=(long double)sumET/(long double)(2.*weights);
			rho[J_SEX][1]=rho[J_SEX][0];
			rho[J_SEX][2]=0.5-rho[J_SEX][0];
			rho[J_SEX][3]=rho[J_SEX][2];
		}
		if(model.zw){
			sumET=0;
			weights=0;
			for (k=0;k<contigs.size();k++) {
				Contig & current_contig = contigs[k];
				if((npolysites=current_contig.snps.size())>0) {
					for (t=0; t<npolysites; t++){
						if(mode==SITE){
							sumET+=(expA[k][t][J_ZW][0]+expA[k][t][J_ZW][1])*expS[k][t][J_ZW];
							weights+=expS[k][t][J_ZW];
						}
						else {
							sumET+=(expA[k][t][J_ZW][0]+expA[k][t][J_ZW][1])*expR[k][J_ZW];
							weights+=expR[k][J_ZW];
						}
					}
				}
			}
			rho[J_ZW][0]=(long double)sumET/(long double)(2.*weights);
			rho[J_ZW][1]=rho[J_ZW][0];
			rho[J_ZW][2]=0.5-rho[J_ZW][0];
			rho[J_ZW][3]=rho[J_ZW][2];
		}
		sumrho=0;
		foreach_l_xy(model,[&](const auto l,const auto jl, const auto j){
			//for(l=0;l<LTYPES;l++){
			sumrho+=rho[j][l];
		});
		foreach_l_zw(model,[&](const auto l,const auto jl, const auto j){
			//for(l=0;l<LTYPES;l++){
			sumrho+=rho[j][l];
		});
		}
		//Calculate new values for e
		if(errormodel==ERRORS){
			for(g=0;g<3;g++){
				//				printf("\n");
				for(gp=0;gp<3;gp++){
					U[g][gp]=0;
					for (k=0;k<contigs.size();k++) {
						Contig & current_contig = contigs[k];
						if((npolysites=current_contig.snps.size())>0) {
							for (t=0; t<npolysites; t++){
								for(s=0;s<2;s++) {
									temp[s]=0;
									auto lfunc = [&](const auto l,const auto jl,const auto j){
										if(mode==SITE){
											temp[s]+=(long double)expTG[k][t][s][jl][g][gp]*(long double)expS[k][t][j]*(long double)expA[k][t][j][l];
											if(isnan(temp[s])){
												fprintf(stderr,"NaN produced (M-step, new value for e): contig %d, site %d, sex %d, type %d, subtype %d (%d): %Le\n",k,current_contig.snps[t].position,s,j,jl,l,temp[s]);
												fprintf(stderr,"%e\t%e\t%e\n",expTG[k][t][s][jl][g][gp],expS[k][t][j],expA[k][t][j][l]);
											}										
										}
										else {
											temp[s]+=(long double)expTG[k][t][s][jl][g][gp]*(long double)expR[k][j]*(long double)expA[k][t][j][l];
											if(isnan(temp[s])){
												fprintf(stderr,"NaN produced (M-step, new value for e): contig %d, site %d, sex %d, type %d, subtype %d (%d): %Le\n",k,current_contig.snps[t].position,s,j,jl,l,temp[s]);
												fprintf(stderr,"%e\t%e\t%e\n",expTG[k][t][s][jl][g][gp],expR[k][j],expA[k][t][j][l]);
											}										
										}
									};
									foreach_j(model,[&](const auto j){
											//for(j=0;j<JTYPES;j++) {
											if(j== J_AUTO || j== J_HAPLOID){
												if(mode==SITE){
													temp[s]+=(long double)expTG[k][t][s][j][g][gp]*(long double)expS[k][t][j];
												}
												else {
													temp[s]+=(long double)expTG[k][t][s][j][g][gp]*(long double)expR[k][j];
												}
												if(isnan(temp[s])){
													fprintf(stderr,"NaN produced (M-step, new value for e): contig %d, site %d, sex %d, type %d: %Le\n",k,current_contig.snps[t].position,s,j,temp[s]);
												}
											}
											if(j==J_PARA){
												foreach_l_para(model,lfunc);
												//if(mode==SITE){
												//	temp[s]+=(long double)(expTG[k][t][s][JL_PARA1][g][gp]+expTG[k][t][s][JL_PARA2][g][gp])*0.5*(long double)expS[k][t][j];
												//}
												//else {
												//	temp[s]+=(long double)(expTG[k][t][s][JL_PARA1][g][gp]+expTG[k][t][s][JL_PARA2][g][gp])*0.5*(long double)expTG[k][t][s][j][g][gp]*(long double)expR[k][j];
												//}
											}
											if(j==J_HEMI){
												if(mode==SITE){
													temp[s]+=(long double)expTG[k][t][s][JL_HEMI][g][gp]*(long double)expS[k][t][j];
												}
												else {
													temp[s]+=(long double)expTG[k][t][s][JL_HEMI][g][gp]*(long double)expR[k][j];
												}
												if(isnan(temp[s])){
													fprintf(stderr,"NaN produced (M-step, new value for e): contig %d, site %d, sex %d, type %d: %Le\n",k,current_contig.snps[t].position,s,j,temp[s]);
												}
											}
											if(j==J_SEX){
												foreach_l_xy(model,lfunc);
											}
											if(j==J_ZHEMI){
												if(mode==SITE){
													temp[s]+=(long double)expTG[k][t][s][JL_ZHEMI][g][gp]*(long double)expS[k][t][j];
												}
												else {
													temp[s]+=(long double)expTG[k][t][s][JL_ZHEMI][g][gp]*(long double)expR[k][j];
												}
												if(isnan(temp[s])){
													fprintf(stderr,"NaN produced (M-step, new value for e): contig %d, site %d, sex %d, type %d: %Le\n",k,current_contig.snps[t].position,s,j,temp[s]);
												}
											}
											if(j==J_ZW){
												foreach_l_zw(model,lfunc);
											}
									});
									if(isnan(temp[s])){
										fprintf(stderr,"NaN produced (M-step, new value for e): contig %d, site %d, sex %d: %Le\n",k,current_contig.snps[t].position,s,temp[s]);
									}
									temp[s]*=(long double)current_contig.snps[t].genotypes_by_sex[g+3*s];
									if(isnan(temp[s])){
										fprintf(stderr,"NaN produced (M-step, new value for e): contig %d, site %d, sex %d: %Le\n",k,current_contig.snps[t].position,s,temp[s]);
										warning=1;
									}
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
		if(model.xy || model.zw || model.haploid || model.paralogs){
			fprintf(stdout,"pi");
			//foreach_j(model,[&](const auto j){
			for(j=0;j<JTYPES;j++) {
				pidelta[j]= pi[j]> minimumvalue ? (pi[j]-oldpi[j])/oldpi[j] : 0.0;
				if (fabsl(pidelta[j]) > fabsl(deltamax)) {
					deltamax=pidelta[j];
				}
				fprintf(stdout," %f (%Le)",pi[j],pidelta[j]);
			}
			//});	
			fprintf(stdout,"; total: %f",sumpi);
			if(model.xy || model.zw){
				fprintf(stdout,"; rho");
				auto lfunc = [&](const auto l,const auto jl,const auto j){
					rhodelta[j][l]= rho[j][l]> minimumvalue ? (rho[j][l]-oldrho[j][l])/oldrho[j][l] : 0.0;
					if (fabsl(rhodelta[j][l]) > fabsl(deltamax)) {
						deltamax=rhodelta[j][l];
					}
					if(l==0 || l==2){
						fprintf(stdout," %f (%Le)",rho[j][l],rhodelta[j][l]);
					}					
				};
				foreach_l_xy(model,lfunc);
				foreach_l_zw(model,lfunc);
				fprintf(stdout,"; total: %f",sumrho);
			}
		}
		if(errormodel==ERRORS){
			fprintf(stdout,"; e: %e",e);
		}
		
		if (fabsl(deltamax) < stop) {
			plateausteps++;
		}
		
		// Calculate conditional probabilities per site
		CondSiteProbs(contigs,model,Q,P,condsiteprob);
		CondSegProbs(contigs,model,rho,condsiteprob,condsegprob);

		//new log-likelihood
		oldloglik=loglik;
		if(mode == SITE){
			loglik=totalsiteloglik(contigs,pi,rho);
		}
		else {
			loglik=totalcontigloglik(contigs,pi,rho);
		}
		fprintf(stdout,", log-likelihood: %Lf\n",loglik);
		if(oldloglik > loglik){
			fprintf(stderr,"Warning: log-likelihood decreases by %Lf; interrupting maximization and proceeding to output\n",oldloglik-loglik);
			warning=1;
		}

//		if(mode == SITE){
//			fprintf(stdout,"M-step: %Lf\n",condsiteexpect(contigs,pi,rho,Q));
//		}
//		else{
//			fprintf(stdout,"M-step: %Lf\n",condcontigexpect(contigs,pi,rho,Q));
//		}
	}
	//End of EM algorithm. Outputting
	
	if(warning){
		fprintf(outfile,"Warning: an error occurred, probably due to numerical problems (see stderr for details). The following output should not be used for other purposes than finding the error.\n");
	}
	
	if((geom_score=(double *)malloc(sizeof(double)*JTYPES))==NULL) {
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}
	if((geom_nopi_score=(double *)malloc(sizeof(double)*JTYPES))==NULL) {
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}
	
	fprintf(outfile,"pi");
	foreach_j(model,[&](const auto j){
		//for(j=0;j<JTYPES;j++) {
		fprintf(outfile," %f ",pi[j]);
	});	
	fprintf(outfile,"; rho");                                      
	auto lfunc= [&](const auto l,const auto jl,const auto j){
		if(l==0 || l==2){
			fprintf(outfile," %f ",rho[j][l]);
		}	
	};
	foreach_l_xy(model,lfunc);
	foreach_l_zw(model,lfunc);
	fprintf(outfile,"; e0: %e",e);
	fprintf(outfile,", log-likelihood: %Lf",loglik);
	npar=model.haploid+model.paralogs+model.xy*3+model.zw*3;
	if(errormodel==ERRORS){
		npar+=1;
	}
	fprintf(outfile,", BIC (sites): %Lf",-2.*loglik+npar*log(totsites));
	fprintf(outfile,", BIC (contigs): %Lf",-2.*loglik+npar*log(ncontigs-nnoncontigs));
	sites_individuals=0;
	for (k=0;k<contigs.size();k++) {
		Contig & current_contig = contigs[k];
		npolysites=current_contig.snps.size();
		for (t=0; t<npolysites; t++){
			for (i=0; i<6; i++) {
				sites_individuals+=current_contig.snps[t].genotypes_by_sex[i];
			}			
		}
	}
	fprintf(outfile,", BIC (sites*individuals): %Lf\n",-2.*loglik+npar*log(sites_individuals));
	
	//output headers
	fprintf(outfile,"#>contig_name\tN_sites\tmean_coverage");
	foreach_j(model,[&](const auto j){
			if(j==J_AUTO){
				fprintf(outfile,"\tj_autosomal posterior_autosomal geometric_autosomal noprior_autosomal");
			}
			if(j==J_HAPLOID){
				fprintf(outfile,"\tj_haploid posterior_haploid geometric_haploid noprior_haploid");
			}
			if(j==J_PARA){
				fprintf(outfile,"\tj_paralog posterior_paralog geometric_paralog noprior_paralog");
			}
			if(j==J_HEMI){
				fprintf(outfile,"\tj_xhemizygote posterior_xhemizygote geometric_xhemizygote noprior_xhemizygote");
			}
			if(j==J_ZHEMI){
				fprintf(outfile,"\tj_zhemizygote posterior_zhemizygote geometric_zhemizygote noprior_zhemizygote");
			}
			if(j==J_SEX){
				fprintf(outfile,"\tj_xy posterior_xy geometric_xy noprior_xy");
			}
			if(j==J_ZW){
				fprintf(outfile,"\tj_zw posterior_zw geometric_zw noprior_zw");
			}
	});
	fprintf(outfile,"\tmax_posterior max_geometric max_noprior\n");
	fprintf(outfile,"#position\talleles\tN11F\tN12F\tN22F\tN11M\tN12M\tN22M\t");
	if(model.xy){
		fprintf(outfile,"fx_max fy_max\t");
		fprintf(outfile,"fx_mean fy_mean\t");
	}
	if(model.zw){
		fprintf(outfile,"fz_max fw_max\t");
		fprintf(outfile,"fz_mean fw_mean\t");		
	}
	foreach_l_xy(model,[&](const auto l,const auto jl,const auto j){
			fprintf(outfile,"subxy_%d\t",l+1);
	});
	foreach_l_zw(model,[&](const auto l,const auto jl,const auto j){
			fprintf(outfile,"subzw_%d\t",l+1);
	});
	if(mode==SITE){
		foreach_j(model,[&](const auto j){
			if(j==J_AUTO){
				fprintf(outfile,"logL_autosomal posterior_autosomal\t");
			}                    
			if(j==J_HAPLOID){    
				fprintf(outfile,"logL_haploid posterior_haploid\t");
			}                    
			if(j==J_PARA){       
				fprintf(outfile,"logL_paralog posterior_paralog\t");
			}                    
			if(j==J_HEMI){       
				fprintf(outfile,"logL_xhemizygote posterior_xhemizygote\t");
			}                    
			if(j==J_ZHEMI){      
				fprintf(outfile,"logL_zhemizygote posterior_zhemizygote\t");
			}                    
			if(j==J_SEX){        
				fprintf(outfile,"logL_xy posterior_xy\t");
			}                    
			if(j==J_ZW){         
				fprintf(outfile,"logL_zw posterior_zw\t");
			}
		});
	}
	else {
		foreach_j(model,[&](const auto j){
			if(j==J_AUTO){
				fprintf(outfile,"logL_autosomal\t");
			}                    
			if(j==J_HAPLOID){    
				fprintf(outfile,"logL_haploid\t");
			}                    
			if(j==J_PARA){       
				fprintf(outfile,"logL_paralog\t");
			}                    
			if(j==J_HEMI){       
				fprintf(outfile,"logL_xhemizygote\t");
			}                    
			if(j==J_ZHEMI){      
				fprintf(outfile,"logL_zhemizygote\t");
			}                    
			if(j==J_SEX){        
				fprintf(outfile,"logL_xy\t");
			}                    
			if(j==J_ZW){         
				fprintf(outfile,"logL_zw\t");
			}
		});
	}
	fprintf(outfile,"j_max\n");
	
	for (k=0;k<contigs.size();k++) {
		Contig & current_contig = contigs[k];
		//calculate likelihoods and posterior probabilities per contig 
			if((npolysites=current_contig.snps.size())>0) {
				if (mode == SITE ) { //calculate mean posterior
					foreach_j(model,[&](const auto j){
					//for(j=0;j<JTYPES;j++) {
						sumET=0;
						for (t=0; t<npolysites; t++){
							sumET+=expS[k][t][j];
							//						sumET+=log(condsegprob[k][t][j]/(condsegprob[k][t][bmin[k][t]]));
						}
						expR[k][j]=sumET/npolysites;
					});
				}
				//calculate geometric mean scores
				pmax=0.0;
				jmax=-1;
				foreach_j(model,[&](const auto j){
				//for(j=0;j<JTYPES;j++) {
					sumL=0;
					for (t=0; t<npolysites; t++){
						sumL+=log(condsegprob[k][t][j]);
					}
					temp[j]=sumL/npolysites;
					if(expl(temp[j])>=pmax){
							pmax=expl(temp[j]);
							jmax=j;
					}
				});
				loghorner(JTYPES,jmax,temp,pi,geom_score);
				//calculate geometric mean scores without priors
				for(j=0;j<JTYPES;j++) {
					dumpi[j]=0.;
				}
				nj=0;
				foreach_j(model,[&](const auto j){
					nj++;
				});
				foreach_j(model,[&](const auto j){
					dumpi[j]=1./nj;
				});
				
				loghorner(JTYPES,jmax,temp,dumpi,geom_nopi_score);
				

//			fprintf(outfile,">%s\t%d",&contig[k*NAME_LEN],npolysites[k]);
			fprintf(outfile,">%s\t%d",contigs[k].name.data(),npolysites);
			ni=0;
			for (t=0; t<npolysites; t++){
				for (i=0; i<6; i++) {
					ni+=current_contig.snps[t].genotypes_by_sex[i];
				}
			}
			fprintf(outfile,"\t%f",(double)ni/(double)npolysites);
				pmax=0.0;
				gmax=0.0;
				rmax=0.0;
				jmax=-1;
				lmax=-1;
				nmax=-1;
			for(j=0;j<JTYPES;j++) {
				if(geom_nopi_score[j] > rmax){
					rmax=geom_nopi_score[j];
					nmax=j;
				}				
				if(geom_score[j] > gmax){
					gmax=geom_score[j];
					lmax=j;
				}				
				if(expR[k][j] > pmax){
					pmax=expR[k][j];
					jmax=j;
				}				
			}
			foreach_j(model,[&](const auto j){
			//for(j=0;j<JTYPES;j++) {
				fprintf(outfile,"\t%d %e %e %e",j+1,expR[k][j],geom_score[j],geom_nopi_score[j]);
			});
//			fprintf(outfile,"\t%f\n",zscore(k,npolysites[k],polysite[k],rho,Q,lmax,1000,P[k],condsegprob[k]));			
			fprintf(outfile,"\t%d %d %d\n",jmax+1,lmax+1,nmax+1);
			
			for (t=0; t<npolysites; t++){
				fprintf(outfile,"%d\t",current_contig.snps[t].position);
				fprintf(outfile,"%c%c\t",int2DNA(current_contig.snps[t].alleles[0]),int2DNA(current_contig.snps[t].alleles[1]));			
				for (i=0; i<6; i++) {
					fprintf(outfile,"%d\t",current_contig.snps[t].genotypes_by_sex[i]);
				}
				//calculate estimated frequencies of allele 1 on X and Y
				//first method: fx and fy correspond to those matching with the most probable XY segregation subtype
				if(model.xy){
					pmax=0.0;
					lmax=-1;
					foreach_l_xy(model,[&](const auto l,const auto jl,const auto j){
							//for(l=0;l<LTYPES;l++){
							if (expA[k][t][j][l]>pmax) {
								lmax=l;
								pmax=expA[k][t][j][l];
							}
					});
					if(lmax==0) {
						fx=1-F[k][t][JL_SEX1];
						fy=1;
					}
					else if(lmax==1) {
						fx=F[k][t][JL_SEX2];
						fy=0;
					}
					else if(lmax==2) {
						fx=1;
						fy=1-F[k][t][JL_SEX3];
					}
					else if(lmax==3) {
						fx=0;
						fy=F[k][t][JL_SEX4];
					}
					fprintf(outfile,"%f %f\t",fx,fy);
					//second method: fx and fy correspond to those matching with the most probable XY segregation subtype
					fx=expA[k][t][J_SEX][0]*(1-F[k][t][JL_SEX1])+expA[k][t][J_SEX][1]*F[k][t][JL_SEX2]+expA[k][t][J_SEX][2];				
					fy=expA[k][t][J_SEX][0]+expA[k][t][J_SEX][2]*(1-F[k][t][JL_SEX3])+expA[k][t][J_SEX][3]*F[k][t][JL_SEX4];				
					fprintf(outfile,"%f %f\t",fx,fy);
				}
				if(model.zw){
					pmax=0.0;
					lmax=-1;
					foreach_l_zw(model,[&](const auto l,const auto jl,const auto j){
							//for(l=0;l<LTYPES;l++){
							if (expA[k][t][j][l]>pmax) {
								lmax=l;
								pmax=expA[k][t][j][l];
							}
					});
					if(lmax==0) {
						fx=1-F[k][t][JL_ZW1];
						fy=1;
					}
					else if(lmax==1) {
						fx=F[k][t][JL_ZW2];
						fy=0;
					}
					else if(lmax==2) {
						fx=1;
						fy=1-F[k][t][JL_ZW3];
					}
					else if(lmax==3) {
						fx=0;
						fy=F[k][t][JL_ZW4];
					}
					fprintf(outfile,"%f %f\t",fx,fy);
					//second method: fx and fy correspond to those matching with the most probable XY segregation subtype
					fx=expA[k][t][J_ZW][0]*(1-F[k][t][JL_ZW1])+expA[k][t][J_ZW][1]*F[k][t][JL_ZW2]+expA[k][t][J_ZW][2];				
					fy=expA[k][t][J_ZW][0]+expA[k][t][J_ZW][2]*(1-F[k][t][JL_ZW3])+expA[k][t][J_ZW][3]*F[k][t][JL_ZW4];				
					fprintf(outfile,"%f %f\t",fx,fy);
				}
				foreach_l_xy(model,[&](const auto l,const auto jl,const auto j){
						fprintf(outfile,"%e\t",expA[k][t][j][l]);
				});
				foreach_l_zw(model,[&](const auto l,const auto jl,const auto j){
						fprintf(outfile,"%e\t",expA[k][t][j][l]);
				});
				
				pmax=0.0;
				jmax=-1;
				if(mode==SITE){
					foreach_j(model,[&](const auto j){
						if (expS[k][t][j]>pmax) {
							jmax=j;
							pmax=expS[k][t][j];
						}
						fprintf(outfile,"%Lf %f\t",logl(condsegprob[k][t][j]),expS[k][t][j]);
					});
					fprintf(outfile,"%d\n",jmax+1);
				}
				else {
					foreach_j(model,[&](const auto j){
						fprintf(outfile,"%Lf\t",logl(condsegprob[k][t][j]));
						if (condsegprob[k][t][j]>pmax) {
							jmax=j;
							pmax=condsegprob[k][t][j];
						}
					});
					fprintf(outfile,"%d\n",jmax+1);
				}
			}
			//		}
		}
	}
	
	fclose(outfile);
	for(j=0;j<JTYPES;j++){
    	free(rho[j]);
    	free(oldrho[j]);
    	free(rhodelta[j]);
    }
    free(rho);
    free(oldrho);
    free(rhodelta);

	freeEM(contigs);
	free(geom_score);
	free(geom_nopi_score);

	return 0;
	
}