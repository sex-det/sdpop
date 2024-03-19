/*
	This file is part of SDpop.	

	SDpop is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SDpop is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

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
#include <cerrno>
#include <string>
#include <vector>

double ***F; //vector containing allele frequencies per segregation type
double *****P; //vector containing P matrices (per sex and segregation type)
long double ***condsiteprob,***condsegprob; //conditional probabilities per site
long double **contigllik; //conditional log likelihood per contig
double ***expS,****expA,**expR,******expTG;

Model model;

extern int NAME_LEN;

void initEM(std::vector<ContigA>& contigs) {
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
		fprintf(stdout,"error in memory allocation\n");
		exit(1);
	}
	if((F=(double ***)calloc((size_t)ncontigs,sizeof(double **)))==NULL) { 
		fprintf(stdout,"error in memory allocation\n");
		exit(1);
	}	
	if((expS=(double ***)calloc((size_t)ncontigs,sizeof(double **)))==NULL) { 
		fprintf(stdout,"error in memory allocation\n");
		exit(1);
	}
	if((expA=(double ****)calloc((size_t)ncontigs,sizeof(double ***)))==NULL) { 
		fprintf(stdout,"error in memory allocation\n");
		exit(1);
	}
	if((condsiteprob=(long double ***)calloc((size_t)ncontigs,sizeof(long double **)))==NULL) { 
		fprintf(stdout,"error in memory allocation\n");
		exit(1);
	}
	if((condsegprob=(long double ***)calloc((size_t)ncontigs,sizeof(long double **)))==NULL) { 
		fprintf(stdout,"error in memory allocation\n");
		exit(1);
	}
	if((expR=(double **)calloc((size_t)ncontigs,sizeof(double *)))==NULL) { 
		fprintf(stdout,"error in memory allocation\n");
		exit(1);
	}
	if((contigllik=(long double **)calloc((size_t)ncontigs,sizeof(long double *)))==NULL) { 
		fprintf(stdout,"error in memory allocation\n");
		exit(1);
	}
	if((expTG=(double ******)calloc((size_t)ncontigs,sizeof(double *****)))==NULL) { 
		fprintf(stdout,"error in memory allocation\n");
		exit(1);
	}
	
	
	for (k=0;k<contigs.size();k++) {
		ContigA & current_contig = contigs[k];
		if((npolysites=current_contig.varsites.size())>0) {
			if((P[k]=(double ****)calloc((size_t)npolysites,sizeof(double ***)))==NULL) { 
				fprintf(stdout,"error in memory allocation\n");
				exit(1);
			}
			if((F[k]=(double **)calloc((size_t)npolysites,sizeof(double *)))==NULL) { 
				fprintf(stdout,"error in memory allocation\n");
				exit(1);
			}
			if((expS[k]=(double **)calloc((size_t)npolysites,sizeof(double *)))==NULL) { 
				fprintf(stdout,"error in memory allocation\n");
				exit(1);
			}
			if((expA[k]=(double ***)calloc((size_t)npolysites,sizeof(double **)))==NULL) { 
				fprintf(stdout,"error in memory allocation\n");
				exit(1);
			}
			if((condsiteprob[k]=(long double **)calloc((size_t)npolysites,sizeof(long double *)))==NULL) { 
				fprintf(stdout,"error in memory allocation\n");
				exit(1);
			}
			if((condsegprob[k]=(long double **)calloc((size_t)npolysites,sizeof(long double *)))==NULL) { 
				fprintf(stdout,"error in memory allocation\n");
				exit(1);
			}
			if((expR[k]=(double *)calloc((size_t)JTYPES,sizeof(double)))==NULL) { 
				fprintf(stdout,"error in memory allocation\n");
				exit(1);
			}
			if((contigllik[k]=(long double *)calloc((size_t)JTYPES,sizeof(long double)))==NULL) { 
				fprintf(stdout,"error in memory allocation\n");
				exit(1);
			}
			if((expTG[k]=(double *****)calloc((size_t)npolysites,sizeof(double ****)))==NULL) { 
				fprintf(stdout,"error in memory allocation\n");
				exit(1);
			}
			for (t=0; t<npolysites; t++){
				if((P[k][t]= (double ***)calloc((size_t)SEXES,sizeof(double **)))==NULL) { 
					fprintf(stdout,"error in memory allocation\n");
					exit(1);
				}
				if((F[k][t]= (double *)calloc((size_t)JLTYPES,sizeof(double)))==NULL) { 
					fprintf(stdout,"error in memory allocation\n");
					exit(1);
				}
				if((expS[k][t]= (double *)calloc((size_t)JTYPES,sizeof(double)))==NULL) { 
					fprintf(stdout,"error in memory allocation\n");
					exit(1);
				}
				if((expA[k][t]= (double **)calloc((size_t)JTYPES,sizeof(double *)))==NULL) { 
					fprintf(stdout,"error in memory allocation\n");
					exit(1);
				}
				if((condsiteprob[k][t]= (long double *)calloc((size_t)JLTYPES,sizeof(long double)))==NULL) { 
					fprintf(stdout,"error in memory allocation\n");
					exit(1);
				}
				if((condsegprob[k][t]= (long double *)calloc((size_t)JTYPES,sizeof(long double)))==NULL) { 
					fprintf(stdout,"error in memory allocation\n");
					exit(1);
				}
				if((expTG[k][t]= (double ****)calloc((size_t)SEXES,sizeof(double ***)))==NULL) { 
					fprintf(stdout,"error in memory allocation\n");
					exit(1);
				}
				for(j=0;j<JTYPES;j++){
					if((expA[k][t][j]= (double *)calloc((size_t)4,sizeof(double)))==NULL) { 
						fprintf(stdout,"error in memory allocation\n");
						exit(1);
					}
				}

				for(s=0;s<SEXES;s++){
					if((P[k][t][s]= (double **)calloc((size_t)JLTYPES,sizeof(double *)))==NULL) { 
						fprintf(stdout,"error in memory allocation\n");
						exit(1);
					}
					if((expTG[k][t][s]= (double ***)calloc((size_t)JLTYPES,sizeof(double **)))==NULL) { 
						fprintf(stdout,"error in memory allocation\n");
						exit(1);
					}
					for(jl=0;jl<JLTYPES;jl++){
						if((P[k][t][s][jl]= (double *)calloc((size_t)3,sizeof(double)))==NULL) { 
							fprintf(stdout,"error in memory allocation\n");
							exit(1);
						}
						if((expTG[k][t][s][jl]= (double **)calloc((size_t)3,sizeof(double *)))==NULL) { 
							fprintf(stdout,"error in memory allocation\n");
							exit(1);
						}
						for(g=0;g<3;g++){
							if((expTG[k][t][s][jl][g]= (double *)calloc((size_t)3,sizeof(double)))==NULL) { 
								fprintf(stdout,"error in memory allocation\n");
								exit(1);
							}
						}
					}
				}
			}							
			
			for (t=0;t<npolysites;t++) {
				n11f=current_contig.varsites[t].genotypes_by_sex[N11F]; //female counts
				n12f=current_contig.varsites[t].genotypes_by_sex[N12F]; 
				n22f=current_contig.varsites[t].genotypes_by_sex[N22F]; 			
				n11m=current_contig.varsites[t].genotypes_by_sex[N11M]; //male counts
				n12m=current_contig.varsites[t].genotypes_by_sex[N12M]; 
				n22m=current_contig.varsites[t].genotypes_by_sex[N22M]; 			
				
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

void freeEM(std::vector<ContigA>& contigs) {
	int k,t,s,jl,j,g;
	int npolysites;
	
	for(k=0; k<contigs.size(); k++){
		ContigA & current_contig = contigs[k];
		npolysites=current_contig.varsites.size();
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

long double totalcontigloglik(std::vector<ContigA>& contigs, double *pi, double **rho) {
	long double loglik=0;
	long double suml,sumj,temp[JTYPES];
	int t,k,l;
	
	for (k=0;k<contigs.size();k++) {
		ContigA & current_contig = contigs[k];
		int npolysites=current_contig.varsites.size();
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

long double totalsiteloglik(std::vector<ContigA>& contigs, double *pi, double **rho) {
	long double loglik=0;
	long double sumj,sum3;
	int t,k;
	
	for (k=0;k<contigs.size();k++) {
		ContigA & current_contig = contigs[k];
		for(t=0;t<current_contig.varsites.size();t++){
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

long double condsiteexpect(std::vector<ContigA>& contigs, double *pi, double **rho, double Q[3][3])
{
	long double tempgs,tempjl,tempgp,logexp=0;
	int k,t,s,g,gp;
	
	for (k=0;k<contigs.size();k++) {
		ContigA & current_contig = contigs[k];
		for(t=0;t<current_contig.varsites.size();t++){
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
				tempgs+=current_contig.varsites[t].genotypes_by_sex[g+3*s]*tempjl;
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
long double condcontigexpect(std::vector<ContigA>& contigs, double *pi, double **rho, double Q[3][3])
{
	long double tempgs,tempjl,tempgp,logexp=0;
	int k,t,s,g,gp;
	
	for (k=0;k<contigs.size();k++) {
		ContigA & current_contig = contigs[k];
		for(t=0;t<current_contig.varsites.size();t++){
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
				tempgs+=current_contig.varsites[t].genotypes_by_sex[g+3*s]*tempjl;
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
//		if(current_contig.varsites.size()>0){
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
	int errormodel=ERRORS,plateausteps,nnoncontigs;
	double stop=0.0001; //relative difference between parameter values used to signal convergence
	double noconv=0.00001; //if a parameter value is below this value, stop evaluating its contribution to convergence
	double minimumvalue=1e-10; //if a parameter value falls below this value, stop evaluating its optimisation
	double Q[3][3],e;
	double maxl,sumgp;
	long double temp[JTYPES],templ[LTYPES],maxgp,sumL;
	double pmax,gmax,rmax,fx,fy;
	int	jmax,lmax,nmax,gpmax,jmax_nopi;
	long double oldloglik,loglik;
	double U[3][3],Usim,Udis,*expS_nopi,pmax_nopi;
	int nj,ni,npi=1,sites_individuals,npar;
	double *geom_score,*geom_nopi_score,mincov;
	int warning=0,expert=0;

	std::vector<ContigA> contigs;

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
	
	if (argc < 8 || argc > 10) {
		fprintf(stdout,"Usage: %s infile outfile errormodel heterogamety ploidy paralogs expert (min_coverage)\n",argv[0]);
		exit(1);
	}

	i=1;
	if((fp=fopen(argv[i],"r"))==NULL){
		fprintf(stdout,"error opening input file %s\n",argv[1]);
		exit(1);
	}
	i++;
	if((outfile=fopen(argv[i],"w"))==NULL){
		fprintf(stdout,"error opening output file %s\n",argv[2]);
		exit(1);
	}
	i++;

	e=0.0001;
	if(strcmp(argv[i],"n")==0 || strcmp(argv[i],"0")==0) {
		errormodel=NONE;
		e=0.;
	}
	else if (strcmp(argv[i],"e")==0 || strcmp(argv[i],"1")==0) {
		errormodel=ERRORS;
	}
	else if (strcmp(argv[i],"f")==0 || strcmp(argv[i],"2")==0) {
		errormodel=FIXED;
		i++;
		char * fin;
		errno=0;
		e = std::strtod(argv[i], &fin);
		if (*fin != '\0' || errno != 0 ) {
			fprintf(stdout,"Failed to convert argument %d (\"%s\") to double\n",i,argv[i]);
			exit(1);
		}		
	}
	else if (strcmp(argv[i],"r")==0 || strcmp(argv[i],"3")==0) {
		errormodel=READ_FROM_CNT_FILE;
	}
	else {
		fprintf(stdout,"Usage: %s infile outfile errormodel heterogamety ploidy paralogs expert (min_coverage)\n",argv[0]);
		fprintf(stdout,"Errormodel should be either \"e\" or \"1\" to estimate errors, \"f\" or \"2\" to use fixed error rates, or \"n\" or \"0\" for no errors\n");
		exit(1);	
	}

	i++;
	if(strcmp(argv[i],"n")==0 || strcmp(argv[i],"0")==0) {
		model.xy=0;
		model.zw=0;
	}
	else if (strcmp(argv[i],"x")==0 || strcmp(argv[i],"1")==0) {
		model.xy=1;
		model.zw=0;
		npi+=2;
	}
	else if (strcmp(argv[i],"z")==0 || strcmp(argv[i],"2")==0) {
		model.xy=0;
		model.zw=1;
		npi+=2;
	}
	else if (strcmp(argv[i],"b")==0 || strcmp(argv[i],"3")==0) {
		model.xy=1;
		model.zw=1;
		npi+=2;
	}
	else {
		fprintf(stdout,"Usage: %s infile outfile errormodel heterogamety ploidy paralogs expert (min_coverage)\n",argv[0]);
		fprintf(stdout,"Heterogamety should be either be \"x\" or \"1\" for XY type, \"z\" or \"2\" for ZW type, \"n\" or \"0\" for no sex chromosomes, or \"b\" or \"3\" for both.\n");
		exit(1);	
	}
	i++;

	if(strcmp(argv[i],"d")==0 || strcmp(argv[i],"0")==0) {
		model.haploid=0;
	}
	else if (strcmp(argv[i],"h")==0 || strcmp(argv[i],"1")==0) {
		model.haploid=1;
		npi++;
	}
	else {
		fprintf(stdout,"Usage: %s infile outfile errormodel heterogamety ploidy paralogs expert (min_coverage)\n",argv[0]);
		fprintf(stdout,"Ploidy should be either be \"d\" or \"0\" for diploid only, \"h\" or \"1\" for including haploid genes\n");
		exit(1);	
	}
	i++;
	if(strcmp(argv[i],"o")==0 || strcmp(argv[i],"0")==0) {
		model.paralogs=0;
	}
	else if (strcmp(argv[i],"p")==0 || strcmp(argv[i],"1")==0) {
		model.paralogs=1;
		npi++;
	}
	else {
		fprintf(stdout,"Usage: %s infile outfile errormodel heterogamety ploidy paralogs expert (min_coverage)\n",argv[0]);
		fprintf(stdout,"Paralogs should be either be \"o\" or \"0\" for orthologs only, \"p\" or \"1\" for including paralogy\n");
		exit(1);	
	}
	i++;
	if(strcmp(argv[i],"s")==0 || strcmp(argv[i],"0")==0) {
		expert=0;
	}
	else if (strcmp(argv[i],"e")==0 || strcmp(argv[i],"1")==0) {
		expert=1;
	}
	else {
		fprintf(stdout,"Usage: %s infile outfile errormodel heterogamety ploidy paralogs expert (min_coverage)\n",argv[0]);
		fprintf(stdout,"expert should be either be \"s\" or \"0\" for simple output mode, \"e\" or \"1\" for expert output mode\n");
		exit(1);	
	}

	i++;
	if (argc == i+1 ) {
		char * fin;
		errno=0;
		mincov = std::strtod(argv[i], &fin);
		if (*fin != '\0' || errno != 0 ) {
			fprintf(stdout,"Failed to convert argument %d (\"%s\") to double\n",i,argv[i]);
			exit(1);
		}
	}
	else {
		fprintf(stdout,"No minimal coverage given; using all contigs with polymorphism\n");
		mincov=0;
	}
	
	foreach_j(model,[&](const auto j){
			pi[j]=1./(double)npi;
	});
	
	if(expert){
		fprintf(stdout,"Some hardcoded or system-dependent values:\n");
		fprintf(stdout,"Minimum positive value for double: %e\n",DBL_MIN);
		fprintf(stdout,"Minimum positive value for double (GSL): %e\n",GSL_DBL_MIN);
		fprintf(stdout,"Minimum positive value for long double: %Le\n",LDBL_MIN);
	}
	fprintf(stdout,"Criterium for convergence: delta < %e\n",stop);
	fprintf(stdout,"\n");
	fprintf(stdout,"Reading...\n");
	if(errormodel==READ_FROM_CNT_FILE){
		e=-1;
		ncontigs=read_cnt_model_error(fp,NAME_LEN,model,contigs,&e);
		if(e<0){
			fprintf(stdout,"Error: reading error rate from cnt file failed. Aborting.\n");
			exit(1);
		}
	}
	else {
		ncontigs=read_cnt_model(fp,NAME_LEN,model,contigs);
	}
	fclose(fp);
	
	fprintf(stdout,"Found %d contigs\n",ncontigs);
	totsites=0;
	nnoncontigs=0;
	for (int k=0;k<contigs.size();k++) {
		ContigA & current_contig = contigs[k];
		npolysites=current_contig.varsites.size();
		totsites+=npolysites;
		if(npolysites==0) {
			nnoncontigs++;
		}
		ni=0;
		for (t=0; t<npolysites; t++){
			for (i=0; i<6; i++) {
				ni+=current_contig.varsites[t].genotypes_by_sex[i];
			}
		}
		current_contig.coverage=(double)ni/(double)npolysites;
	}
	fprintf(stdout,"...and %d polymorphic sites\n",totsites);
	std::vector<ContigA>::iterator kcontig = contigs.begin();
	while (kcontig != contigs.end()) {
		if (kcontig->coverage < mincov || kcontig->varsites.size() ==0) {
        	kcontig = contigs.erase(kcontig);
        }
        else {
        	kcontig++;
        }
	}
	fprintf(stdout,"%d contigs remaining after filtering for coverage\n",contigs.size());
	if(contigs.size()==0){
		fprintf(stdout,"No contigs left; nothing to do.\n");
		exit(0);
	}
	if(totsites==0){
		fprintf(stdout,"No polymorphic sites found; nothing to do.\n");
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
    	fprintf(stdout,"error in memory allocation\n");
    	exit(1);
    }
    if((oldrho=(double **)calloc((size_t)JTYPES,sizeof(double *)))==NULL) { 
    	fprintf(stdout,"error in memory allocation\n");
    	exit(1);
    }
    if((rhodelta=(long double **)calloc((size_t)JTYPES,sizeof(long double *)))==NULL) { 
    	fprintf(stdout,"error in memory allocation\n");
    	exit(1);
    }
    for(j=0;j<JTYPES;j++){
    	if((rho[j]=(double *)calloc((size_t)4,sizeof(double)))==NULL) { 
    		fprintf(stdout,"error in memory allocation\n");
    		exit(1);
    	}
    	if((oldrho[j]=(double *)calloc((size_t)4,sizeof(double)))==NULL) { 
    		fprintf(stdout,"error in memory allocation\n");
    		exit(1);
    	}
    	if((rhodelta[j]=(long double *)calloc((size_t)4,sizeof(long double)))==NULL) { 
    		fprintf(stdout,"error in memory allocation\n");
    		exit(1);
    	}
    }
    
	foreach_l_xy(model,[&](const auto l,const auto jl,const auto j){
		rho[J_SEX][l]=1./4;
	});		
	foreach_l_zw(model,[&](const auto l,const auto jl,const auto j){
		rho[J_ZW][l]=1./4;
	});		
	foreach_l_para(model,[&](const auto l,const auto jl,const auto j){
		rho[J_PARA][l]=1./2;
	});		
	calcQ(Q,e,e,e);		
	
	// Calculate conditional probabilities per site
	CondSiteProbs(contigs,model,Q,P,condsiteprob);
	CondSegProbs(contigs,model,rho,condsiteprob,condsegprob);	
	
	//initial likelihood
	loglik=totalsiteloglik(contigs,pi,rho);

	fprintf(stdout,"Initial values:\npi:");
	for(j=0;j<JTYPES;j++) {
		fprintf(stdout," %f",pi[j]);
	}
	if(model.xy || model.zw){
		fprintf(stdout,"; rho");
		foreach_l_xy(model,[&](const auto l,const auto jl,const auto j){
			if(l==0 || l==2) {
				fprintf(stdout," %f (%Le)",rho[j][l],rhodelta[j][l]);
			}
		});
		foreach_l_zw(model,[&](const auto l,const auto jl,const auto j){
			if(l==0 || l==2) {
				fprintf(stdout," %f (%Le)",rho[j][l],rhodelta[j][l]);
			}
		});
	}
	fprintf(stdout,"; e: %e",e);				
	fprintf(stdout,", log-likelihood: %Lf\n",loglik);
	
	it=0;
	plateausteps=0;
//	warning=1;
	while(plateausteps<10 && warning==0){
		it++;
		fprintf(stdout,"Iteration %d: ",it);

		//E-step
		for (k=0;k<contigs.size();k++) {
			ContigA & current_contig = contigs[k];
			if((npolysites=current_contig.varsites.size())>0) {
				// para, XY and ZW detailed types
					for(t=0;t<npolysites;t++){
						if(model.xy && pi[J_SEX]>minimumvalue){
							maxl=-INFINITY;
							lmax=-1;
							foreach_l_xy(model,[&](const auto l,const auto jl,const auto j){
									templ[l]=logl(condsiteprob[k][t][jl]);
									if(templ[l]>=maxl){
										maxl=templ[l];
										lmax=l;
									}
							});
							loghorner(4,lmax,templ,rho[J_SEX],expA[k][t][J_SEX]);
							foreach_l_xy(model,[&](const auto l,const auto jl,const auto j){
									if(isnan(expA[k][t][j][l])){
										fprintf(stdout,"NaN produced (E-step, new value for A): contig %d, site %d, type %d, subtype %d (%d): %e\n",k,current_contig.varsites[t].position,j,jl,l,expA[k][t][j][l]);
										for (i=0; i<6; i++) {
											fprintf(stdout,"%d\t",current_contig.varsites[0].genotypes_by_sex[i]);
										}
										fprintf(stdout,"\n");
										foreach_jl(model,[&](const auto jl){
										fprintf(stdout,"%Le\t",condsiteprob[k][t][jl]);
										});
										fprintf(stdout,"\n");
										foreach_j(model,[&](const auto j){
										fprintf(stdout,"%d %f %Le\t",j,pi[j],condsegprob[k][t][j]);
										});
										fprintf(stdout,"\n");
										warning=1;
									}
									else if(isinf(expA[k][t][j][l])){
										fprintf(stdout,"Inf produced (expA): contig %d, site %d, type %d, subtype %d (%d): %e\n",k,current_contig.varsites[t].position,j,jl,l,expA[k][t][j][l]);
									}
							});
						}
						if(model.zw && pi[J_ZW]>minimumvalue){
							maxl=-INFINITY;
							lmax=-1;
							foreach_l_zw(model,[&](const auto l,const auto jl,const auto j){
									templ[l]=logl(condsiteprob[k][t][jl]);
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
									templ[l]=logl(condsiteprob[k][t][jl]);
									if(templ[l]>=maxl){
										maxl=templ[l];
										lmax=l;
									}
							});
							loghorner(2,lmax,templ,rho[J_PARA],expA[k][t][J_PARA]);
						}
					}
					//segregation types
						for(t=0;t<npolysites;t++){
							//Expectations of segregation types
							jmax=SITEprobs(JTYPES,condsegprob[k][t],temp);
							loghorner(JTYPES,jmax,temp,pi,expS[k][t]);			
							foreach_j(model,[&](const auto j){
									if(isnan(expS[k][t][j])){
										fprintf(stdout,"NaN produced (expS): contig %d, site %d, type %d: %f\n",k,current_contig.varsites[t].position,j,expS[k][t][j]);
										for (i=0; i<6; i++) {
											fprintf(stdout,"%d\t",current_contig.varsites[0].genotypes_by_sex[i]);
										}
										fprintf(stdout,"\n");
										foreach_jl(model,[&](const auto jl){
										fprintf(stdout,"%Le\t",condsiteprob[k][t][jl]);
										});
										fprintf(stdout,"\n");
										foreach_j(model,[&](const auto j){
										fprintf(stdout,"%d %f %Le\t",j,pi[j],condsegprob[k][t][j]);
										});
										fprintf(stdout,"\n");
										warning=1;
									}
									else if(isinf(expS[k][t][j])){
										fprintf(stdout,"Inf produced (expS): contig %d, site %d, type %d: %f\n",k,current_contig.varsites[t].position,j,expS[k][t][j]);
									}
							});
						}
			//Expectations of true genotypes (independent of site- or contig-wise mode)
			if(errormodel==ERRORS) {
				for(t=0;t<npolysites;t++){
					foreach_jl(model,[&](const auto jl){
						for(s=0;s<2;s++){
							for(g=0;g<3;g++){
								maxgp=-INFINITY;
								gpmax=-1;
								sumgp=0;
								for(gp=0;gp<3;gp++){
									temp[gp]=logl((long double)P[k][t][s][jl][gp]);
									if(temp[gp]>=maxgp){
										maxgp=temp[gp];
										gpmax=gp;
									}
									sumgp+=temp[gp];
								}
								if(isnan(sumgp)){
									fprintf(stdout,"NaN produced (expTG): k=%d t=%d s=%d jl=%d g=%d gp=%d: %f %f\n",k,t,s,jl,g,gp,expR[k][j],temp[gp]);
									for(gp=0;gp<3;gp++){
										fprintf(stdout,"temp[%d]: %Le\t",gp,temp[gp]);
									}
									fprintf(stdout,"\nn =");
									for(i=0;i<6;i++){
										fprintf(stdout," %d",contigs[k].varsites[t].genotypes_by_sex[i]);
									}
									fprintf(stdout,"\n");
									warning=1;
								}
								loghorner(3,gpmax,temp,Q[g],expTG[k][t][s][jl][g]);						
								for(gp=0;gp<3;gp++){
									if(isnan(expTG[k][t][s][jl][g][gp])){
										fprintf(stdout,"NaN produced (expTG): k=%d t=%d s=%d jl=%d g=%d gp=%d: %f %f\n",k,t,s,jl,g,gp,expR[k][j],temp[gp]);
										fprintf(stdout,"Q:%e\t",Q[g][gp]);
										for(i=0;i<3;i++){
											fprintf(stdout,"temp[%d]: %Le\t",i,temp[i]);
										}
										fprintf(stdout,"\n");
										warning=1;
									}
									else if(isinf(expTG[k][t][s][jl][g][gp])){
										fprintf(stdout,"Inf produced (expTG): k=%d t=%d s=%d jl=%d g=%d gp=%d: %f\n",k,t,s,jl,g,gp,expR[k][j]);
									}
								}
							}
						}
					});
				}
			}
			}
		}
		
		//Calculate new values for pi
		sumpi=0;
		foreach_j(model,[&](const auto j){
				oldpi[j]=pi[j];
				sumET=0;
				totsites=0;
				for (k=0;k<contigs.size();k++) {
					ContigA & current_contig = contigs[k];
					if((npolysites=current_contig.varsites.size())>0) {
						for (t=0; t<npolysites; t++){
							sumET+=expS[k][t][j];
						}
						totsites+=npolysites;
					}
				}
				pi[j]=(long double)sumET/(long double)totsites;
				sumpi+=pi[j];
		});
		if (model.xy || model.zw) {
			//Calculate new values for rho
			foreach_l_xy(model,[&](const auto l,const auto jl, const auto j){
			oldrho[j][l]=rho[j][l];
		});
		foreach_l_zw(model,[&](const auto l,const auto jl, const auto j){
			oldrho[j][l]=rho[j][l];
		});

		if(model.xy){
			sumET=0;
			weights=0;
			for (k=0;k<contigs.size();k++) {
				ContigA & current_contig = contigs[k];
				if((npolysites=current_contig.varsites.size())>0) {
					for (t=0; t<npolysites; t++){
							sumET+=(expA[k][t][J_SEX][0]+expA[k][t][J_SEX][1])*expS[k][t][J_SEX];
							weights+=expS[k][t][J_SEX];
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
				ContigA & current_contig = contigs[k];
				if((npolysites=current_contig.varsites.size())>0) {
					for (t=0; t<npolysites; t++){
							sumET+=(expA[k][t][J_ZW][0]+expA[k][t][J_ZW][1])*expS[k][t][J_ZW];
							weights+=expS[k][t][J_ZW];
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
			sumrho+=rho[j][l];
		});
		foreach_l_zw(model,[&](const auto l,const auto jl, const auto j){
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
						ContigA & current_contig = contigs[k];
						if((npolysites=current_contig.varsites.size())>0) {
							for (t=0; t<npolysites; t++){
								for(s=0;s<2;s++) {
									temp[s]=0;
									auto lfunc = [&](const auto l,const auto jl,const auto j){
										temp[s]+=(long double)expTG[k][t][s][jl][g][gp]*(long double)expS[k][t][j]*(long double)expA[k][t][j][l];
										if(isnan(temp[s])){
											fprintf(stdout,"NaN produced (M-step, new value for e): contig %d, site %d, sex %d, type %d, subtype %d (%d): %Le\n",k,current_contig.varsites[t].position,s,j,jl,l,temp[s]);
											fprintf(stdout,"%e\t%e\t%e\n",expTG[k][t][s][jl][g][gp],expS[k][t][j],expA[k][t][j][l]);
											warning=1;
										}										
									};
									foreach_j(model,[&](const auto j){
											if(j== J_AUTO || j== J_HAPLOID){
												temp[s]+=(long double)expTG[k][t][s][j][g][gp]*(long double)expS[k][t][j];
												if(isnan(temp[s])){
													fprintf(stdout,"NaN produced (M-step, new value for e): contig %d, site %d, sex %d, type %d: %Le\n",k,current_contig.varsites[t].position,s,j,temp[s]);
												warning=1;
												}
											}
											if(j==J_PARA){
												foreach_l_para(model,lfunc);
											}
											if(j==J_HEMI){
												temp[s]+=(long double)expTG[k][t][s][JL_HEMI][g][gp]*(long double)expS[k][t][j];
												if(isnan(temp[s])){
													fprintf(stdout,"NaN produced (M-step, new value for e): contig %d, site %d, sex %d, type %d: %Le\n",k,current_contig.varsites[t].position,s,j,temp[s]);
												warning=1;
												}
											}
											if(j==J_SEX){
												foreach_l_xy(model,lfunc);
											}
											if(j==J_ZHEMI){
												temp[s]+=(long double)expTG[k][t][s][JL_ZHEMI][g][gp]*(long double)expS[k][t][j];
												if(isnan(temp[s])){
													fprintf(stdout,"NaN produced (M-step, new value for e): contig %d, site %d, sex %d, type %d: %Le\n",k,current_contig.varsites[t].position,s,j,temp[s]);
												warning=1;
												}
											}
											if(j==J_ZW){
												foreach_l_zw(model,lfunc);
											}
									});
									if(isnan(temp[s])){
										fprintf(stdout,"NaN produced (M-step, new value for e): contig %d, site %d, sex %d: %Le\n",k,current_contig.varsites[t].position,s,temp[s]);
												warning=1;
									}
									temp[s]*=(long double)current_contig.varsites[t].genotypes_by_sex[g+3*s];
									if(isnan(temp[s])){
										fprintf(stdout,"NaN produced (M-step, new value for e): contig %d, site %d, sex %d: %Le\n",k,current_contig.varsites[t].position,s,temp[s]);
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
			//printf("\n");
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
			//printf("Udis: %f Usim: %f\n",Udis,Usim);
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
				if (pi[j]> noconv && fabsl(pidelta[j]) > fabsl(deltamax)) {
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
					if (rho[j][l] > noconv && fabsl(rhodelta[j][l]) > fabsl(deltamax)) {
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
		loglik=totalsiteloglik(contigs,pi,rho);
		fprintf(stdout,", log-likelihood: %Lf\n",loglik);
		if(oldloglik > loglik){
			fprintf(stdout,"Warning: log-likelihood decreases by %Lf; interrupting maximization and proceeding to output\n",oldloglik-loglik);
			warning=2;
		}

	}
	//End of EM algorithm. Outputting
	
	if(warning==1){
		fprintf(outfile,"An error occurred, probably due to numerical problems (see stdout for details). The current output is provided only for debugging purposes.\n");
	}
	else if(warning==2){
		fprintf(outfile,"Warning: there have been some numerical problems. These might just be due to the limits of precision, so please re-run the program and check if the optimized values are close to the current values. If they are, it's probably safe to interpret them.\n");
	}
	
	if((geom_score=(double *)malloc(sizeof(double)*JTYPES))==NULL) {
		fprintf(stdout,"error in memory allocation\n");
		exit(1);
	}
	if((geom_nopi_score=(double *)malloc(sizeof(double)*JTYPES))==NULL) {
		fprintf(stdout,"error in memory allocation\n");
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
		ContigA & current_contig = contigs[k];
		npolysites=current_contig.varsites.size();
		for (t=0; t<npolysites; t++){
			for (i=0; i<6; i++) {
				sites_individuals+=current_contig.varsites[t].genotypes_by_sex[i];
			}			
		}
	}
	fprintf(outfile,", BIC (sites*individuals): %Lf\n",-2.*loglik+npar*log(sites_individuals));
	
	//output headers
	fprintf(outfile,"#>contig_name\tN_sites\tmean_coverage");
	if(expert){
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
		foreach_j(model,[&](const auto j){
				if(j==J_AUTO){
					fprintf(outfile,"logL_autosomal posterior_autosomal noprior_autosomal\t");
				}                    
				if(j==J_HAPLOID){    
					fprintf(outfile,"logL_haploid posterior_haploid noprior_haploid\t");
				}                    
				if(j==J_PARA){       
					fprintf(outfile,"logL_paralog posterior_paralog noprior_paralog\t");
				}                    
				if(j==J_HEMI){       
					fprintf(outfile,"logL_xhemizygote posterior_xhemizygote noprior_xhemizygote\t");
				}                    
				if(j==J_ZHEMI){      
					fprintf(outfile,"logL_zhemizygote posterior_zhemizygote noprior_zhemizygote\t");
				}                    
				if(j==J_SEX){        
					fprintf(outfile,"logL_xy posterior_xy noprior_xy\t");
				}                    
				if(j==J_ZW){         
					fprintf(outfile,"logL_zw posterior_zw noprior_zw\t");
				}
		});
		fprintf(outfile,"j_max j_max_noprior\n");
	}
	else { //simple output mode
		foreach_j(model,[&](const auto j){
				if(j==J_AUTO){
					fprintf(outfile,"\tj_autosomal noprior_autosomal");
				}
				if(j==J_HAPLOID){
					fprintf(outfile,"\tj_haploid noprior_haploid");
				}
				if(j==J_PARA){
					fprintf(outfile,"\tj_paralog noprior_paralog");
				}
				if(j==J_HEMI){
					fprintf(outfile,"\tj_xhemizygote noprior_xhemizygote");
				}
				if(j==J_ZHEMI){
					fprintf(outfile,"\tj_zhemizygote noprior_zhemizygote");
				}
				if(j==J_SEX){
					fprintf(outfile,"\tj_xy noprior_xy");
				}
				if(j==J_ZW){
					fprintf(outfile,"\tj_zw noprior_zw");
				}
		});
		fprintf(outfile,"\tmax\n");
		fprintf(outfile,"#position\talleles\tN11F\tN12F\tN22F\tN11M\tN12M\tN22M\t");
		if(model.xy){
			fprintf(outfile,"fx_mean fy_mean\t");
		}
		if(model.zw){
			fprintf(outfile,"fz_mean fw_mean\t");
		}
		foreach_l_xy(model,[&](const auto l,const auto jl,const auto j){
				fprintf(outfile,"subxy_%d\t",l+1);
		});
		foreach_l_zw(model,[&](const auto l,const auto jl,const auto j){
				fprintf(outfile,"subzw_%d\t",l+1);
		});
		foreach_j(model,[&](const auto j){
				if(j==J_AUTO){
					fprintf(outfile,"noprior_autosomal\t");
				}                    
				if(j==J_HAPLOID){    
					fprintf(outfile,"noprior_haploid\t");
				}                    
				if(j==J_PARA){       
					fprintf(outfile,"noprior_paralog\t");
				}                    
				if(j==J_HEMI){       
					fprintf(outfile,"noprior_xhemizygote\t");
				}                    
				if(j==J_ZHEMI){      
					fprintf(outfile,"noprior_zhemizygote\t");
				}                    
				if(j==J_SEX){        
					fprintf(outfile,"noprior_xy\t");
				}                    
				if(j==J_ZW){         
					fprintf(outfile,"noprior_zw\t");
				}
		});
		fprintf(outfile,"j_max_noprior\n");
	}

	if((expS_nopi= (double *)calloc((size_t)JTYPES,sizeof(double)))==NULL) { 
		fprintf(stdout,"error in memory allocation\n");
		exit(1);
	}

	for (k=0;k<contigs.size();k++) {
		ContigA & current_contig = contigs[k];
		//calculate likelihoods and posterior probabilities per contig 
		if((npolysites=current_contig.varsites.size())>0) {
			foreach_j(model,[&](const auto j){
					//for(j=0;j<JTYPES;j++) {
					sumET=0;
					for (t=0; t<npolysites; t++){
						sumET+=expS[k][t][j];
						//						sumET+=log(condsegprob[k][t][j]/(condsegprob[k][t][bmin[k][t]]));
					}
					expR[k][j]=sumET/npolysites;
			});
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
			
			
			fprintf(outfile,">%s\t%d",current_contig.name.data(),npolysites);
			fprintf(outfile,"\t%f",current_contig.coverage);
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
			if(expert){
				foreach_j(model,[&](const auto j){
						fprintf(outfile,"\t%d %e %e %e",j+1,expR[k][j],geom_score[j],geom_nopi_score[j]);
				});
				fprintf(outfile,"\t%d %d %d\n",jmax+1,lmax+1,nmax+1);
			}
			else {
				foreach_j(model,[&](const auto j){
						fprintf(outfile,"\t%d %e",j+1,geom_nopi_score[j]);
				});
				fprintf(outfile,"\t%d\n",nmax+1);
				
			}
			
			for (t=0; t<npolysites; t++){
				fprintf(outfile,"%d\t",current_contig.varsites[t].position);
				fprintf(outfile,"%s,%s\t",current_contig.varsites[t].alleles[0].data(),current_contig.varsites[t].alleles[1].data());			
				for (i=0; i<6; i++) {
					fprintf(outfile,"%d\t",current_contig.varsites[t].genotypes_by_sex[i]);
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
					if(expert){
						fx=expA[k][t][J_SEX][0]*(1-F[k][t][JL_SEX1])+expA[k][t][J_SEX][1]*F[k][t][JL_SEX2]+expA[k][t][J_SEX][2];				
						fy=expA[k][t][J_SEX][0]+expA[k][t][J_SEX][2]*(1-F[k][t][JL_SEX3])+expA[k][t][J_SEX][3]*F[k][t][JL_SEX4];				
						fprintf(outfile,"%f %f\t",fx,fy);
					}
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
					if(expert){
						fx=expA[k][t][J_ZW][0]*(1-F[k][t][JL_ZW1])+expA[k][t][J_ZW][1]*F[k][t][JL_ZW2]+expA[k][t][J_ZW][2];				
						fy=expA[k][t][J_ZW][0]+expA[k][t][J_ZW][2]*(1-F[k][t][JL_ZW3])+expA[k][t][J_ZW][3]*F[k][t][JL_ZW4];				
						fprintf(outfile,"%f %f\t",fx,fy);
					}
				}
				foreach_l_xy(model,[&](const auto l,const auto jl,const auto j){
						fprintf(outfile,"%e\t",expA[k][t][j][l]);
				});
				foreach_l_zw(model,[&](const auto l,const auto jl,const auto j){
						fprintf(outfile,"%e\t",expA[k][t][j][l]);
				});
				
				//calculate posteriors without priors
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
				jmax=SITEprobs(JTYPES,condsegprob[k][t],temp);
				loghorner(JTYPES,jmax,temp,dumpi,expS_nopi);			

				pmax=0.0;
				jmax=-1;
				pmax_nopi=0.0;
				jmax_nopi=-1;
				
				foreach_j(model,[&](const auto j){
						if (expS_nopi[j]>pmax_nopi) {
							jmax_nopi=j;
							pmax_nopi=expS_nopi[j];
						}
						if (expS[k][t][j]>pmax) {
							jmax=j;
							pmax=expS[k][t][j];
						}
						if(expert){
							fprintf(outfile,"%Lf %f %f\t",logl(condsegprob[k][t][j]),expS[k][t][j],expS_nopi[j]);
						}
						else {
							fprintf(outfile,"%f\t",expS_nopi[j]);
						}
				});
				if(expert){
					fprintf(outfile,"%d %d\n",jmax+1,jmax_nopi+1);
				}
				else {
					fprintf(outfile,"%d\n",jmax_nopi+1);
				}
			}
		}
	}
	
	fclose(outfile);
	for(j=0;j<JTYPES;j++){
    	free(rho[j]);
    	free(oldrho[j]);
    	free(rhodelta[j]);
    }
    free(expS_nopi);
    free(rho);
    free(oldrho);
    free(rhodelta);
    
	freeEM(contigs);
	free(geom_score);
	free(geom_nopi_score);

	return 0;
	
}