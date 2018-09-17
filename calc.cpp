#include "types.h"
#include "calc.h"

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

void CondSiteProbs(std::vector<Contig>& contigs, double Q[3][3], double *****P, long double ***condsiteprob)
{
	int k,t,jl,s,g,gp;
//	unsigned int *nmulti;
//	double *pmulti, *emulti;
	double tempgp;
	
	for (k=0;k<contigs.size();k++) {
		Contig & current_contig = contigs[k];
		for(t=0;t<current_contig.snps.size();t++){
			for(jl=0;jl<JLTYPES;jl++) {
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
						condsiteprob[k][t][jl]*=intpow(tempgp,current_contig.snps[t].genotypes_by_sex[g+3*s]);
					}
				}
				if(isnan(condsiteprob[k][t][jl])){
					fprintf(stdout,"NaN produced (CondSiteProbs): contig %d, site %d, type %d: %Lf\n",k,current_contig.snps[t].position,jl,condsiteprob[k][t][jl]);
					exit(1);
				}
				else if(isinf(condsiteprob[k][t][jl])){
					fprintf(stdout,"Inf produced (CondSiteProbs): contig %d, site %d, type %d: %Lf\n",k,current_contig.snps[t].position,jl,condsiteprob[k][t][jl]);
					for(s=0;s<2;s++){
						fprintf(stdout,"s=%d\n",s);
						for(g=0;g<3;g++){
							for(gp=0;gp<3;gp++){
								fprintf(stdout,"%e\t",P[k][t][s][jl][gp]);
							}
							fprintf(stdout,"\n");
						}
					}
				}
			}
		}
	}
}

void CondSegProbs(std::vector<Contig>& contigs, double *rho, long double ***condsiteprob, long double ***condsegprob)
{
	int k,t,l,j;
	
	for (k=0;k<contigs.size();k++) {
		Contig & current_contig = contigs[k];
		for(t=0;t<current_contig.snps.size();t++){
			condsegprob[k][t][J_AUTO]=condsiteprob[k][t][JL_AUTO];
			condsegprob[k][t][J_HAPLOID]=condsiteprob[k][t][JL_HAPLOID];
			condsegprob[k][t][J_PARA]=(long double)0.5*(condsiteprob[k][t][JL_PARA1]+condsiteprob[k][t][JL_PARA2]);	
			condsegprob[k][t][J_HEMI]=condsiteprob[k][t][JL_HEMI];
			condsegprob[k][t][J_SEX]=0;
			for(l=0;l<LTYPES;l++) {
				condsegprob[k][t][J_SEX]+=(long double)rho[l]*condsiteprob[k][t][JL_SEX1+l];
			}
			for(j=0;j<JTYPES;j++){
			if(isnan(condsegprob[k][t][j])){
				fprintf(stdout,"NaN produced: contig %d, site %d, type %d: %Lf\n",k,t,j,condsegprob[k][t][j]);
				exit(1);
			}
			else if(isinf(condsegprob[k][t][j])){
				fprintf(stdout,"Inf produced: contig %d, site %d, type %d: %Lf\n",k,t,j,condsegprob[k][t][j]);
			}
			}
		}
	}
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
	int i,n;
	//	long double temp[3];
	
	//check if non-zero likelihoods exist
	tot=0;
	n=0;
	for(i=0;i<length;i++){
		if(param[i]>LDBL_MIN){
			tot+=expl(logarray[i]);
			n++;
		}
	}
	if(tot<LDBL_MIN){
		for(i=0;i<length;i++){
			if(param[i]>LDBL_MIN){
				result[i]=1./(double)n;
			}
		}
		return;
	}
	
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

double zscore(int k, int npolysites, int **polysite, double *rho, double Q[3][3], int lmax, int iterations, double ****P, long double **condsegprob) {				
	int i,t,s,g,gp,type,it,l,jl;		
	int n11f,n12f,n22f,n11m,n12m,n22m,nfem,nmal,ntot;
	double pmulti[3],*emulti;
	double rtype,gtype;
	long double *ccp,csp;
	double tempgp,tempcsp,ccpmean,ccpsum,ccpsd,test,z;
	int permsite[7],poly,nnonpoly;		
	double f,PP[2][JLTYPES][3];
	
	if((ccp=(long double *)calloc((size_t)iterations,sizeof(long double)))==NULL) { 
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}
	
	//calculate test value
	test=0;
	for(t=0;t<npolysites;t++){
		test+=logl(condsegprob[t][J_SEX]);
//		for(i=0;i<7;i++){
//			fprintf(stdout,"%d\t",polysite[t][i]);
//		}
//		for(i=0;i<LTYPES;i++){
//			fprintf(stdout,"%Lf\t",condsiteprob[k][t][J_SEX+i]);
//		}
//		fprintf(stdout,"%Lf\n",condsegprob[k][t][J_SEX]);		
	}

	ccpmean=0;
	for(it=0;it<iterations;it++){
		ccp[it]=0;
		for (t=0; t<npolysites; t++){
			n11f=polysite[t][N11F]; //female counts
			n12f=polysite[t][N12F]; 
			n22f=polysite[t][N22F]; 			
			n11m=polysite[t][N11M]; //male counts
			n12m=polysite[t][N12M]; 
			n22m=polysite[t][N22M]; 			
			
			nfem=n11f+n12f+n22f;
			nmal=n11m+n12m+n22m;
			ntot=nfem+nmal;
			
			//random rho
			rtype=(double)rand()/RAND_MAX;
			if (rtype <= rho[0]) {
				type=0;
			}
			else if (rtype <= rho[0]+rho[1]) {
				type=1;
			}
			else if (rtype <= rho[0]+rho[1]+rho[2]) {
				type=2;
			}
			else {
				type=3;
			}
			
			permsite[0]=t;
			for(i=1;i<7;i++){
				permsite[i]=0;
			}
			poly=0;
			nnonpoly=0;
			while(poly==0){				
				//homozygous sex (females)
				pmulti[0]=P[t][FEMALE][JL_SEX1+type][0];
				pmulti[1]=P[t][FEMALE][JL_SEX1+type][1];
				pmulti[2]=P[t][FEMALE][JL_SEX1+type][2];
				emulti=errormult(Q, pmulti);
				for(i=0;i<nfem;i++){
					gtype=(double)rand()/RAND_MAX;
					if(gtype<emulti[0]){
						permsite[N11F]++;
					}
					else if(gtype<emulti[0]+emulti[1]){
						permsite[N12F]++;
					}
					else {
						permsite[N22F]++;
					}	
				}
				//heterozygous sex (males)
				pmulti[0]=P[t][MALE][JL_SEX1+type][0];
				pmulti[1]=P[t][MALE][JL_SEX1+type][1];
				pmulti[2]=P[t][MALE][JL_SEX1+type][2];
				emulti=errormult(Q, pmulti);
				for(i=0;i<nmal;i++){
					gtype=(double)rand()/RAND_MAX;
					if(gtype<emulti[0]){
						permsite[N11M]++;
					}
					else if(gtype<emulti[0]+emulti[1]){
						permsite[N12M]++;
					}
					else {
						permsite[N22M]++;
					}	
				}
				//test for polymorphism
				if((permsite[N12F]+permsite[N12M])>0){
					poly=1;
				}
				else if((permsite[N11F]+permsite[N11M])>0 && (permsite[N22F]+permsite[N22M])>0){
					poly=1;
				}
				else {
					nnonpoly++;
//					fprintf(stdout,"%d\t%d\t%d\t",k,t,it);
//					for(i=1;i<7;i++){
//						fprintf(stdout,"%d %d\t",polysite[t][i],permsite[i]);
//					}
//					fprintf(stdout,"\n");
					if(nnonpoly>1000){ //if we're really stuck, draw a new segregation type
						rtype=(double)rand()/RAND_MAX;
						if (rtype <= rho[0]) {
							type=0;
						}
						else if (rtype <= rho[0]+rho[1]) {
							type=1;
						}
						else if (rtype <= rho[0]+rho[1]+rho[2]) {
							type=2;
						}
						else {
							type=3;
						}						
					}
				}
			}
//			fprintf(stdout,"%d\t%d\t%f\t%f\t",k,it,rtype,gtype);
//			for(i=1;i<7;i++){
//				fprintf(stdout,"%d\t",permsite[i]);
//			}
			
			//new frequencies
			n11f=permsite[N11F]; //female counts
			n12f=permsite[N12F]; 
			n22f=permsite[N22F]; 			
			n11m=permsite[N11M]; //male counts
			n12m=permsite[N12M]; 
			n22m=permsite[N22M]; 			
			
			nfem=n11f+n12f+n22f;
			nmal=n11m+n12m+n22m;
			ntot=nfem+nmal;

			//x/y snsp ; x-polymorphism
				//allele 1 is fixed on Y; f21 is the frequency of allele 2 on X
				jl=JL_SEX1;
				f=(double)(2*n22f+n12f+n12m)/(double)(2*nfem+n11m+n12m);
				PP[FEMALE][jl][0]=(1.-f)*(1.-f);
				PP[FEMALE][jl][1]=2.*f*(1.-f);
				PP[FEMALE][jl][2]=f*f;
				PP[MALE][jl][0]=1.-f;
				PP[MALE][jl][1]=f;
				PP[MALE][jl][2]=0.;
				//allele 2 is fixed on Y; f22 is the frequency of allele 1 on X
				jl=JL_SEX2;
				f=(double)(2*n11f+n12f+n12m)/(double)(2*nfem+n12m+n22m);
				PP[FEMALE][jl][0]=f*f;
				PP[FEMALE][jl][1]=2.*f*(1.-f);
				PP[FEMALE][jl][2]=(1.-f)*(1.-f);
				PP[MALE][jl][0]=0.;
				PP[MALE][jl][1]=f;
				PP[MALE][jl][2]=1.-f;
				//x/y snsp ; y-polymorphism. f3 is the frequency, on Y, of the allele that is not fixed on X
				//case 1: allele 1 is fixed on X; f3 is the frequency of allele 2 on Y
				jl=JL_SEX3;
				//f=(double)(n12m)/(double)(n11m+n12m);
				f=(double)(n12m)/(double)(nmal);
				PP[FEMALE][jl][0]=1.;
				PP[FEMALE][jl][1]=0.;
				PP[FEMALE][jl][2]=0.;
				PP[MALE][jl][0]=1.-f;
				PP[MALE][jl][1]=f;
				PP[MALE][jl][2]=0.;
				//case 2: allele 2 is fixed on X; f3 is the frequency of allele 1 on Y
				jl=JL_SEX4;
				//f=(double)(n12m)/(double)(n22m+n12m);
				f=(double)(n12m)/(double)(nmal);
				PP[FEMALE][jl][0]=0.;
				PP[FEMALE][jl][1]=0.;
				PP[FEMALE][jl][2]=1.;
				PP[MALE][jl][0]=0.;
				PP[MALE][jl][1]=f;
				PP[MALE][jl][2]=1.-f;

			//conditional site probabilities
			csp=0;
			for(l=0;l<LTYPES;l++) {
				tempcsp=log(rho[l]);
				for(s=0;s<2;s++){
					for(g=0;g<3;g++){
						tempgp=0.;
						for(gp=0;gp<3;gp++){
							tempgp+=PP[s][JL_SEX1+l][gp]*Q[g][gp];
						}
						tempcsp+=log(tempgp)*permsite[1+3*s+g];
					}
				}
				csp+=expl((long double)tempcsp);
			}
			//			fprintf(stdout,"%f\n",log(csp));
			if(!isfinite(logl(csp))){
				fprintf(stdout,"%d\t%d\n",k,t);
			}
			ccp[it]+=logl(csp);		
		}
//		if(isfinite(ccp[it])){
			ccpmean+=ccp[it];
//		}
//		else {
//			it--;
//		}
	}
	
	//calculate mean ccp and standard deviation
	ccpmean/=iterations;
	ccpsum=0;
	for(it=0;it<iterations;it++){
		ccpsum+=intpow((double)(ccp[it]-ccpmean),2);
	}
	ccpsd=sqrt(ccpsum/iterations);
	//compare with log real value
	z=(test-ccpmean)/ccpsd;
	fprintf(stdout,"%f\t%f\t%f\t%f\n",test,ccpmean,ccpsd,z);
	
	free(ccp);
	
	return z;
}
