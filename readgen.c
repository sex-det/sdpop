//#include <R.h>
//#include <Rinternals.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_machine.h>

#define FEMALE 1
#define MALE 2
#define NTYPES 4 //number of site segregation types
#define NRTYPES 3 //number of contig segregation types
#define MINLOGL -10000 //log-likelihood corresponding to likelihood=0 ; defined to avoid -Inf
#define CONTIG 1
#define SITE 2

double* vectorize_d(double x0, double x1, double x2) {
	double vector[3];
	vector[0]=x0;
	vector[1]=x1;
	vector[2]=x2;
	return vector;
}

unsigned int* vectorize_ui(unsigned int x0, unsigned int x1, unsigned int x2) {
	unsigned int vector[3];
	vector[0]=x0;
	vector[1]=x1;
	vector[2]=x2;
	return vector;
}

double* errormult(double ematrix[3][3], double x[3]) {
	double vector[3];
	int i,j;
	
	for (i=0;i<3;i++) {
		vector[i]=0.;
		for (j=0;j<3;j++) {
			vector[i]+=x[j]*ematrix[i][j];
		}
	}
	return vector;
}

double log2exp(double logvalue) 
{
	double expvalue;
//	if (fabs(logvalue) <= MINLOGL+GSL_DBL_MIN) { //check whether the value is 0, depending on machine precision
//		expvalue=0;
//	}
//	else {
		expvalue=exp(logvalue);
//	}
	return expvalue;
}

double exp2log(double expvalue) 
{
	double logvalue;
//	if (fabs(expvalue) <= GSL_DBL_MIN) { //check whether the value is 0, depending on machine precision
//		logvalue=MINLOGL;
//	}
//	else {
		logvalue=log(expvalue);
//	}
	return logvalue;
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

double ***initEM(int npolysites, int **polysite) {
	int n11f,n12f,n22f,n11m,n12m,n22m,nfem,nmal,ntot;
	int i,j,type;
	double f0,f1,f21,f22,f3,ProbS1,ProbS2;
	double lPS[NTYPES],PS[NTYPES];
	double pmulti3[3];
	unsigned int nmulti3[3];
	double ***siteprobs;
	double PStot;
	double epsilon=0.05;
	
	//siteprobs : a npolysites x 3 x 4 matrix, containing frequencies, log-likelihoods and likelihoods per segregation type
	
	if (npolysites > 0) {
	if((siteprobs=(double ***)calloc((size_t)npolysites,sizeof(double **)))==NULL) { 
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}
	for (i=0; i<npolysites; i++){
		if((siteprobs[i]= (double **)calloc((size_t)3,sizeof(double *)))==NULL) { 
			fprintf(stderr,"error in memory allocation\n");
			exit(1);
		}
		for(j=0;j<3;j++){
			if((siteprobs[i][j]= (double *)calloc((size_t)NTYPES,sizeof(double)))==NULL) { 
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
			}
		}
	}

	
	for (i=0;i<npolysites;i++) {
		n11f=polysite[i][1]; //female counts
		n12f=polysite[i][2]; 
		n22f=polysite[i][3]; 			
		n11m=polysite[i][4]; //male counts
		n12m=polysite[i][5]; 
		n22m=polysite[i][6]; 			
		
		nfem=n11f+n12f+n22f;
		nmal=n11m+n12m+n22m;
		ntot=nfem+nmal;
						
		//Conditional probabilities per site
		//using gsl function double gsl_ran_multinomial_lnpdf (size_t K, const double p[], const unsigned int n[]) 
/*		
//      ****without errors****		
		//autosomal snps. f0 is the frequency of allele 1 (symmetry)
		f0=(double)(2*n11f+2*n11m+n12f+n12m)/(double)(2*ntot);
//		f0=0.95;
		nmulti3[0]=n11f;
		nmulti3[1]=n12f;
		nmulti3[2]=n22f;
		pmulti3[0]=f0*f0;
		pmulti3[1]=2.*f0*(1.-f0);
		pmulti3[2]=(1.-f0)*(1.-f0);
		PS[0]=gsl_ran_multinomial_pdf(3,pmulti3,nmulti3);
//		lPS[0]=gsl_ran_multinomial_lnpdf(3,pmulti3,nmulti3);
		nmulti3[0]=n11m;
		nmulti3[1]=n12m;
		nmulti3[2]=n22m;
		PS[0]*=gsl_ran_multinomial_pdf(3,pmulti3,nmulti3);
//		lPS[0]+=gsl_ran_multinomial_lnpdf(3,pmulti3,nmulti3);
		lPS[0]=exp2log(PS[0]);
		//x-hemizygous snps. f1 is the frequency of allele 1 (symmetry)
		f1=(double)(2*n11f+n12f+n11m)/(double)(2*nfem+nmal);
//		f1=0.95;
		PS[1]=0;
		if (n12m == 0) {
			nmulti3[0]=n11f;
			nmulti3[1]=n12f;
			nmulti3[2]=n22f;
			pmulti3[0]=f1*f1;
			pmulti3[1]=2.*f1*(1.-f1);
			pmulti3[2]=(1.-f1)*(1.-f1);
			PS[1]=gsl_ran_multinomial_pdf(3,pmulti3,nmulti3);
			PS[1]*=gsl_ran_binomial_pdf(n11m,f1,n11m+n22m);
		}
		lPS[1]=exp2log(PS[1]);
		//x/y snsp ; x-polymorphism
		ProbS1=ProbS2=0;
		//allele 1 is fixed on Y; f21 is the frequency of allele 2 on X
		f21=(double)(2*n22f+n12f+n12m)/(double)(2*nfem+nmal);
//		f21=0.95;
		if (n22m == 0 && n22f+n12f+n12m != 0 && n12f+n11f+n11m != 0) {
			nmulti3[0]=n11f;
			nmulti3[1]=n12f;
			nmulti3[2]=n22f;
			pmulti3[0]=(1.-f21)*(1.-f21);
			pmulti3[1]=2.*f21*(1.-f21);
			pmulti3[2]=f21*f21;
			ProbS1=gsl_ran_multinomial_pdf(3,pmulti3,nmulti3);
			ProbS1*=gsl_ran_binomial_pdf(n11m,1.-f21,n11m+n12m);
		}
		else {
			ProbS1=-100;
		}
		//allele 2 is fixed on Y; f22 is the frequency of allele 1 on X
		f22=(double)(2*n11f+n12f+n12m)/(double)(2*nfem+nmal);
//		f22=0.95;
		if (n11m == 0 && n11f+n12f+n12m != 0 && n12f+n22f+n22m != 0) {
			nmulti3[0]=n11f;
			nmulti3[1]=n12f;
			nmulti3[2]=n22f;
			pmulti3[0]=f22*f22;
			pmulti3[1]=2.*f22*(1.-f22);
			pmulti3[2]=(1.-f22)*(1.-f22);
			ProbS2=gsl_ran_multinomial_pdf(3,pmulti3,nmulti3);
			ProbS2*=gsl_ran_binomial_pdf(n12m,f22,n22m+n12m);
		}
		else {
			ProbS2=-100;
		}
		if (ProbS1<-1 && ProbS2 > 0) {
			PS[2]=0.5*ProbS2;
			siteprobs[i][0][2]=f22;
		}
		else if (ProbS2 < -1 && ProbS1 > 0) {
			PS[2]=0.5*ProbS1;
			siteprobs[i][0][2]=f21;
		}
		else if ( ProbS1 > 0 && ProbS2 > 0) {
			PS[2]=0.5*(ProbS1+ProbS2);
			siteprobs[i][0][2]=-1;
		}
		else {
			PS[2]=0;
			siteprobs[i][0][2]=-2;
		}
		lPS[2]=exp2log(PS[2]);
		//x/y snsp ; y-polymorphism. f3 is the frequency, on Y, of the allele that is not fixed on X
		f3=(nmal==0) ? 0 : (double)(n12m)/(double)(nmal);
//		f3=0.95;
		ProbS1=ProbS2=0;
		//case 1: allele 1 is fixed on X; f3 is the frequency of allele 2 on Y
		if (n12f+n22f+n22m == 0) {
			ProbS1=gsl_ran_binomial_pdf(n11m,1.-f3,n11m+n12m);
		}
		//case 2: allele 2 is fixed on X; f3 is the frequency of allele 1 on Y
		if (n12f+n11f+n11m == 0) {
			ProbS2=gsl_ran_binomial_pdf(n12m,f3,n22m+n12m);
		}
		PS[3]=0.5*(ProbS1+ProbS2);
		lPS[3]=exp2log(PS[3]);
//		****END without errors****
*/
/*
//      ****simple errors****		
		//autosomal snps. f0 is the frequency of allele 1 (symmetry)
		f0=(double)(2*n11f+2*n11m+n12f+n12m)/(double)(2*ntot);
//		f0=0.95;
		nmulti3[0]=n11f;
		nmulti3[1]=n12f;
		nmulti3[2]=n22f;
		pmulti3[0]=f0*f0;
		pmulti3[1]=2.*f0*(1.-f0);
		pmulti3[2]=(1.-f0)*(1.-f0);
		PS[0]=gsl_ran_multinomial_pdf(3,pmulti3,nmulti3);
//		lPS[0]=gsl_ran_multinomial_lnpdf(3,pmulti3,nmulti3);
		nmulti3[0]=n11m;
		nmulti3[1]=n12m;
		nmulti3[2]=n22m;
		PS[0]*=gsl_ran_multinomial_pdf(3,pmulti3,nmulti3);
//		lPS[0]+=gsl_ran_multinomial_lnpdf(3,pmulti3,nmulti3);
		lPS[0]=exp2log(PS[0]);
		//x-hemizygous snps. f1 is the frequency of allele 1 (symmetry)
		f1=(double)(2*n11f+n12f+n11m)/(double)(2*nfem+nmal);
//		f1=0.95;
		PS[1]=0;
			nmulti3[0]=n11f;
			nmulti3[1]=n12f;
			nmulti3[2]=n22f;
			pmulti3[0]=f1*f1;
			pmulti3[1]=2.*f1*(1.-f1);
			pmulti3[2]=(1.-f1)*(1.-f1);
			PS[1]=gsl_ran_multinomial_pdf(3,pmulti3,nmulti3);
			nmulti3[0]=n11m;
			nmulti3[1]=n12m;
			nmulti3[2]=n22m;
			pmulti3[0]=f1;
			pmulti3[1]=epsilon;
			pmulti3[2]=(1.-f1);
			PS[1]*=gsl_ran_multinomial_pdf(3,pmulti3,nmulti3);
		lPS[1]=exp2log(PS[1]);
		//x/y snsp ; x-polymorphism
		ProbS1=ProbS2=0;
		//allele 1 is fixed on Y; f21 is the frequency of allele 2 on X
		f21=(double)(2*n22f+n12f+n12m)/(double)(2*nfem+nmal);
//		f21=0.95;
		if (n22f+n12f+n12m != 0 && n12f+n11f+n11m != 0) {
			nmulti3[0]=n11f;
			nmulti3[1]=n12f;
			nmulti3[2]=n22f;
			pmulti3[0]=(1.-f21)*(1.-f21);
			pmulti3[1]=2.*f21*(1.-f21);
			pmulti3[2]=f21*f21;
			ProbS1=gsl_ran_multinomial_pdf(3,pmulti3,nmulti3);
			nmulti3[0]=n11m;
			nmulti3[1]=n12m;
			nmulti3[2]=n22m;
			pmulti3[0]=(1.-f21);
			pmulti3[1]=f21;
			pmulti3[2]=epsilon;
			ProbS1*=gsl_ran_multinomial_pdf(3,pmulti3,nmulti3);
		}
		else {
			ProbS1=-100;
		}
		//allele 2 is fixed on Y; f22 is the frequency of allele 1 on X
		f22=(double)(2*n11f+n12f+n12m)/(double)(2*nfem+nmal);
//		f22=0.95;
		if (n11f+n12f+n12m != 0 && n12f+n22f+n22m != 0) {
			nmulti3[0]=n11f;
			nmulti3[1]=n12f;
			nmulti3[2]=n22f;
			pmulti3[0]=f22*f22;
			pmulti3[1]=2.*f22*(1.-f22);
			pmulti3[2]=(1.-f22)*(1.-f22);
			ProbS2=gsl_ran_multinomial_pdf(3,pmulti3,nmulti3);
			nmulti3[0]=n11m;
			nmulti3[1]=n12m;
			nmulti3[2]=n22m;
			pmulti3[0]=epsilon;
			pmulti3[1]=f22;
			pmulti3[2]=(1.-f22);
			ProbS2*=gsl_ran_multinomial_pdf(3,pmulti3,nmulti3);
		}
		else {
			ProbS2=-100;
		}
		if (ProbS1<-1 && ProbS2 > 0) {
			PS[2]=0.5*ProbS2;
			siteprobs[i][0][2]=f22;
		}
		else if (ProbS2 < -1 && ProbS1 > 0) {
			PS[2]=0.5*ProbS1;
			siteprobs[i][0][2]=f21;
		}
		else if ( ProbS1 > 0 && ProbS2 > 0) {
			PS[2]=0.5*(ProbS1+ProbS2);
			siteprobs[i][0][2]=-1;
		}
		else {
			PS[2]=0;
			siteprobs[i][0][2]=-2;
		}
		lPS[2]=exp2log(PS[2]);
		//x/y snsp ; y-polymorphism. f3 is the frequency, on Y, of the allele that is not fixed on X
		f3=(nmal==0) ? 0 : (double)(n12m)/(double)(nmal);
//		f3=0.95;
		ProbS1=ProbS2=0;
		//case 1: allele 1 is fixed on X; f3 is the frequency of allele 2 on Y
			nmulti3[0]=n11f;
			nmulti3[1]=n12f;
			nmulti3[2]=n22f;
			pmulti3[0]=1.-2.*epsilon;
			pmulti3[1]=epsilon;
			pmulti3[2]=epsilon;
			ProbS1=gsl_ran_multinomial_pdf(3,pmulti3,nmulti3);
			nmulti3[0]=n11m;
			nmulti3[1]=n12m;
			nmulti3[2]=n22m;
			pmulti3[0]=1.-f3;
			pmulti3[1]=f3;
			pmulti3[2]=epsilon;
			ProbS1*=gsl_ran_multinomial_pdf(3,pmulti3,nmulti3);
		//case 2: allele 2 is fixed on X; f3 is the frequency of allele 1 on Y
			nmulti3[0]=n11f;
			nmulti3[1]=n12f;
			nmulti3[2]=n22f;
			pmulti3[0]=epsilon;
			pmulti3[1]=epsilon;
			pmulti3[2]=1.-2.*epsilon;
			ProbS2=gsl_ran_multinomial_pdf(3,pmulti3,nmulti3);
			nmulti3[0]=n11m;
			nmulti3[1]=n12m;
			nmulti3[2]=n22m;
			pmulti3[0]=epsilon;
			pmulti3[1]=f3;
			pmulti3[2]=1.-f3;
			ProbS2*=gsl_ran_multinomial_pdf(3,pmulti3,nmulti3);
		PS[3]=0.5*(ProbS1+ProbS2);
		lPS[3]=exp2log(PS[3]);
//		****END simple errors****
*/		
//      ****with errors****		

		errorm0[0]=vectorize_d(1.-2.*e-e*e,e,e*e);
		errorm0[1]=vectorize_d(2.*e,1-2.*e,2.*e);
		errorm0[2]=vectorize_d(e*e,e,1.-2.*e-e*e);
		errorm1[0]=vectorize_d(1.-e,0.,e);
		errorm1[1]=vectorize_d(0.,0.,0.);
		errorm1[2]=vectorize_d(e,0.,1.-e);
		//autosomal snps. f0 is the frequency of allele 1 (symmetry)
		f0=(double)(2*n11f+2*n11m+n12f+n12m)/(double)(2*ntot);
//		f0=0.95;
		nmulti3=vectorize_ui(n11f,n12f,n22f);
		pmulti3=vectorize_d(f0*f0, 2.*f0*(1.-f0), (1.-f0)*(1.-f0));
		emulti3=errormult(errorm0, pmulti3);
		PS[0]=gsl_ran_multinomial_pdf(3,emulti3,nmulti3);
//		lPS[0]=gsl_ran_multinomial_lnpdf(3,pmulti3,nmulti3);
		nmulti3=vectorize_ui(n11m,n12m,n22m);
		PS[0]*=gsl_ran_multinomial_pdf(3,emulti3,nmulti3);
//		lPS[0]+=gsl_ran_multinomial_lnpdf(3,pmulti3,nmulti3);
		lPS[0]=exp2log(PS[0]);
		//x-hemizygous snps. f1 is the frequency of allele 1 (symmetry)
		f1=(double)(2*n11f+n12f+n11m)/(double)(2*nfem+n11m+n22m);
//		f1=0.95;
		PS[1]=0;
		nmulti3=vectorize_ui(n11f,n12f,n22f);
		pmulti3=vectorize_d(f1*f1, 2.*f1*(1.-f1), (1.-f1)*(1.-f1));
		emulti3=errormult(errorm0, pmulti3);
		PS[1]=gsl_ran_multinomial_pdf(3,emulti3,nmulti3);
		nmulti3=vectorize_ui(n11m,n12m,n22m);
		pmulti3=vectorize_d(f1,0.,(1.-f1));
		emulti3=errormult(errorm1, pmulti3);
		PS[1]*=gsl_ran_multinomial_pdf(3,emulti3,nmulti3);
		lPS[1]=exp2log(PS[1]);
		//x/y snsp ; x-polymorphism
		ProbS1=ProbS2=0;
		//allele 1 is fixed on Y; f21 is the frequency of allele 2 on X
		f21=(double)(2*n22f+n12f+n12m)/(double)(2*nfem+n11m+n12m);
//		f21=0.95;
//		if (n22f+n12f+n12m != 0 && n12f+n11f+n11m != 0) {
		nmulti3=vectorize_ui(n11f,n12f,n22f);
		pmulti3=vectorize_d((1.-f21)*(1.-f21), 2.*f21*(1.-f21), f21*f21);
		emulti3=errormult(errorm0, pmulti3);
		ProbS1=gsl_ran_multinomial_pdf(3,emulti3,nmulti3);
		nmulti3=vectorize_ui(n11m,n12m,n22m);
		pmulti3=vectorize_d(1.-f21,f21,0.);
		emulti3=errormult(errorm0, pmulti3);
		ProbS1*=gsl_ran_multinomial_pdf(3,emulti3,nmulti3);
//		}
		//allele 2 is fixed on Y; f22 is the frequency of allele 1 on X
		f22=(double)(2*n11f+n12f+n12m)/(double)(2*nfem+n12m+n22m);
//		if (n11f+n12f+n12m != 0 && n12f+n22f+n22m != 0) {
		nmulti3=vectorize_ui(n11f,n12f,n22f);
		pmulti3=vectorize_d(f22*f22, 2.*f22*(1.-f22), (1.-f22)*(1.-f22));
		emulti3=errormult(errorm0, pmulti3);
		ProbS2=gsl_ran_multinomial_pdf(3,emulti3,nmulti3);
		nmulti3=vectorize_ui(n11m,n12m,n22m);
		pmulti3=vectorize_d(0.,f22,1.-f22);
		emulti3=errormult(errorm0, pmulti3);
		ProbS2*=gsl_ran_multinomial_pdf(3,emulti3,nmulti3);
//		}
		PS[2]=0.5*(ProbS1+ProbS2);
		siteprobs[i][0][2]=-1;
		lPS[2]=exp2log(PS[2]);
		//x/y snsp ; y-polymorphism. f3 is the frequency, on Y, of the allele that is not fixed on X
		ProbS1=ProbS2=0;
		//case 1: allele 1 is fixed on X; f3 is the frequency of allele 2 on Y
		f3=(double)(n12m)/(double)(n11m+n12m);
		nmulti3=vectorize_ui(n11f,n12f,n22f);
		pmulti3=vectorize_d(1., 0., 0.);
		emulti3=errormult(errorm0, pmulti3);
		ProbS1=gsl_ran_multinomial_pdf(3,emulti3,nmulti3);
		nmulti3=vectorize_ui(n11m,n12m,n22m);
		pmulti3=vectorize_d(1.-f3,f3,0.);
		emulti3=errormult(errorm0, pmulti3);
		ProbS1*=gsl_ran_multinomial_pdf(3,emulti3,nmulti3);
		//case 2: allele 2 is fixed on X; f3 is the frequency of allele 1 on Y
		f3=(double)(n12m)/(double)(n22m+n12m);
		nmulti3=vectorize_ui(n11f,n12f,n22f);
		pmulti3=vectorize_d(0., 0., 1.);
		emulti3=errormult(errorm0, pmulti3);
		ProbS2=gsl_ran_multinomial_pdf(3,emulti3,nmulti3);
		nmulti3=vectorize_ui(n11m,n12m,n22m);
		pmulti3=vectorize_d(0.,f3,1.-f3);
		emulti3=errormult(errorm0, pmulti3);
		ProbS2*=gsl_ran_multinomial_pdf(3,emulti3,nmulti3);
		PS[3]=0.5*(ProbS1+ProbS2);
		lPS[3]=exp2log(PS[3]);
//      ****END with errors****		

		PStot=0.0;
		siteprobs[i][0][0]=f0;
		siteprobs[i][0][1]=f1;
//		siteprobs[i][0][2]=f21;
		siteprobs[i][0][3]=f3;
		for(type=0;type<NTYPES;type++) {
			//			siteprobs[i][0][type]=f[type];
//			if (isnan(PS[type])) {
//				fprintf(stderr,"NaN observed at position %d, type %d\n",i,type);
//			}
			siteprobs[i][1][type]=lPS[type];
			siteprobs[i][2][type]=PS[type];
			PStot+=PS[type];
		}
		if (PStot<=GSL_DBL_MIN) {
			fprintf(stderr,"O total probability at position %d\n",i);
			for(type=0;type<NTYPES;type++) {
				fprintf(stderr,"%e (%f) ",PS[type],lPS[type]);
			}
			fprintf(stderr,"\n");
			fprintf(stderr,"%f\t%f\t%f\t%f\t%f\n",f0,f1,f21,f22,f3);
			exit(1);
		}
	}		
	return siteprobs;
	}
	else {
		return NULL;
	}
}

void freeEM(int npolysites, double ***siteprobs) {
	int i,j;
		
	for (i=0; i<npolysites; i++){
		for(j=0;j<3;j++){
			free(siteprobs[i][j]);
		}
		free(siteprobs[i]);
	}
	free(siteprobs);
	
}

long double totalsiteloglik(int ncontigs, int *npolysites, long double *pi, double ****siteprobs, int *piorder) {
	long double loglik=0;
	long double pnppi,pmax;
	int t,j,k,i,jmax;
	for (k=0;k<ncontigs;k++) {
		for(t=0;t<npolysites[k];t++){
			pnppi=0;
			// Horner decomposition (could be faster by storing jmax, which doesn't change throughout optimisation)
			pmax=0.0;
			jmax=-1;
			for(i=0;i<NTYPES;i++) {
				j=piorder[i];
				if (siteprobs[k][t][2][j]>pmax) {
					jmax=j;
					pmax=siteprobs[k][t][2][j];
				}
			}
			pnppi=pi[jmax];
			for(i=0;i<NTYPES;i++) {
				j=piorder[i];
				if (j!=jmax) {
					pnppi+=pi[j]*(long double)siteprobs[k][t][2][j]/pmax;
				}
			}
			pnppi*=(long double)pmax;
			//			for(i=0;i<NTYPES;i++) {
			//				j=piorder[i];
			//				pnppi+=siteprobs[k][t][2][j]*pi[j];
			//			}
			loglik+=logl(pnppi);
		}
	}
	return loglik;
}
long double totalcontigloglik(int ncontigs, int *npolysites, long double *pi, double ****siteprobs) {
	long double loglik=0;
	long double contiglik,logG[NRTYPES],logGmax;
	int t,k,l,lmax;

	for (k=0;k<ncontigs;k++) {
		if (npolysites[k]>0) {
			for(l=0;l<2;l++) {
				logG[l]=0.;
				for(t=0;t<npolysites[k];t++){
					logG[l]+=siteprobs[k][t][1][l];
				}
			}
			logG[2]=0.;
			for(t=0;t<npolysites[k];t++){
				logG[2]+=log(pi[3]*siteprobs[k][t][2][2]+(1.-pi[3])*siteprobs[k][t][2][3]);
			}
			logGmax=log(0);
			for(l=0;l<NRTYPES;l++){
				if(logG[l]>logGmax){
					logGmax=logG[l];
					lmax=l;
				}
			}
			contiglik=pi[lmax];
			for(l=0;l<NRTYPES;l++) {
				if(l!=lmax){
					contiglik+=exp(log(pi[l])+logG[l]-logG[lmax]);
				}
			}
			loglik+=logG[lmax]+log(contiglik);
		}
	}
	return loglik;
}
			
int **polyfilter(int npos, int nind, char *genotypes, int *sex, int *npolysites) {//number of positions and individuals,
	// genotype matrix, pointer to the number of polymorphic sites
	// function treats one contig
	char nucvec[4] = {0};
	int i,j,randbit;
	int div,n11f,n12f,n22f,n11m,n12m,n22m,nA,nT,nG,nC,nN;
	int **polysite=NULL;
	char pos1,pos2,nuc1,nuc2;
	int n;
		
	n=0;
	for (j=0; j<npos; j++){
		//		printf("%d\t",j+1);
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
		//		if (div<=1) { //no polymorphism
		//		}
		if (div == 2) { //simple polymorphism
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
				else { //if MALE
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
			//			else if (div > 2) { //more than two alleles observed : what to do ? 
			//			}
		}
		//		for (i=0; i<nind; i++){
		//			printf("%c%c\t",genotypes[j*nind*2+2*i],genotypes[j*nind*2+2*i+1]);
		//		}
		//		printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\n",div,n11f,n12f,n22f,n11m,n12m,n22m);
		if (div == 2 ) {
			if (n==0) {
				if((polysite=(int **)malloc(sizeof(int *)))==NULL){
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}
				if((polysite[n]=(int *)malloc(sizeof(int)*7))==NULL){
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}
			}
			else {
				if((polysite=realloc(polysite,sizeof(int *)*(n+1)))==NULL){
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}
				if((polysite[n]=(int *)malloc(sizeof(int)*7))==NULL){
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}
			}
			polysite[n][0]=j+1; //position
			polysite[n][1]=n11f; //female counts
			polysite[n][2]=n12f; 
			polysite[n][3]=n22f; 			
			polysite[n][4]=n11m; //male counts
			polysite[n][5]=n12m; 
			polysite[n][6]=n22m; 			
			n++;
		}
	}

	*npolysites=n;
	if(n==0) {
		return NULL;
	}
	else {
		return polysite;
	}
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

int main(int argc, char *argv[]) {
	//read the reads2snp genotype file and return, for each gene, each polymorphic position with its counts
	//arguments : file name, comma-separated list of females names, comma-separated list of male names
	
	FILE *fp,*outfile;
	int CUR_MAX = 4095;
	int NAME_LEN = 100;
	int NCONTIG_BATCH = 100;
	int ncontigs_allocated;
	char *line = calloc((size_t)CUR_MAX,sizeof(char)); // allocate buffer.
	char *tmpline = calloc((size_t)CUR_MAX,sizeof(char));
	char *contig;
	int count = 0; 
	int length = 0;
	char ch,c;
	int i,j,k,l,ni,firstcontig,pos,t,jmax,lmax,it;
	int nI=0,nJ; //number of individuals
	char **name,**femname,**malname;
	char tmpname[NAME_LEN];
	double true_prob=0,totgenprob=0;
	int gensite=0;
	char *genotypes;
	int ***polysite;
	int *npolysites;
	int pipepos,ichar,nfem,nmal,ifem,imal,ncontigs,totsites;
	int *sex;
	double ****siteprobs,***expS,**expR,sumET,weights,sumnxy,nxy;
	long double pi[NTYPES],oldpi[NTYPES],pidelta[NTYPES],pideltamax,sumpi;
	long double G[NRTYPES],logG[NRTYPES],PStot,logGmax;
	int piorder[NTYPES];
	long double loglik,oldloglik;
	double increase,pmax,oldexpS;
	int mode=CONTIG,plateausteps,nnoncontigs;
	int simul=0;
	double stop=0.000001;

	for(i=0;i<argc;i++) {
		fprintf(stdout,"%s ",argv[i]);
	}
	fprintf(stdout,"\n");
	
	if (argc == 4) {
		simul=1;
		fprintf(stdout,"No individuals names given ; supposing simulated data\n");
	}
	else if (argc != 6) {
		fprintf(stderr,"Usage: %s infile outfile mode namelist1 namelist2\n",argv[0]);
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
	else if (strcmp(argv[3],"s")==0 || strcmp(argv[3],"2")==0) {
		mode=SITE;
	}
	else {
		fprintf(stderr,"Usage: %s infile outfile mode namelist1 namelist2\n",argv[0]);
		fprintf(stderr,"Mode should be either \"c\" or \"1\" for contig-mode, or \"s\" or \"2\" for site-wise optimisation\n");
		exit(1);	
	}
	
	if (!simul) {
		femname=parsenames(argv[4],&nfem,NAME_LEN);
		malname=parsenames(argv[5],&nmal,NAME_LEN);
	}
	
	//genotype file reading
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
	
	while (ch != EOF) { //loop through the genotype file
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
		//		printf("%d\t%d\t%d\n",l,count,length);
		//		printf("%s\n",line);
		//cases : 
		//line starts with ">" : name of contig
		//line starts with "position" (directly after name line) : header specifying the individuals
		//line starts with a number : actual data
		
		
		if ( line[0] == '>' ){
			if (firstcontig==0) {	//Filtering polymorphisms for the last read contig
				if (!simul) {
					nJ=pos;
					polysite[k]=polyfilter(nJ,ni,genotypes,sex,&npolysites[k]);
					free(genotypes);
				}
				else {
					npolysites[k]=t;
					t=0;
				}
				k++;
				if(k % 5000 == 0){
					fprintf(stdout,"Still reading...\n");
				}
			}
					
			if (firstcontig==1 && simul ) {
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
		else if ( line[0] == 'p' && !simul ){
			if (firstcontig) {
				//Important : individuals can be lacking in some contigs !
				//we'll have to count the number of individuals first : the number of words separated by tabs, minus one
				ni=0; //we'll have "position\tname1\tname2\t...\tnameN\0 : the number of tabs is the number of individuals.
				i=0;
				while ((c=line[i]) != '\0'){
					if (c == '\t') 
						ni++;
					i++;
				}
				if(ni>nfem+nmal){
					fprintf(stderr,"Contig %s seems to have more observations than I expected\n",&contig[k*NAME_LEN]);
					exit(1);
				}
				if((name=(char **)calloc((size_t)(nfem+nmal),sizeof(char *)))==NULL) { 
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}
				//find names iteratively ; shorten the line for each name found :
				sscanf(line,"position%[^\n]",tmpline);
				strcpy(line,tmpline);
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
				if((sex=(int *)calloc((size_t)(nfem+nmal),sizeof(int)))==NULL) { 
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}
				for (i=0; i<ni; i++){			
					sex[i]=0;
					for (ifem=0;ifem<nfem;ifem++) {
						if (strcmp(name[i],femname[ifem])==0) {
							sex[i]=FEMALE;
						}
					}
					for (imal=0;imal<nmal;imal++) {
						if (strcmp(name[i],malname[imal])==0) {
							sex[i]=MALE;
						}
					}
					if ( sex[i]==0 ) {
						fprintf(stderr,"Error: no sex found for individual %s\n",name[i]);
						exit(1);
					}		
				}
				firstcontig=0;
			}
			else { //if not the first contig
				//check if the number of individuals is the same
				ni=0; //we'll have "position\tname1\tname2\t...\tnameN\0 : the number of tabs is the number of individuals.
				i=0;
				while ((c=line[i]) != '\0' && line[i] != '\n'){
					if (c == '\t') 
						ni++;
					i++;
				}
				if(ni>nfem+nmal){
					fprintf(stderr,"Contig %s seems to have more observations than I expected\n",&contig[k*NAME_LEN]);
					exit(1);
				}
//				if (ni != nI) {
//					fprintf(stderr,"Problem reading the genotype file: different number of individuals in line %d\n",l);
//					fprintf(stderr,"problematic line :\n");
//					fprintf(stderr,"%s\n",line);
//					fprintf(stderr,"counted %d, expected %d\n",ni,nI);
//					exit(1);					
//				}
				//Find out the sex of the individuals.
				for (i=0; i<ni; i++){			
					sex[i]=0;
					for (ifem=0;ifem<nfem;ifem++) {
						if (strcmp(name[i],femname[ifem])==0) {
							sex[i]=FEMALE;
						}
					}
					for (imal=0;imal<nmal;imal++) {
						if (strcmp(name[i],malname[imal])==0) {
							sex[i]=MALE;
						}
					}
					if ( sex[i]==0 ) {
						fprintf(stderr,"Error: no sex found for individual %s\n",name[i]);
						exit(1);
					}		
				}
				//Check if the order is the same
//				sscanf(line,"position%[^\n]",tmpline);
//				strcpy(line,tmpline);
//				for (i=0; i<nI; i++){
//					sscanf(line,"\t%s%[^\n]",tmpname,tmpline);
//					strcpy(line,tmpline);
//					//strip the names unnecessary prefixes (cut after last pipe)
//					ichar=0;
//					pipepos=0;
//					while (tmpname[ichar]!='\0') {
//						if (tmpname[ichar]=='|') {
//							pipepos=ichar;
//						}
//						ichar++;
//						if(ichar == NAME_LEN) {
//							fprintf(stderr,"Error: one of the individual names seems to have reached the maximum length of %d characters (including prefixes).\n",NAME_LEN);
//							exit(1);
//						}
//					}
//					if(pipepos>0) {
//						for(ichar=pipepos+1;ichar<NAME_LEN;ichar++) {
//							tmpname[ichar-pipepos-1]=tmpname[ichar];
//						}
//					}
//					if (strcmp(tmpname,name[i])!=0){
//						fprintf(stderr,"Problem reading the genotype file: the individuals seem to be in a different order at line %d\n",l);
//						fprintf(stderr,"I didn't think this would happen. Check if this is indeed the case and contact me.\n");
//						exit(1);					
//					}
//				}
			}
			j=0;
		}
		else { //line contains genotype data (or counts in the case of simulated counts)
			if (!simul) {
				sscanf(line,"%d%[^\n]",&pos,tmpline);
				strcpy(line,tmpline);
				if(pos != j+1){
					fprintf(stderr,"Error in input file line number %d: expecting to read \"%d ...\"\n",l,j+1);
					fprintf(stderr,"read \"%s\" instead\n",line);
					exit(1);
				}
				if (j==0) {
					if((genotypes=(char *)malloc(sizeof(char)*ni*2))==NULL){
						fprintf(stderr,"error in memory allocation\n");
						exit(1);
					}
				}
				else {
					if((genotypes=realloc(genotypes,sizeof(char)*ni*2*pos))==NULL){
						fprintf(stderr,"error in memory allocation\n");
						exit(1);
					}
				}
				for (i=0; i<ni; i++){
					sscanf(line,"\t%c%c|%lf%[^\n]",&genotypes[j*ni*2+2*i],&genotypes[j*ni*2+2*i+1],&true_prob,tmpline);
					if(true_prob>2*GSL_DBL_MIN) {
						totgenprob+=true_prob;
						gensite++;
					}
					strcpy(line,tmpline);
				}
				j++;
			}
			else { //simulated data
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
				polysite[k][t][0]=t;
				sscanf(line,"%*d\t%*f\t%d\t%d\t%d\t%d\t%d\t%d\t",&polysite[k][t][1],&polysite[k][t][2],&polysite[k][t][3],&polysite[k][t][4],&polysite[k][t][5],&polysite[k][t][6]);
				t++;
			}
		}
		//		printf("%s\n",genotypes);
		memset(line, '\0', strlen(line)*(sizeof line));
	}
	fprintf(stdout,"%d sites (individuals * positions), error rate (reads2snp) = %f\n",gensite,1.-totgenprob/gensite);
	nJ=pos;
	
//	if (!simul) {
//		for (i=0;i<nI;i++) {
//			printf("%s\t",name[i]);
//		}
//		printf("\n");
//	}
	//Filtering polymorphisms for the last contig
	
	if (!simul) {
		polysite[k]=polyfilter(nJ,ni,genotypes,sex,&npolysites[k]);
	}
	else {
		npolysites[k]=t;
	}
	ncontigs=k+1;
	
//	else { //read simulated data
//		l=0;
//		ch='a';
//		if((npolysites=(int *)malloc(sizeof(int)))==NULL) {
//			fprintf(stderr,"error in memory allocation\n");
//			exit(1);
//		}
//		if((polysite=(int ***)malloc(sizeof(int **)))==NULL) {
//			fprintf(stderr,"error in memory allocation\n");
//			exit(1);
//		}
//		ncontigs_allocated=1;
//		
//		while (ch != EOF) { //loop through the genotype file
//		ch='a';
//		count = 0;
//		length = 0;
//		while ( (ch != '\n') && (ch != EOF) ) { //loop through the line
//			if(count == CUR_MAX) { // time to expand (for unexepectedly large line lengths) ?
//				CUR_MAX *= 2; 
//				count = 0;
//				line = realloc(line, sizeof(char) * CUR_MAX); 
//				tmpline = realloc(tmpline, sizeof(char) * CUR_MAX); 
//			}
//			ch = getc(fp); // read from stream.
//			line[length] = ch;
//			length++;
//			count++;
//		}
//		line[length] = '\0';
//		if (length <= 1) { //empty line : suppose it's the end of the file
//			break;
//		}
//		//We've read one line
//		if (l==0) {
//			if((polysite[0]=(int **)malloc(sizeof(int *)))==NULL){
//				fprintf(stderr,"error in memory allocation\n");
//				exit(1);
//			}
//			if((polysite[0][l]=(int *)malloc(sizeof(int)*7))==NULL){
//				fprintf(stderr,"error in memory allocation\n");
//				exit(1);
//			}
//		}
//		else {
//			if((polysite[0]=realloc(polysite[0],sizeof(int *)*(l+1)))==NULL){
//				fprintf(stderr,"error in memory allocation\n");
//				exit(1);
//			}
//			if((polysite[0][l]=(int *)malloc(sizeof(int)*7))==NULL){
//				fprintf(stderr,"error in memory allocation\n");
//				exit(1);
//			}
//		}
//		polysite[0][l][0]=l;
//		sscanf(line,"%*d\t%*f\t%d\t%d\t%d\t%d\t%d\t%d\t",&polysite[0][l][1],&polysite[0][l][2],&polysite[0][l][3],&polysite[0][l][4],&polysite[0][l][5],&polysite[0][l][6]);
//		l++;
//		}
//		npolysites[0]=l;
//		ncontigs=1;
//	}
	
	//All reading has been done : clear memory.
	free(tmpline);
	free(line);
	if (!simul) {
		for (i=0; i<nI; i++){
			free(name[i]);
		}
		free(name);
		for (i=0; i<nfem; i++){
			free(femname[i]);
		}
		free(femname);
		for (i=0; i<nmal; i++){
			free(malname[i]);
		}
		free(malname);
		free(genotypes);
		free(sex);
	}
	fclose(fp);

	//end of data reading and filtering
	
	fprintf(stdout,"Found %d contigs:\n",ncontigs);
	totsites=0;
	nnoncontigs=0;
	for (k=0;k<ncontigs;k++) {
		totsites+=npolysites[k];
		if(npolysites[k]==0) {
			nnoncontigs++;
		}
//		fprintf(stdout,"%d: %s (%d polymorphic positions)\n",k+1,&contig[k*NAME_LEN],npolysites[k]);
	}
	fprintf(stdout,"...and %d polymorphic sites\n",totsites);
	if(totsites==0){
		fprintf(stdout,"No polymorphic sites found; nothing to do.\n");
		exit(0);
	}

	//EM algorithm
	
	// initialization
	if((siteprobs=(double ****)calloc((size_t)ncontigs,sizeof(double ***)))==NULL) { 
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}
	for (k=0;k<ncontigs;k++) {
		siteprobs[k]=initEM(npolysites[k],polysite[k]);	
	}
	if (mode == SITE) {
		for(j=0;j<NTYPES;j++) {
			pi[j]=0.25;
		}
		pisort(NTYPES,pi,piorder);
	}
	else {
		for(j=0;j<NRTYPES;j++) {
			pi[j]=1./3.;
		}
		pi[NRTYPES]=0.5;
		pisort(NRTYPES,pi,piorder);
	}		
	
	if((expS=(double ***)calloc((size_t)ncontigs,sizeof(double **)))==NULL) { 
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}
	for (k=0;k<ncontigs;k++) {
		if (npolysites[k]>0) {
			if((expS[k]=(double **)calloc((size_t)npolysites[k],sizeof(double *)))==NULL) { 
				fprintf(stderr,"error in memory allocation\n");
				exit(1);
			}
			for (t=0; t<npolysites[k]; t++){
				if((expS[k][t]= (double *)calloc((size_t)NTYPES,sizeof(double)))==NULL) { 
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}
			}
		}
	}
	if((expR=(double **)calloc((size_t)ncontigs,sizeof(double *)))==NULL) { 
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}
	for (k=0;k<ncontigs;k++) {
		if (npolysites[k]>0) {
			if((expR[k]=(double *)calloc((size_t)NRTYPES,sizeof(double)))==NULL) { 
				fprintf(stderr,"error in memory allocation\n");
				exit(1);
			}
		}
	}
		
	
	for (t=0; t<npolysites[0]; t++){
		for (i=0; i<7; i++){
			printf("%d ",polysite[0][t][i]);
		}
		for(j=0;j<NTYPES;j++){
		printf("%f\t%f\t%f\t",siteprobs[0][t][0][j],siteprobs[0][t][1][j],siteprobs[0][t][2][j]);
		}
		printf("\n");
		if (t>500) {
			break;
		}
	}
	
	//initial likelihood
	if (mode == SITE) {
		loglik=totalsiteloglik(ncontigs,npolysites,pi,siteprobs,piorder);
	}
	else {
		loglik=totalcontigloglik(ncontigs,npolysites,pi,siteprobs);
	}

	fprintf(stdout,"Initial log-likelihood: %Lf\n",loglik);
	
	//Iteration starts
	increase=1;
	it=0;
	plateausteps=0;
	while(plateausteps<10){
		//	while(l<10){
		it++;
		if (mode==SITE) {		//site-wise
			//E-step
//			loglik=0;
			for (k=0;k<ncontigs;k++) {
				for(t=0;t<npolysites[k];t++){
					//				printf("pos %d\t",i);
					// Horner decomposition (could be faster by storing jmax, which doesn't change throughout optimisation)
					pmax=0.0;
					jmax=-1;
					for(i=0;i<NTYPES;i++) {
						j=piorder[i];
						if (siteprobs[k][t][2][j]>pmax) {
							jmax=j;
							pmax=siteprobs[k][t][2][j];
						}
					}
					PStot=pi[jmax];
					//							if (isnan(PStot)) {
					//								fprintf(stderr,"Nan observed at contig %d, site %d\n",k,t);
					//							}
					for(i=0;i<NTYPES;i++) {
						j=piorder[i];
						if (j!=jmax) {
							PStot+=pi[j]*siteprobs[k][t][2][j]/pmax;
							//							if (isnan(PStot)) {
							//								fprintf(stderr,"Nan observed at contig %d, site %d, type %d, PS=%e\n",k,t,j,siteprobs[k][t][1][j]);
							//							}
						}
					}
					PStot*=pmax;
					for(i=0;i<NTYPES;i++) {
						j=piorder[i];
						oldexpS=expS[k][t][j];
						expS[k][t][j]=pi[j]*siteprobs[k][t][2][j]/PStot;
						//					printf("%f\t",expS[k][i][type]);
//						if(plateausteps>2){
//							printf("%f\t%f\t%e\t%e\n",expS[k][t][j],oldexpS,expS[k][t][j]-oldexpS,oldexpS-nextafter(oldexpS,oldexpS+1.0));
//						}
					}
//					loglik+=log(PStot);
					//				printf("\n");
				}
			}
//					fprintf(stdout,",End of E-step : log-likelihood: %f\n",loglik);
			//M-step
			sumpi=0;
			for(i=0;i<NTYPES;i++) {
				j=piorder[i];
				oldpi[j]=pi[j];
				sumET=0;
				totsites=0;
				for (k=0;k<ncontigs;k++) {
					for(t=0;t<npolysites[k];t++){
						sumET+=expS[k][t][j];
					}
					totsites+=npolysites[k];
				}
				pi[j]=(long double)sumET/(long double)totsites;
				sumpi+=pi[j];
			}
//			fprintf(stdout,"Iteration %d: pi",l);
//			for(j=0;j<NTYPES;j++) {
//				fprintf(stdout," %Lf (%Le)",pi[j],logl(pi[j]));
//			}
//			fprintf(stdout,"; total: %Lf",sumpi);
			//reordering the labels
			pisort(NTYPES,pi,piorder);
//			fprintf(stdout,"\n");
//			for (i=0; i<NTYPES; i++) {
//				fprintf(stdout,"%d\t%f\t",piorder[i],pi[piorder[i]]);
//			}
//			fprintf(stdout,"\n");

			//		printf(", sums %f (%f)",pi[0]+pi[1]+pi[2]+pi[3],exp2log(pi[0])+exp2log(pi[1])+exp2log(pi[2])+exp2log(pi[3]));
			// 
			pideltamax=0.0;
			for(i=0;i<NTYPES;i++) {
				j=piorder[i];
				pidelta[j]= pi[j]>stop ? (pi[j]-oldpi[j])/pi[j] : 0.0;
				if (fabsl(pidelta[j]) > fabsl(pideltamax)) {
					pideltamax=pidelta[j];
				}
			}	
			if (fabsl(pideltamax) < stop) {
				plateausteps++;
			}
			fprintf(stdout,"Iteration %d: pi",it);
			for(j=0;j<NTYPES;j++) {
				fprintf(stdout," %Lf (%Le)",pi[j],pidelta[j]);
			}
			fprintf(stdout,"; total: %Lf",sumpi);
			//new log-likelihood
			oldloglik=loglik;
			loglik=totalsiteloglik(ncontigs,npolysites,pi,siteprobs,piorder);
			fprintf(stdout,", log-likelihood: %Lf ",loglik);

			
			increase=loglik-oldloglik;
//			if (increase<=0.000001) {
//				plateausteps++;
//			}
		}
		else if (mode == CONTIG ) { //contig-wise
			//E-step
			for (k=0;k<ncontigs;k++) {
				if (npolysites[k]>0) {
					//autosomal and X-hemizygous: calculate product of site probabilities
					//do this in log, as the products rapidly shrink when the number of sites increases
					for(j=0;j<2;j++) {
						logG[j]=0.;
						for(t=0;t<npolysites[k];t++){
							logG[j]+=siteprobs[k][t][1][j];
						}
					}
					//XY: calculate product of weighted average of site probabilities of both types
					//s=pi[3]
					logG[2]=0.;
					for(t=0;t<npolysites[k];t++){
						logG[2]+=log(pi[3]*siteprobs[k][t][2][2]+(1.-pi[3])*siteprobs[k][t][2][3]);
					}
					//expectations
					PStot=0;
					//Horner decomposition
					logGmax=log(0);
					for(l=0;l<NRTYPES;l++){
						if(logG[l]>logGmax){
							logGmax=logG[l];
							lmax=l;
						}
					}
					PStot=pi[lmax];
					for(l=0;l<NRTYPES;l++){
						if(l!=lmax){
							PStot+=exp(log(pi[l])+logG[l]-logG[lmax]);
						}
					}
					//PStot will be in log from here on !!
					PStot=log(PStot);
					PStot+=logG[lmax];
					for(l=0;l<NRTYPES;l++){
						expR[k][l]=exp(log(pi[l])+logG[l]-PStot);
						if (isnan(expR[k][l])) {
							fprintf(stderr,"Nan observed at contig %d\n",k);
							fprintf(stderr,"type %d, PStot %Le, pi %Lf\n",l,PStot,pi[l]);
							fprintf(stderr,"contig has %d polymorphic sites\n",npolysites[k]);
							for(l=0;l<NRTYPES;l++){
								fprintf(stderr,"logG[%d] : %Le\t",l,logG[l]);
							}
							fprintf(stderr,"\n");
							exit(1);
						}
					}
				}
			}
			//M-step
			// 1) rho-values (pi[0]-pi[2])
			for(j=0;j<NTYPES;j++) {
				oldpi[j]=pi[j];
			}
			sumpi=0;
			for(l=0;l<NRTYPES;l++){
				sumET=0;
				for (k=0;k<ncontigs;k++) {
					if (npolysites[k]>0) {					
						sumET+=expR[k][l];
						if (isnan(sumET)) {
							fprintf(stderr,"Nan observed at contig %d\n",k);
							fprintf(stderr,"type %d, expectation %f\n",l,expR[k][l]);
							exit(1);
						}
					}
				}
				pi[l]=sumET/(ncontigs-nnoncontigs);
				sumpi+=pi[l];
			}
			// 2) s (pi[3])
			if ( pi[2] > LDBL_MIN ){
				sumET=0;
				weights=0;
				for (k=0;k<ncontigs;k++) {
					if (npolysites[k]>0) {					
						sumnxy=0;
						for(t=0;t<npolysites[k];t++){
							if(siteprobs[k][t][2][2]+siteprobs[k][t][2][3] > GSL_DBL_MIN ) {
								sumnxy+=pi[3]*siteprobs[k][t][2][2]/(pi[3]*siteprobs[k][t][2][2]+(1.-pi[3])*siteprobs[k][t][2][3]);
								//sumnxy+=1./(1.+(1.-pi[3])*siteprobs[k][t][2][3]/(pi[3]*siteprobs[k][t][2][2]));
								//printf("%Le\t%e\t%e\n",pi[3],siteprobs[k][t][2][2],siteprobs[k][t][2][3]);
								if (isnan(sumnxy)) {
									fprintf(stderr,"Nan observed at contig %d, site %d\n",k,t);
									exit(1);
								}
							}
							else {
								sumnxy+=pi[3];
							}
						}
						sumET+=expR[k][2]*sumnxy;
						weights+=expR[k][2]*npolysites[k];
					}
				}
				pi[3]=sumET/weights;
				if (isnan(pi[3])) {
					printf("%Le\t%e\t%e\n",pi[3],sumET,weights);
				}
			}
			else {
				pi[3]=0.5;
			}

			//calculate improvement
			pideltamax=0.0;
			for(j=0;j<NTYPES;j++) {
				pidelta[j]= pi[j]>stop ? (pi[j]-oldpi[j])/pi[j] : 0.0;
				if (fabsl(pidelta[j]) > fabsl(pideltamax)) {
					pideltamax=pidelta[j];
				}
			}	
			if (fabsl(pideltamax) < stop) {
				plateausteps++;
			}
			fprintf(stdout,"Iteration %d: pi",it);
			for(j=0;j<NTYPES;j++) {
				fprintf(stdout," %Lf (%Le)",pi[j],pidelta[j]);
			}
			fprintf(stdout,"; total: %Lf",sumpi);
			//new log-likelihood
			oldloglik=loglik;
			loglik=totalcontigloglik(ncontigs,npolysites,pi,siteprobs);
			fprintf(stdout,", log-likelihood: %Lf ",loglik);

			
			increase=loglik-oldloglik;
			
			
		}
		printf("increase: %e (precision: %Le)\n",increase,loglik-nextafterl(loglik,loglik+100000.0));
	}
	
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
				expR[k][j]=sumET/(t+1);
			}
			//XY: calculate average of summed probabilities of both types
			sumET=0;
			for (t=0; t<npolysites[k]; t++){
				for(j=2;j<NTYPES;j++) {
					sumET+=expS[k][t][j];
				}
			}
			expR[k][2]=sumET/(t+1);
			}
		}
	}
//	else {
//		//calculate posterior probabilities per site
//		for (k=0;k<ncontigs;k++) {
//			for (t=0; t<npolysites[k]; t++){
//				sumET=0;
//				for(l=0;l<2;l++){
//					sumET+=expR[k][l]*siteprobs[k][t][2][l];
//				}
//				sumET+=expR[k][2]*(pi[3]*siteprobs[k][t][2][2]+(1.-pi[3])*siteprobs[k][t][2][3]);
//				for(j=0;j<2;j++){
//					expS[k][t][j]=expR[k][j]*siteprobs[k][t][2][j]/sumET;
//				}
//				expS[k][t][2]=expR[k][2]*pi[3]*siteprobs[k][t][2][2]/sumET;
//				expS[k][t][3]=expR[k][2]*(1.-pi[3])*siteprobs[k][t][2][3]/sumET;
//			}
//		}
//	}
	
	//calculate average Fst per contig
	//Fst for each site
//	for (k=0;k<ncontigs;k++) {
//		for(t=0;t<npolysites[k];t++) {
//			n11f=polysite[i][1]; //female counts
//			n12f=polysite[i][2]; 
//			n22f=polysite[i][3]; 			
//			n11m=polysite[i][4]; //male counts
//			n12m=polysite[i][5]; 
//			n22m=polysite[i][6]; 			
//			nfem=n11f+n12f+n22f;
//			nmal=n11m+n12m+n22m;
//			ntot=nfem+nmal;
//			if (nfem !=0 && nmal !=0) { 
//			f0f=(double)(2*n11f+n12f)/(double)(2*nfem);
//			f0m=(double)(2*n11m+n12m)/(double)(2*nmal);
//			fst=(f0f*f0f+f0m*f0m)/2.-
//			}
//			else {
//				fst=-999;
//			}
//		}
//	}
	
	
	
	fprintf(outfile,"log-likelihood: %Lf, pi:",loglik);
	for(j=0;j<NTYPES;j++) {
		fprintf(outfile,"\t%Lf",pi[j]);
	}
	fprintf(outfile,"\n");
	
	for (k=0;k<ncontigs;k++) {
		fprintf(outfile,">%s\t%d",&contig[k*NAME_LEN],npolysites[k]);
		if(npolysites[k]>0) {
		pmax=0.0;
		lmax=-1;
		for(l=0;l<NRTYPES;l++) {
			if (expR[k][l]>pmax) {
				lmax=l;
				pmax=expR[k][l];
			}
			fprintf(outfile,"\tl : %d; pp : %f",l,expR[k][l]);
		}
		fprintf(outfile,"\t%d\n",lmax);
//		if(mode==SITE){
			for (t=0; t<npolysites[k]; t++){
				for (i=0; i<7; i++) {
					fprintf(outfile,"%d\t",polysite[k][t][i]);
				}
				pmax=0.0;
				jmax=-1;
				for(j=0;j<NTYPES;j++) {
					if (expS[k][t][j]>pmax) {
						jmax=j;
						pmax=expS[k][t][j];
					}
					if(mode==SITE){
						fprintf(outfile,"%f\t",expS[k][t][j]);
					}
					else {
						fprintf(outfile,"%f\t",siteprobs[k][t][2][j]);
					}
					fprintf(outfile,"(%f)\t",siteprobs[k][t][0][j]);
				}
				fprintf(outfile,"%d\n",jmax);
			}
//		}
		}
		else {
			fprintf(outfile,"\n");
		}
	}
	
	fclose(outfile);
	
	for (k=0;k<ncontigs;k++) {	
		freeEM(npolysites[k], siteprobs[k]);
		for (t=0; t<npolysites[k]; t++){
			free(expS[k][t]);
			free(polysite[k][t]);
		}
		if(npolysites[k]>0){
			free(polysite[k]);
			free(expS[k]);
			free(expR[k]);
		}
	}
	free(expS);
	free(expR);
	free(polysite);
	free(npolysites);
	free(siteprobs);
	free(contig);
	
	return 0;
}