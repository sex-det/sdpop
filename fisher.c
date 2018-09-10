#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "reading.h"
#include "calc.h"


//Fisher's exact test
//From : Korkuc & Walther (2016) The identification of cis-regulatory sequence motifs in gene promoters based on SNP information.  Methods in Mol. Biol.
//Additional information and software utilities. http://bioinformatics.mpimp-golm.mpg.de/research-projects-publications/supplementary-data/walther/go-term-enrichment-analysis-1/fisher-exact-test-c-code

double factorInc(int chi11, int chi12, int chi21, int chi22)
{
  double factor_inc;
  factor_inc = (double)chi12 * chi21;
  factor_inc /= (double)(chi11 + 1) * (chi22 + 1);
  return factor_inc;
}

double factorDec(int chi11, int chi12, int chi21, int chi22)
{
  double factor_dec = (double)chi11 * chi22;
  factor_dec /= (double)(chi21 + 1) * (chi12 + 1);
  return factor_dec;
}
 
double gammln(double xx)
{
  static double cof[6] = { 76.18009172947146, -86.50532032941677, 
			  24.01409824083091, -1.231739572450155,
                          0.1208650973866179e-2, -0.5395239384953e-5 };
  double x, tmp, ser;
  int j;

  x = xx - 1.0;
 
  tmp = x + 5.5;
  tmp -= (x + 0.5) * log(tmp);
  ser = 1.0;
  for (j = 0; j <= 5; j++) {
    x += 1.0;
    ser += cof[j] / x;
  }
  return -tmp + log(2.50662827465 * ser);
}
 
double factln0(int n)
{
  static double pi=3.1415926536; 
 
  if (n < 0) {
     //nrerror("Negative factorial in routine FACTLN");
     return 0.0;
  }
  if (n <= 1) return 0.0;
  return 0.5*log((2.0*n+1./3.)*pi)+n*(log(1.0*n)-1.0);
}
 
double factln(int n)
{
  static double a[101];
 
 
  if (n < 0) {
     //nrerror("Negative factorial in routine FACTLN");
     return 0.0;
  }
  if (n <= 1) return 0.0;
  if (n <= 100) return a[n] ? a[n] : (a[n] = gammln((double)(n + 1.0)));
  else return gammln((double)(n + 1.0));
}
 
double binomialln(int n, int k)
{
  return (factln(n) - factln(k) - factln(n - k));
 
}
 
double calc_hypergeom(int chi11, int chi12, int chi21, int chi22)
{
   static double choose1, choose2, choose3,choose0;
   static double b1, b2, b3;
   static int total;
 
   total = chi11 + chi12 + chi21 + chi22;
 
   b1 = binomialln(chi11 + chi12, chi11);
   b2 = binomialln(chi21 + chi22, chi21);
   b3 = binomialln(total, chi11 + chi21);
 
//   choose1 = exp(b1);
//   choose2 = exp(b2);
//   choose3 = exp(b3);
 
   choose0=exp(b1+b2-b3);
 
// Return the hypergeometric
//   return (choose1 * choose2 / choose3);
   return (choose0);
}
 
double fisher_exact_tiss(int chi11, int chi12, int chi21, int chi22)
{
   register int co_occ, total_libs;
   static int min_co_occ, max_co_occ;
   static int gene_a, gene_b;
   static double factor_inc, factor_dec;
   static double base_p, curr_p;
   int sign=1;
 
   co_occ = chi11;
 
   gene_a = chi11 + chi12;
   gene_b = chi11 + chi21;
   total_libs = chi11 + chi12 + chi21 + chi22;
 
   // If the two genes occur few enough times, the minimum number of
   // co-occurrences is 0.  If the total number of times they occur
   // exceeds the number of libraries (say by N), they must overlap
   // at least N times.
   min_co_occ = 0;
   if (gene_a + gene_b > total_libs) {
      min_co_occ = gene_a + gene_b - total_libs;
   }
 
   // Maximum number of co-occurrences is at most the number of times
   // the rarer gene occurs in the library :
 
   if (gene_a < gene_b) { max_co_occ = gene_a; }
   else { max_co_occ = gene_b; }
 
   // Calculate the first hypergeometric value
 
   base_p = calc_hypergeom(chi11, chi12, chi21, chi22);
 
// printf("base_p=%e\n",base_p);
 
   // If co-occurrences at max possible, then this is our p-value,
   // Also if co-occurrences at min possible, this is our p-value.
 
   if (co_occ == max_co_occ || co_occ == min_co_occ) 
     {
       if(co_occ == max_co_occ) {
	 sign=1;
       }
       else {
	 sign=-1;
       }
     }
   else {
      // Need to add in the other possible p-values.
      factor_inc = factorInc(chi11, chi12, chi21, chi22);
      factor_dec = factorDec(chi11, chi12, chi21, chi22);
 
// printf("%e,%e:%e\n",factor_inc,factor_dec,base_p);
 
      // Start out with the current p-value
      curr_p = base_p;
 
      // Want to sum the probabilites in the direction of decreasing P
      if (factor_dec < factor_inc) {  sign=-1;
         // Loop down over co-occurrences
         do {
            // Determine P-value for current chi^2 matrix
            curr_p *= factor_dec;
 
            // Add to probability based on recurrence factor
            base_p += curr_p;
            co_occ--;
 
            // Alter chi^2 matrix to reflect number of co-occurrences
            chi11--; chi22--;
            chi12++; chi21++;
 
            // Get the next value for the recurrence factor
            factor_dec = factorDec(chi11, chi12, chi21, chi22);
         } while (co_occ > min_co_occ);
 
      }
      else if (factor_inc < factor_dec) { sign=1;
         // Loop up over co-occurrences
         do {
            // Determine P-value for chi^2 matrix from recurrence factor
            curr_p *= factor_inc;
 
            // Add to probability based on recurrence factor
            base_p += curr_p;
            co_occ++;
 
            // Alter chi^2 matrix to reflect number of co-occurrences
            chi11++; chi22++;
            chi12--; chi21--;
 
            // Get the next value for the recurrence factor
            factor_inc = factorInc(chi11, chi12, chi21, chi22);
// printf("[%d %d %d %d]\t%e\t%e\t%e\n",chi11, chi12, chi21, chi22, curr_p,base_p,factor_inc);
         } while (co_occ < max_co_occ);
      }
      else {
         // We are on a saddle point, which means p-value is 1.
         base_p = 1.0;
      }
   }
// remove strange convention to return negative p-values
//   return base_p*sign;
   return base_p;
}

//End Fisher's exact test

int read_cnt0(FILE *fp, 	int NAME_LEN, char **contig_p, int **npolysites_p, int ****polysite_p) 
{
	char *contig;
	int *npolysites,***polysite;
	int CUR_MAX = 4095;
	int NCONTIG_BATCH = 100;
	int ncontigs_allocated,ncontigs;
	char *line = calloc((size_t)CUR_MAX,sizeof(char));
	char *tmpline = calloc((size_t)CUR_MAX,sizeof(char));
	int count = 0; 
	int length = 0;
	char ch;
	int i,k,l,firstcontig,t;
	int ni,nI=0;


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
//				if(chromosomes==XY){
					sscanf(line,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t",&polysite[k][t][0],&polysite[k][t][1],&polysite[k][t][2],&polysite[k][t][3],&polysite[k][t][4],&polysite[k][t][5],&polysite[k][t][6]);
//				}
//				else {
//					sscanf(line,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t",&polysite[k][t][0],&polysite[k][t][4],&polysite[k][t][5],&polysite[k][t][6],&polysite[k][t][1],&polysite[k][t][2],&polysite[k][t][3]);
//				}
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
	
	*npolysites_p=npolysites;
	*polysite_p=polysite;
	*contig_p=contig;
	
	return ncontigs;
}

int main(int argc, char *argv[]) {

	FILE *fp,*outfile;
	int NAME_LEN = 900;
	char *contig;
	int i,j,k,l,jl,t,it,s,g,gp,ni,nI;
	int ***polysite;
	int *npolysites;
	double **fisher_p,**fst;
	double pmin,pval;
	double f,fm,ff,a,b,alpham,alphaf,suma,sumab;
	int n11f,n12f,n22f,n11m,n12m,n22m,n,nm,nf;
	int ncontigs,nnoncontigs,totsites;

	for(i=0;i<argc;i++) {
		fprintf(stdout,"%s ",argv[i]);
	}
	fprintf(stdout,"\n");
	
	if (argc != 3) {
		fprintf(stdout,"Usage: %s infile outfile\n",argv[0]);
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

//	ncontigs=read_cnt(fp,NAME_LEN,&contig,&npolysites,&polysite);
	ncontigs=read_cnt2(fp,NAME_LEN,2,&contig,&npolysites,&polysite);
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
	
	if((fisher_p=(double **)calloc((size_t)ncontigs,sizeof(double *)))==NULL) { 
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}
	for (k=0;k<ncontigs;k++) {
		if (npolysites[k] > 0) {
			if((fisher_p[k]=(double *)calloc((size_t)npolysites[k],sizeof(double)))==NULL) { 
				fprintf(stderr,"error in memory allocation\n");
				exit(1);
			}
		}
	}
	if((fst=(double **)calloc((size_t)ncontigs,sizeof(double *)))==NULL) { 
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}
	for (k=0;k<ncontigs;k++) {
		if (npolysites[k] > 0) {
			if((fst[k]=(double *)calloc((size_t)npolysites[k],sizeof(double)))==NULL) { 
				fprintf(stderr,"error in memory allocation\n");
				exit(1);
			}
		}
	}
	
	//Fisher test & Fst
	
	for (k=0;k<ncontigs;k++) {
		fprintf(outfile,">%s\t%d",&contig[k*NAME_LEN],npolysites[k]);
		if(npolysites[k]>0) {
			pmin=1.;
			suma=0;
			sumab=0;
			for (t=0; t<npolysites[k]; t++){
				n11f=polysite[k][t][1]; //female counts
				n12f=polysite[k][t][2]; 
				n22f=polysite[k][t][3]; 			
				n11m=polysite[k][t][4]; //male counts
				n12m=polysite[k][t][5]; 
				n22m=polysite[k][t][6];
				nf=n11f+n12f+n22f;
				nm=n11m+n12m+n22m;
				n=nf+nm;
				//Fst (equations 1-4 in Fumagalli et al., 2013, Genetics 195: 979--992)
				ff=(2.*n11f+n12f)/(2.*nf);
				fm=(2.*n11m+n12m)/(2.*nm);
				f=(2.*n11f+2.*n11m+n12f+n12m)/(2.*nf+2.*nm);
				alphaf=2.*ff*(1.-ff);
				alpham=2.*fm*(1.-fm);
				b=(nf*alphaf+nm*alpham)/(nf+nm-1.);
//				b=(nf*alphaf+nm*alpham)/(nf+nm);
				a=(4.*nf*intpow(ff-f,2)+4.*nm*intpow(fm-f,2)-b)/(2.*((2.*nf*nm)/(nf+nm)));
				fst[k][t]=a/(a+b);
//				for (i=0; i<7; i++) {
//					printf("%d\t",polysite[k][t][i]);
//				}
//				printf("%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",nf,nm,n,ff,fm,f,alphaf,alpham,b,a,fst[k][t]);
				suma+=a;
				sumab+=(a+b);
				//fisher
				pval=fisher_exact_tiss(n11f+n22f,n12f,n11m+n22m,n12m);
				if(pval<pmin) {
					pmin=pval;
				}
				fisher_p[k][t]=pval;
			}
			fprintf(outfile,"\t%e\t%f\n", 1. - intpow(1.-pmin,npolysites[k]), suma/sumab);
			for (t=0; t<npolysites[k]; t++){
				for (i=0; i<7; i++) {
					fprintf(outfile,"%d\t",polysite[k][t][i]);
				}
				fprintf(outfile,"%e\t%f\n",fisher_p[k][t],fst[k][t]);
			}
		}
		else {
			fprintf(outfile,"\t%e\n",1.);
		}
	}
	
	fclose(outfile);
	free(contig);
	for(k=0;k<ncontigs;k++) {
		if(npolysites[k]>0) {
			for(t=0; t<npolysites[k]; t++){
				free(polysite[k][t]);
			}
			free(polysite[k]);
			free(fisher_p[k]);
			free(fst[k]);
		}
	}
	free(fisher_p);
	free(fst);
	free(polysite);
	free(npolysites);

	
	return 0;
}