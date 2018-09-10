#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define N11 1
#define N12 2
#define N22 3

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

int **polyfilter(int npos, int nind, char *genotypes, int *npolysites, int *n3, int *n4, double *theta) {//number of positions and individuals,
	// genotype matrix, pointer to the number of polymorphic sites
	// function treats one contig
	char nucvec[4] = {0};
	int i,j,randbit,tempi;
	int div,n11,n12,n22,nA,nT,nG,nC,nN;
	int **polysite=NULL;
	char pos1,pos2,nuc1,nuc2;
	int n,n2;
	double a;
	
	*n3=0;
	*n4=0;
	*theta=0.;
	
	n=0;
	for (j=0; j<npos; j++){
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
		if (div == 2 || div==1 ) { //simple polymorphism
			if (n==0) {
				if((polysite=(int **)malloc(sizeof(int *)))==NULL){
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}
				if((polysite[n]=(int *)malloc(sizeof(int)*(nind+4)))==NULL){
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}
			}
			else {
				if((polysite=realloc(polysite,sizeof(int *)*(n+1)))==NULL){
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}
				if((polysite[n]=(int *)malloc(sizeof(int)*(nind+4)))==NULL){
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}
			}
			polysite[n][0]=j+1; //position
			for (i=0; i<nind+3; i++){
				polysite[n][i+1]=0;				
			}
		}
		if (div == 1) {
			for (i=0; i<nind; i++){ //loop through all chromosomes
				pos1=genotypes[j*nind*2+2*i];
				pos2=genotypes[j*nind*2+2*i+1];
				if (pos1 != 'n' && pos1 != 'N'){
					polysite[n][i+1]=N11;				
					polysite[n][nind+N11]++;
				}
			}
			n++;
		}
		else if (div == 2) { //simple polymorphism
			//randomize 
			nuc1=nucvec[0];
			nuc2=nucvec[1];
			randbit=rand()%2;
			nucvec[randbit]=nuc1;
			nucvec[!randbit]=nuc2;
			n11=0;
			n12=0;
			n22=0;
			for (i=0; i<nind; i++){ //loop through all chromosomes
				pos1=genotypes[j*nind*2+2*i];
				pos2=genotypes[j*nind*2+2*i+1];
				if (pos1 == pos2 && pos1 == nucvec[0]) { //homozygous
					n11++;
					polysite[n][i+1]=N11;
				}
				else if (pos1 == pos2 && pos1 == nucvec[1]) { //homozygous
					n22++;
					polysite[n][i+1]=N22;
				}
				else if (pos1 != pos2) { //heterozygous
					n12++;
					polysite[n][i+1]=N12;
				}
			}
			//Watterson's estimate of theta
			a=1.;
			for (i=2; i<2*(n11+n22+n12); i++){
				a+=1./(double)i;
			}
			*theta+=1./a;
			n2++;		
			polysite[n][nind+N11]=n11;
			polysite[n][nind+N12]=n12; 
			polysite[n][nind+N22]=n22; 			
			n++;
		}
		else if (div ==3) {
			*n3++;
		}
		else if (div ==4) {
			*n4++;
		}
	}

	*theta/=(double)npos;
	*npolysites=n;
	if(n==0) {
		return NULL;
	}
	else {
		return polysite;
	}
}

int main(int argc, char *argv[]) 
{
	FILE *fp,*outfile;
	int CUR_MAX = 4095;
	int NAME_LEN = 100;
	int NCONTIG_BATCH = 100;
	int ncontigs_allocated;
	char *line = calloc((size_t)CUR_MAX,sizeof(char));
	char *tmpline = calloc((size_t)CUR_MAX,sizeof(char));
	char *contig;
	int ***polysite;
	int *npolysites,*n3,*n4;
	double *theta;
	int ncontigs,totsites=0;
	int count = 0; 
	int length = 0;
	int i,j,k,l,ni,fi,firstcontig,pos,t,itemp;
	int pipepos,ichar,nnames,iname,nplus=0;
	char ch,c;
	char **name,**names;
	int *found;
	int nJ,nfound; //number of individuals
	char *genotypes;
	int nsites,nhetsites;

	for(i=0;i<argc;i++) {
		fprintf(stdout,"%s ",argv[i]);
	}
	fprintf(stdout,"\n");
	
	if (argc != 4) {
		fprintf(stdout,"Usage: %s infile outfile namelist\n",argv[0]);
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
	
	names=parsenames(argv[3],&nnames,NAME_LEN);

	//genotype file reading
	l=0; //line number
	ch='a';
	if((npolysites=(int *)malloc(sizeof(int)*NCONTIG_BATCH))==NULL) {
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}
	if((n3=(int *)malloc(sizeof(int)*NCONTIG_BATCH))==NULL) {
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}
	if((n4=(int *)malloc(sizeof(int)*NCONTIG_BATCH))==NULL) {
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}
	if((theta=(double *)malloc(sizeof(double)*NCONTIG_BATCH))==NULL) {
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
		
		if ( line[0] == '>' ){
			if (firstcontig==0) {	//Filtering polymorphisms for the last read contig
				nJ=pos;
				polysite[k]=polyfilter(nJ,nfound,genotypes,&npolysites[k],&n3[k],&n4[k],&theta[k]);
				free(genotypes);
				totsites+=npolysites[k];
				k++;
				if(k % 5000 == 0){
					fprintf(stdout,"%d contigs, %d polymorphic sites, and still reading...\n",k,totsites);
				}
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
				if((npolysites=(int *)realloc(npolysites,sizeof(int)*ncontigs_allocated))==NULL) {
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}		
				if((n3=(int *)realloc(n3,sizeof(int)*ncontigs_allocated))==NULL) {
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}		
				if((n4=(int *)realloc(n4,sizeof(int)*ncontigs_allocated))==NULL) {
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}		
				if((theta=(double *)realloc(theta,sizeof(double)*ncontigs_allocated))==NULL) {
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
		else if ( line[0] == 'p' ){
			if (firstcontig) {
				firstcontig=0;
			}
			else {
				for (i=0; i<ni; i++){
					free(name[i]);
				}
				free(name);
				free(found);
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
			if((found=(int *)calloc((size_t)(ni),sizeof(int)))==NULL) {
				fprintf(stderr,"error in memory allocation\n");
				exit(1);
			}
			for (i=0;i<ni;i++) {
				found[i]=0;
			}
			
			//find names iteratively ; shorten the line for each name found :
			sscanf(line,"position%[^\n]",tmpline);
			strcpy(line,tmpline);
			nfound=0;
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
				for(iname=0;iname<nnames;iname++) {
					if (strcmp(name[i],names[iname])==0) {
						found[i]=1;
						nfound++;
						break;
					}
				}
			}

			if(ni>nfound){
				if(nplus==0) {
					fprintf(stdout,"Contig %s seems to have more observations than names given:\n",&contig[k*NAME_LEN]);
					fprintf(stdout,"found %d individuals in contig, and %d names in command line\n",ni,nfound);
					fprintf(stdout,"Individual(s) not given in command line: ");
					for (i=0; i<ni; i++){
						if (found[i]==0){
							fprintf(stdout,"%s ",name[i]);
						}
					}
					fprintf(stdout,"\n");
					fprintf(stdout,"Suppressing warnings about following contigs\n");
					nplus=1;
				}
			}
			j=0;
		}
		else { //line contains genotype data
				sscanf(line,"%d%[^\n]",&pos,tmpline);
				strcpy(line,tmpline);
				if(pos != j+1){
					fprintf(stderr,"Error in input file line number %d: expecting to read \"%d ...\"\n",l,j+1);
					fprintf(stderr,"read \"%s\" instead\n",line);
					exit(1);
				}
				if (j==0) {
					if((genotypes=(char *)malloc(sizeof(char)*(nfound)*2))==NULL){
						fprintf(stderr,"error in memory allocation\n");
						exit(1);
					}
				}
				else {
					if((genotypes=realloc(genotypes,sizeof(char)*(nfound)*2*pos))==NULL){
						fprintf(stderr,"error in memory allocation\n");
						exit(1);
					}
				}
				fi=0;
				for (i=0; i<ni; i++){
					if (found[i]>0) {
						sscanf(line,"\t%c%c|%*f%[^\n]",&genotypes[j*nfound*2+2*fi],&genotypes[j*nfound*2+2*fi+1],tmpline);
						fi++;
					}
					else {
						sscanf(line,"\t%*c%*c|%*f%[^\n]",tmpline);
					}
//					if(true_prob>2*GSL_DBL_MIN) {
//						totgenprob+=true_prob;
//						gensite++;
//					}
					strcpy(line,tmpline);
				}
				j++;
			
		}
//		memset(line, '\0', strlen(line)*(sizeof line));
	}
	//	fprintf(stdout,"%d sites (individuals * positions), error rate (reads2snp) = %f\n",gensite,1.-totgenprob/gensite);
	nJ=pos;
	
	//Filtering polymorphisms for the last contig	
	polysite[k]=polyfilter(nJ,nfound,genotypes,&npolysites[k],&n3[k],&n4[k],&theta[k]);
	free(genotypes);
	totsites+=npolysites[k];
	ncontigs=k+1;
	fprintf(stdout,"Read %d contigs and %d polymorphic sites\n",ncontigs,totsites);

	for (iname=0;iname<nnames;iname++) {
		if(found[iname]==0){
			fprintf(stdout,"Individual %s was not found in any of the contigs\n",name[iname]);
		}
	}

	//All reading has been done : clear memory.
	free(tmpline);
	free(line);
	fclose(fp);

	
	//outputting
	for (k=0;k<ncontigs;k++) {
		fprintf(outfile,">%s\t%d\t%d\t%d\t%f\n",&contig[k*NAME_LEN],npolysites[k],n3[k],n3[k],theta[k]);
		for (iname=0;iname<nnames;iname++) {
			for(i=0;i<ni;i++){
				if (strcmp(name[i],names[iname])==0) {
					nsites=0;
					nhetsites=0;
					for (t=0; t<npolysites[k]; t++){
						if(polysite[k][t][i+1]) {
							nsites++;
						}
						if(polysite[k][t][i+1]==N12) {
							nhetsites++;
						}
					}
					fprintf(outfile,"%s: %f\t",names[iname],(double)nhetsites/(double)nsites);
				}
			}
		}
		fprintf(outfile,"\n");
		if(npolysites[k]>0) {
			for (t=0; t<npolysites[k]; t++){
				for (i=0; i<nnames+4; i++) {
					fprintf(outfile,"%d\t",polysite[k][t][i]);
				}
				fprintf(outfile,"\n");
			}
		}
	}
	
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
	free(found);
	fclose(outfile);
	
	return 0;

}