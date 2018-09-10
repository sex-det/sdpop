#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define FEMALE 0
#define MALE 1
#define SEXES 2
#define N11 1
#define N12 2
#define N22 3
#define N11F 1
#define N12F 2
#define N22F 3
#define N11M 4
#define N12M 5
#define N22M 6
#define READS2SNP 1
#define GBS_FILTERED 2

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

int **genfilter(int npos, int nind, char *genotypes, int *sex, int *npolysites, int *n3, int *n4, double *theta) {//number of positions and individuals,
	// genotype matrix, pointer to the number of polymorphic sites
	// function treats one contig
	char nucvec[4] = {0};
	int i,j,randbit,tempi;
	int div,n11f,n12f,n22f,n11m,n12m,n22m,nA,nT,nG,nC,nN;
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
		if (div > 1) {
			a=1.;
			for (i=2; i<2*nind; i++){
				a+=1./(double)i;
			}
			*theta+=1./a;
		}
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
//			printf("%d\t%d\t%d\t%d\t%d\t%d\n",n11f,n12f,n22f,n11m,n12m,n22m);
//			if(n22m+n12m == 0 && n22f+n12f == 0){ //are the polymorphisms in the individuals we study ?
//				div=-1;
//			}				
//			if(n11m+n12m == 0 && n11f+n12f == 0){ 
//				div=-1;
//			}
			if(div == 2) {
				n2++;
			}
			if(n11m+n22m+n12m == 0 || n11f+n22f+n12f == 0){ //if we only have data on one sex, don't consider the site
				div=-1;
			}				
		}
		else if (div ==3) {
			*n3++;
		}
		else if (div ==4) {
			*n4++;
		}
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
			polysite[n][N11F]=n11f; //female counts
			polysite[n][N12F]=n12f; 
			polysite[n][N22F]=n22f; 			
			polysite[n][N11M]=n11m; //male counts
			polysite[n][N12M]=n12m; 
			polysite[n][N22M]=n22m; 	
			n++;
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
	int NAME_LEN = 9000;
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
	int pipepos,ichar,nfem,nmal,ifem,imal,nplus=0;
	int *sex,*foundsex;
	char ch,c;
	char **name,**femname,**malname,***sitename;
	int *ffound,*mfound;
	int nI=0,nJ,nf,nm,nfound; //number of individuals
	char *genotypes;
	int randomise=0,ri,tempsex,datafmt;

	for(i=0;i<argc;i++) {
		fprintf(stdout,"%s ",argv[i]);
	}
	fprintf(stdout,"\n");
	
	if (argc < 6 || argc > 7) {
		fprintf(stdout,"Usage: %s infile outfile type namelist1 namelist2 (r)\n",argv[0]);
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
	
	if(strcmp(argv[3],"r")==0 || strcmp(argv[3],"1")==0) {
		datafmt=READS2SNP;
	}
	else if (strcmp(argv[3],"g")==0 || strcmp(argv[3],"2")==0) {
		datafmt=GBS_FILTERED;
	}
	else {
		fprintf(stdout,"Usage: %s infile outfile type namelist1 namelist2 (r)\n",argv[0]);
		fprintf(stderr,"type should be either \"r\" or \"1\" for reads2snp genotype file, or \"g\" or \"2\" for filtered GBS data.\n");
		exit(1);	
	}
	
	femname=parsenames(argv[4],&nfem,NAME_LEN);
	malname=parsenames(argv[5],&nmal,NAME_LEN);	
	if (argc == 7) {
		if (strcmp(argv[6],"r")==0) {
			fprintf(stdout,"Randomisation chosen\n");
			randomise=1;
		}
		else {
			fprintf(stderr,"6th argument should be \"r\" if you want randomisation, or be nothing at all\n");
			exit(1);
		}
	}		
	
	if((ffound=(int *)calloc((size_t)(nfem),sizeof(int)))==NULL) {
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}
	if((mfound=(int *)calloc((size_t)(nmal),sizeof(int)))==NULL) {
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}
	for (ifem=0;ifem<nfem;ifem++) {
		ffound[ifem]=0;
	}
	for (imal=0;imal<nmal;imal++) {
		mfound[imal]=0;
	}

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
	if ( datafmt==GBS_FILTERED ) {
		if((sitename=(char ***)malloc(sizeof(char **)*NCONTIG_BATCH))==NULL) {
			fprintf(stderr,"error in memory allocation\n");
			exit(1);
		}
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
		
		if ( (datafmt==READS2SNP && line[0] == '>') || (datafmt==GBS_FILTERED && l==1)){
			if (firstcontig==0) {	//Filtering polymorphisms for the last read contig
				nJ=pos;
				if(randomise){ //Fisher-Yates shuffling
					for (i=0;i<nfound-1;i++) {
						ri=i+rand()/(RAND_MAX/(nfound-i)+1);
						tempsex=foundsex[i];
						foundsex[i]=foundsex[ri];
						foundsex[ri]=tempsex;
					}
//					for (i=0;i<nfound;i++) {
//						fprintf(stdout,"%d\t",foundsex[i]);
//					}
//					fprintf(stdout,"\n");
				}
				polysite[k]=polyfilter(nJ,nfound,genotypes,foundsex,&npolysites[k],&n3[k],&n4[k],&theta[k]);
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
				if ( datafmt==GBS_FILTERED ) {
					if((sitename=(char ***)realloc(sitename,sizeof(char **)*ncontigs_allocated))==NULL) {
						fprintf(stderr,"error in memory allocation\n");
						exit(1);
					}
				}		
				if((polysite=(int ***)realloc(polysite,sizeof(int **)*ncontigs_allocated))==NULL) {
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}		
				if((npolysites=(int *)realloc(npolysites,sizeof(int **)*ncontigs_allocated))==NULL) {
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}		
				if((n3=(int *)realloc(n3,sizeof(int **)*ncontigs_allocated))==NULL) {
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}		
				if((n4=(int *)realloc(n4,sizeof(int **)*ncontigs_allocated))==NULL) {
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}		
				if((theta=(double *)realloc(theta,sizeof(double **)*ncontigs_allocated))==NULL) {
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}		
			}
			if ( strlen(line) >= NAME_LEN ) {
				fprintf(stderr,"Error: a contig name is too long (more than %d characters) on line %d.\n",NAME_LEN-1,l);
				exit(1);
			}
			if ( datafmt==READS2SNP ) {
				sscanf(line,">%s",&contig[k*NAME_LEN]);
				//			printf("%s\n",&contig[k*NAME_LEN]);
			}
			else if ( datafmt==GBS_FILTERED ) { //only one "contig", and names are in the first line
				sprintf(&contig[k*NAME_LEN],"GBS_data");				
				//we'll have to count the number of individuals first : 11 fields separated by tabs, then the names of individuals
				ni=0; //The number of tabs-10 is the number of individuals.
				i=0;
				while ((c=line[i]) != '\0'){
					if (c == '\t') 
						ni++;
					i++;
				}
				ni-=10;
				if((name=(char **)calloc((size_t)(ni),sizeof(char *)))==NULL) { 
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}
				//find names iteratively ; shorten the line for each name found :
				for (i=0;i<11;i++) { //first remove the non-names
					sscanf(line,"%*s\t%[^\n]",tmpline);
					strcpy(line,tmpline);
				}
				for (i=0; i<ni; i++){
					if((name[i]= (char *)malloc(sizeof(char) * NAME_LEN))==NULL) { 
						fprintf(stderr,"error in memory allocation\n");
						exit(1);
					}
					sscanf(line,"%s%[^\n]",name[i],tmpline);
					strcpy(line,tmpline);
					fprintf(stdout,"%s\t",name[i]);
				}
				//Find out the sex of the individuals.
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
//				if(nm+nf<ni){
//					fprintf(stdout,"No sex found for individual(s): ");
//					for (i=0; i<ni; i++){
//						if ( sex[i]==-1 ) {
//							fprintf(stdout,"%s ",name[i]);
//						}		
//					}
//					fprintf(stdout,"\n");
//				}
				firstcontig=0;
			}
		}
		else if ( datafmt==READS2SNP && line[0] == 'p' ){
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
				if((name=(char **)calloc((size_t)(ni),sizeof(char *)))==NULL) { 
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
//				if(nm+nf<ni){
//					fprintf(stdout,"No sex found for individual(s): ");
//					for (i=0; i<ni; i++){
//						if ( sex[i]==-1 ) {
//							fprintf(stdout,"%s ",name[i]);
//						}		
//					}
//					fprintf(stdout,"\n");
//				}
				firstcontig=0;
			}
			else { //if not the first contig
				//check if the number of individuals is the same
				ni=0; //we'll have "position\tname1\tname2\t...\tnameN\0 : the number of tabs is the number of individuals.
				i=0;
				for (ifem=0;ifem<nfem;ifem++) {
					ffound[ifem]=0;
				}
				for (imal=0;imal<nmal;imal++) {
					mfound[imal]=0;
				}
				while ((c=line[i]) != '\0' && line[i] != '\n'){
					if (c == '\t') 
						ni++;
					i++;
				}
				if(ni>nfem+nmal){
					if(nplus==0) {
					fprintf(stdout,"Contig %s seems to have more observations than names given\n",&contig[k*NAME_LEN]);
					fprintf(stdout,"Suppressing warnings about following contigs\n");
					nplus=1;
					}
				}
				//Find out the sex of the individuals.
				//sex has the sex of all individuals in the *gen file (sex can be male, female, or absent)
				//foundsex has the sex of the individuals we want to study, in the order of the *gen file
				//genotypes only has the genotypes of the individuals we want to study, in the order of the *gen file
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
//				if(nm+nf<ni){
//					fprintf(stdout,"No sex found for individual(s): ");
//					for (i=0; i<ni; i++){
//						if ( sex[i]==-1 ) {
//							fprintf(stdout,"%s ",name[i]);
//						}		
//					}
//					fprintf(stdout,"\n");
//				}
			}
			j=0;
		}
		else { //line contains genotype data
			if(datafmt==READS2SNP){
				sscanf(line,"%d%[^\n]",&pos,tmpline);
				strcpy(line,tmpline);
				if(pos != j+1){
					fprintf(stderr,"Error in input file line number %d: expecting to read \"%d ...\"\n",l,j+1);
					fprintf(stderr,"read \"%s\" instead\n",line);
					exit(1);
				}
			}
			else if (datafmt==GBS_FILTERED){
				if(j==0){
					if((sitename[k]=(char **)malloc(sizeof(char *)))==NULL){
						fprintf(stderr,"error in memory allocation\n");
						exit(1);
					}
				} 
				else {
					if((sitename[k]=realloc(sitename[k],sizeof(char *)*(j+1)))==NULL){
						fprintf(stderr,"error in memory allocation\n");
						exit(1);
					}
				}
				if((sitename[k][j]=(char *)malloc(sizeof(char)*NAME_LEN))==NULL){
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}
				sscanf(line,"%s%[^\n]",sitename[k][j],tmpline);
				strcpy(line,tmpline);
				for (i=0;i<10;i++) { //remove the non-used fields
					sscanf(line,"\t%*s%[^\n]",tmpline);
					strcpy(line,tmpline);
				}
			}
			if (j==0) {
				//					if((genotypes=(char *)malloc(sizeof(char)*ni*2))==NULL){
				if((genotypes=(char *)malloc(sizeof(char)*(nfound)*2))==NULL){
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}
			}
			else {
				//					if((genotypes=realloc(genotypes,sizeof(char)*ni*2*pos))==NULL){
				//					if((genotypes=realloc(genotypes,sizeof(char)*(nfound)*2*pos))==NULL){
				if((genotypes=realloc(genotypes,sizeof(char)*(nfound)*2*(j+1)))==NULL){
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}
			}
				fi=0;
				for (i=0; i<ni; i++){
//					sscanf(line,"\t%c%c|%lf%[^\n]",&genotypes[j*ni*2+2*i],&genotypes[j*ni*2+2*i+1],&true_prob,tmpline);
//					sscanf(line,"\t%c%c|%*f%[^\n]",&genotypes[j*ni*2+2*i],&genotypes[j*ni*2+2*i+1],tmpline);
					if (sex[i]>=0) {
//						sscanf(line,"\t%c%c|%*f%[^\n]",&genotypes[j*nfound*2+2*fi],&genotypes[j*nfound*2+2*fi+1],tmpline);
						sscanf(line,"\t%c|%c%[^\n]",&genotypes[j*nfound*2+2*fi],&genotypes[j*nfound*2+2*fi+1],tmpline);
//						fprintf(stdout,"%c%c\t",genotypes[j*nfound*2+2*fi],genotypes[j*nfound*2+2*fi+1]);
						fi++;
					}
					else {
//						sscanf(line,"\t%*c%*c|%*f%[^\n]",tmpline);
						sscanf(line,"\t%*c|%*c%[^\n]",tmpline);
					}
//					if(true_prob>2*GSL_DBL_MIN) {
//						totgenprob+=true_prob;
//						gensite++;
//					}
					strcpy(line,tmpline);
				}
//					fprintf(stdout,"\n");
				j++;
			
		}
//		memset(line, '\0', strlen(line)*(sizeof line));
	}
	//	fprintf(stdout,"%d sites (individuals * positions), error rate (reads2snp) = %f\n",gensite,1.-totgenprob/gensite);
	nJ=j;
	
	//Filtering polymorphisms for the last contig	
	polysite[k]=polyfilter(nJ,nfound,genotypes,foundsex,&npolysites[k],&n3[k],&n4[k],&theta[k]);
	free(genotypes);
	totsites+=npolysites[k];
	ncontigs=k+1;
	fprintf(stdout,"Read %d contigs and %d polymorphic sites\n",ncontigs,totsites);

	for (ifem=0;ifem<nfem;ifem++) {
		if(ffound[ifem]==0){
			fprintf(stdout,"Individual %s was not found in any of the contigs\n",femname[ifem]);
		}
	}
	for (imal=0;imal<nmal;imal++) {
		if (mfound[imal]==0) {
			fprintf(stdout,"Individual %s was not found in any of the contigs\n",malname[imal]);
		}
	}

	//All reading has been done : clear memory.
	free(tmpline);
	free(line);
	fclose(fp);

	
	//outputting
	for (k=0;k<ncontigs;k++) {
		fprintf(outfile,">%s\t%d\t%d\t%d\t%f\n",&contig[k*NAME_LEN],npolysites[k],n3[k],n3[k],theta[k]);
		if(npolysites[k]>0) {
			for (t=0; t<npolysites[k]; t++){
				if (datafmt==GBS_FILTERED){
					fprintf(outfile,"%s\t",sitename[k][polysite[k][t][0]-1]);
					for (i=1; i<7; i++) {
						fprintf(outfile,"%d\t",polysite[k][t][i]);
					}
				}
				else if (datafmt==READS2SNP){
					for (i=0; i<7; i++) {
						fprintf(outfile,"%d\t",polysite[k][t][i]);
					}
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
	free(ffound);
	free(mfound);
	fclose(outfile);
	
	return 0;

}