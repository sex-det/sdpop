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

int **zygosity(int npos, int nind, char *genotypes, int *nsites, int *n2, int *n3, int *n4, double *theta) {
	int i,j,tempi;
	int div,nA,nT,nG,nC,nN;
	int **zyg=NULL;
	char pos1,pos2;
	int n,ndata;
	double a;
	
	*n2=0;
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
		if (div > 0 ) { //data
			if (n==0) {
				if((zyg=(int **)malloc(sizeof(int *)))==NULL){
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}
				if((zyg[n]=(int *)malloc(sizeof(int)*(nind+1)))==NULL){
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}
			}
			else {
				if((zyg=realloc(zyg,sizeof(int *)*(n+1)))==NULL){
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}
				if((zyg[n]=(int *)malloc(sizeof(int)*(nind+1)))==NULL){
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}
			}
			zyg[n][0]=j+1; //position
			for (i=0; i<nind; i++){
				zyg[n][i+1]=0;				
			}
			
			ndata=0;
			for (i=0; i<nind; i++){ //loop through all chromosomes
				pos1=genotypes[j*nind*2+2*i];
				pos2=genotypes[j*nind*2+2*i+1];
				if (pos1 != 'n' && pos1 != 'N'){
					if (pos1 == pos2) { //homozygous
						zyg[n][i+1]=1;
					}
					else if (pos1 != pos2) { //heterozygous
						zyg[n][i+1]=2;
					}
					ndata++;
				}
			}
			//Watterson's estimate of theta
			if(div>1){
				a=1.;
				for (i=2; i<2*(ndata); i++){
					a+=1./(double)i;
				}
				*theta+=1./a;
			}
			if (div==2){
				(*n2)++;
			}
			else if (div ==3) {
				(*n3)++;
			}
			else if (div ==4) {
				(*n4)++;
			}
			n++;
		}
	}

	*theta/=(double)n;
	*nsites=n;
	if(n==0) {
		return NULL;
	}
	else {
		return zyg;
	}
}

int main(int argc, char *argv[]) 
{
	FILE *fp,*outfile;
	int CUR_MAX = 4095;
	int NAME_LEN = 500;
	int NCONTIG_BATCH = 100;
	int ncontigs_allocated;
	char *line = calloc((size_t)CUR_MAX,sizeof(char));
	char *tmpline = calloc((size_t)CUR_MAX,sizeof(char));
	char *contig;
	int ***zyg;
	int *nsites,*n2,*n3,*n4;
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
	int nhom,nhet,**nhetsites,**nhomsites,all=0;

	for(i=0;i<argc;i++) {
		fprintf(stdout,"%s ",argv[i]);
	}
	fprintf(stdout,"\n");
	
	if (argc > 4 || argc < 3) {
		fprintf(stdout,"Usage: %s infile outfile (namelist)\n",argv[0]);
		fprintf(stdout,"Without a list of names, all individuals are taken into account to calculate statistics\n");
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
	
	if (argc == 4) {
		names=parsenames(argv[3],&nnames,NAME_LEN);
	}
	else {
		all=1;
	}
	
	//genotype file reading
	l=0; //line number
	ch='a';
	if((nsites=(int *)malloc(sizeof(int)*NCONTIG_BATCH))==NULL) {
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}
	if((zyg=(int ***)malloc(sizeof(int **)*NCONTIG_BATCH))==NULL) {
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}
	if((n2=(int *)malloc(sizeof(int)*NCONTIG_BATCH))==NULL) {
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
	if((contig=(char *)malloc(sizeof(char)*NCONTIG_BATCH*NAME_LEN))==NULL) { 
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}
	ncontigs_allocated=NCONTIG_BATCH;
	
	firstcontig=1;
	k=0;
	l=1;
	
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
				zyg[k]=zygosity(nJ,nnames,genotypes,&nsites[k],&n2[k],&n3[k],&n4[k],&theta[k]);
				free(genotypes);
				totsites+=nsites[k];
				k++;
				if(k % 5000 == 0){
					fprintf(stdout,"%d contigs, %d sites, and still reading...\n",k,totsites);
				}
			}
			
			//Start preparing to read a new contig
			if ( k >= ncontigs_allocated ) {
				ncontigs_allocated+=NCONTIG_BATCH;
				if((contig=(char *)realloc(contig,sizeof(char)*ncontigs_allocated*NAME_LEN))==NULL) { 
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}	
				if((nsites=(int *)realloc(nsites,sizeof(int)*ncontigs_allocated))==NULL) {
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}		
				if((n2=(int *)realloc(n2,sizeof(int)*ncontigs_allocated))==NULL) {
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
				if((zyg=(int ***)realloc(zyg,sizeof(int **)*ncontigs_allocated))==NULL) {
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
			if (!firstcontig) {
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
				found[i]=-1;
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
				if(! (all && firstcontig) ){
					for(iname=0;iname<nnames;iname++) {
						if (strcmp(name[i],names[iname])==0) {
							found[i]=iname;
							nfound++;
							break;
						}
					}
				}
			}
			if(firstcontig) {
				firstcontig=0;
				if(all){ // Here, I suppose that all individuals are present in the first contig
					nnames=ni;
					if((names=(char **)calloc((size_t)nnames,sizeof(char *)))==NULL) { 
						fprintf(stderr,"error in memory allocation\n");
						exit(1);
					}
					for (i=0; i<nnames; i++){
						if((names[i]= (char *)malloc(sizeof(char) * NAME_LEN))==NULL) { 
							fprintf(stderr,"error in memory allocation\n");
							exit(1);
						}
					}
					for(i=0;i<nnames;i++) {
						strcpy(names[i],name[i]);						
						found[i]=i;
						nfound++;
					}						
				}
			}

			
			if(ni>nfound){
				if(nplus==0) {
					fprintf(stdout,"Contig %s seems to have more observations than names given:\n",&contig[k*NAME_LEN]);
					fprintf(stdout,"found %d individuals in contig, of which %d in command line\n",ni,nfound);
					fprintf(stdout,"Individual(s) not given in command line: ");
					for (i=0; i<ni; i++){
						if (found[i]==-1){
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
					if((genotypes=(char *)malloc(sizeof(char)*(nnames)*2))==NULL){
						fprintf(stderr,"error in memory allocation\n");
						exit(1);
					}
				}
				else {
					if((genotypes=(char *)realloc(genotypes,sizeof(char)*(nnames)*2*pos))==NULL){
						fprintf(stderr,"error in memory allocation\n");
						exit(1);
					}
				}
				for(iname=0;iname<nnames;iname++) {
					genotypes[j*nnames*2+2*iname]='N';
					genotypes[j*nnames*2+2*iname+1]='N';					
				}
				fi=0;
				for (i=0; i<ni; i++){
					if (found[i]>=0) {
						sscanf(line,"\t%c%c|%*f%[^\n]",&genotypes[j*nnames*2+2*found[i]],&genotypes[j*nnames*2+2*found[i]+1],tmpline);
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
	
	//Filtering polymorphisms for the last contig	
	nJ=pos;
	zyg[k]=zygosity(nJ,nnames,genotypes,&nsites[k],&n2[k],&n3[k],&n4[k],&theta[k]);
	free(genotypes);
	totsites+=nsites[k];
	ncontigs=k+1;
	fprintf(stdout,"Read %d contigs and %d sites\n",ncontigs,totsites);

	//All reading has been done : clear memory.
	free(tmpline);
	free(line);
	fclose(fp);

	if((nhomsites=(int **)malloc(sizeof(int*)*ncontigs))==NULL) {
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}
	if((nhetsites=(int **)malloc(sizeof(int*)*ncontigs))==NULL) {
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}
	for (k=0;k<ncontigs;k++) {
		if((nhomsites[k]=(int *)malloc(sizeof(int)*nnames))==NULL) {
			fprintf(stderr,"error in memory allocation\n");
			exit(1);
		}
		if((nhetsites[k]=(int *)malloc(sizeof(int)*nnames))==NULL) {
			fprintf(stderr,"error in memory allocation\n");
			exit(1);
		}
	}

	//counting total number of homozygous and heterozygous sites per individual
	for (k=0;k<ncontigs;k++) {
		if(nsites[k]>0) {
			for (iname=0;iname<nnames;iname++) {
				nhomsites[k][iname]=0;
				nhetsites[k][iname]=0;
				for (t=0; t<nsites[k]; t++){
					if(zyg[k][t][iname+1]==1) {
						nhomsites[k][iname]++;
					}
					if(zyg[k][t][iname+1]==2) {
						nhetsites[k][iname]++;
					}
				}
			}
		}
	}
	
	//outputting
	fprintf(outfile,"Summary per individual:\n");
	for (iname=0;iname<nnames;iname++) {
		nhom=0;
		nhet=0;
		for (k=0;k<ncontigs;k++) {
			if(nsites[k]>0){
			nhom+=nhomsites[k][iname];
			nhet+=nhetsites[k][iname];
			}
		}				
		fprintf(outfile,"%s: %d, %d, %f\n",names[iname],nhom,nhet,(double)nhet/(double)(nhom+nhet));
	}			

	fprintf(outfile,"\n");
	fprintf(outfile,"Contig-wise statistics:\n");

	for (k=0;k<ncontigs;k++) {
		fprintf(outfile,">%s\t%d\t%d\t%d\t%d\t%f\n",&contig[k*NAME_LEN],nsites[k],n2[k],n3[k],n4[k],theta[k]);
		if(nsites[k]>0) {
			for (iname=0;iname<nnames;iname++) {
				if(nhomsites[k][iname]+nhetsites[k][iname]>0){
					fprintf(outfile,"%s: %f\t",names[iname],(double)nhetsites[k][iname]/(double)(nhomsites[k][iname]+nhetsites[k][iname]));
				}
			}
			fprintf(outfile,"\n");
		}
	}
	
	
	free(contig);
	for(k=0;k<ncontigs;k++) {
		if(nsites[k]>0) {
			for(t=0; t<nsites[k]; t++){
				free(zyg[k][t]);
			}
			free(zyg[k]);
		}
	}
	free(zyg);
	free(nsites);
	free(found);
	fclose(outfile);
	
	return 0;

}