#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
//#include "reading.h"

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
#define VCF 3

int NAME_LEN = 9000;

void writegenheader(FILE *outfile, char *contig, char **name, int *conserve, int ni, int nnames)
{
	int fi,i;
	fprintf(outfile,">%s\n",contig);			
	fprintf(outfile,"position\t");
	fi=0;
	for (i=0; i<ni; i++){
		if (conserve[i]>0) {
			fprintf(outfile,"sp|%s",name[i]);
			fi++;
			if(fi<nnames){
				fprintf(outfile,"\t");
			}
			else {
				fprintf(outfile,"\n");
			}	
		}
	}
}


int vcfsnp(char *line,char *nuc, int *nnuc)
{
	char *tmpline = calloc((size_t)strlen(line),sizeof(char));
	char allele1[NAME_LEN],allele2[NAME_LEN];
	int maxnuc=5; //A, T, G, C, or N
	int i,len;
	char c;
	int mnuc;
	
	//find out if this is a SNP, and what the observed alleles are
	sscanf(line,"\t%*s\t%s\t%s%[^\n]",allele1,allele2,tmpline);
	strcpy(line,tmpline);
//	fprintf(stdout,"%s\t%s\n",allele1,allele2);
	if(strlen(allele1)>1){
		return(1); //not a SNP
	}	
//	fprintf(stdout,"%s\t%s\t%d\t%d\t",allele1,allele2,(int)strlen(allele1),(int)strlen(allele2));
	nuc[0]=allele1[0];
	//allele2 might be a comma-separated list
	if((len=(int)strlen(allele2))>1){
		i=0;
		mnuc=1;
		while(i<=len-1){
			c=allele2[i];
			if(i==len-1 || allele2[i+1]==','){
				if(mnuc>maxnuc){
					fprintf(stderr,"Error: too many different nucleotides given\n");
					free(tmpline);
					return(-1);
				}
				nuc[mnuc]=c;
				mnuc++;
			}
			else {
				free(tmpline);
				return(1); //not a SNP
			}
			i=i+2;
		}
	}
	else { //strlen(allele2)==1
		nuc[1]=allele2[0];
		mnuc=2;
	}
//	fprintf(stdout,"--%d--",mnuc);
	*nnuc=mnuc;
	free(tmpline);
	return(0);
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
		exit(56);
	}
	for (i=0; i<n; i++){
		if((names[i]= (char *)malloc(sizeof(char) * maxlength))==NULL) { 
			fprintf(stderr,"error in memory allocation\n");
			exit(57);
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

int main(int argc, char *argv[]) 
{
	FILE *fp,*outfile;
	int CUR_MAX = 4095;
	char *line = calloc((size_t)CUR_MAX,sizeof(char));
	char *tmpline = calloc((size_t)CUR_MAX,sizeof(char));
	char *contig;
	int count = 0; 
	int length = 0;
	int i,ii,j,k,l,ni,fi,firstcontig;
	int iind,nplus=0;
	char ch,c,oldc;
	char **name,**names;
	int *found;
	int nfound,nnames,*conserve; //number of individuals
	char genotypes[2];
	int val,nnuc,maxnuc=5; //A, T, G, C, or N
	char nuc[maxnuc];
	char chrom[NAME_LEN],chrpos[NAME_LEN],vcfformatstring[NAME_LEN],genotypestring[NAME_LEN],contigpos[NAME_LEN];
	int all=0,maxsize=0;

	for(i=0;i<argc;i++) {
		fprintf(stdout,"%s ",argv[i]);
	}
	fprintf(stdout,"\n");	
	
	if (argc < 4 || argc > 5) {
		fprintf(stdout,"Usage: %s input_vcf_file output_gen_file [a|namelist] (maxsize)\n",argv[0]);
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
	if(strcmp(argv[3],"a")==0){
		all=1;
		fprintf(stdout,"All samples will be parsed\n");
	}
	else {
		names=parsenames(argv[3],&nnames,NAME_LEN);
	}
	if(argc==5){
		maxsize=atoi(argv[4]);
		fprintf(stdout,"Maximum size of contigs in output file: %d\n",maxsize);
	}
	
	if((found=(int *)calloc((size_t)(nnames),sizeof(int)))==NULL) {
		fprintf(stderr,"error in memory allocation\n");
		exit(6);
	}	
	for (iind=0;iind<nnames;iind++) {
		found[iind]=0;
	}

	if((contig=(char *)malloc(sizeof(char)*NAME_LEN))==NULL) { 
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}	
	
	//vcf file reading
	l=0; //line number
	ch='a';
		
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
		
		if( strncmp(line,"#CHROM",6)==0 ) { //only one line with individual names in vcf file				
			ni=0;
			i=0;
			while ((c=line[i]) != '\0'){
				//					fprintf(stdout,"%c,",c);
				if ( c == '\t' || c == ' ' ){
					if ( oldc == '\t' || oldc == ' '){
						i++;							
						continue;
					}						
					ni++;
					//						fprintf(stdout,"ni=%d\n",ni);
				}
				if(ni<=8){
					strcpy(tmpline,&line[1]);
					strcpy(line,tmpline);
				}
				else {
					i++;
				}
				oldc=c;
			}
			ni=ni-8; //first 8 tabs separate information fields
			fprintf(stdout,"Treating header; %d names found...\n",ni);
			if((name=(char **)calloc((size_t)(ni),sizeof(char *)))==NULL) { 
				fprintf(stderr,"error in memory allocation\n");
				exit(36);
			}
			//find names iteratively ; shorten the line for each name found :
			//				sscanf(line,"position%[^\n]",tmpline);
			//				strcpy(line,tmpline);
			for (i=0; i<ni; i++){
				if((name[i]= (char *)malloc(sizeof(char) * NAME_LEN))==NULL) { 
					fprintf(stderr,"error in memory allocation\n");
					exit(37);
				}
				sscanf(line,"\t%s%[^\n]",name[i],tmpline);
				strcpy(line,tmpline);
				//					fprintf(stdout,"%s\t",name[i]);
			}
			if((conserve=(int *)calloc((size_t)(ni),sizeof(int)))==NULL) { 
				fprintf(stderr,"error in memory allocation\n");
				exit(38);
			}
			if(all==1){
				for (i=0; i<ni; i++){			
					conserve[i]=1;
				}
				nnames=ni;
			}
			else {
				nfound=0;
				for (i=0; i<ni; i++){			
					conserve[i]=-1;
					//fprintf(stdout,"%s",name[i]);
					for (iind=0;iind<nnames;iind++) {
						if (strcmp(name[i],names[iind])==0) {
							//fprintf(stdout,": female\t");
							conserve[i]=1;
							found[iind]=1;
							nfound++;
						}
					}
				}
				if(ni>nnames){
					if(nplus==0) {
						fprintf(stdout,"There seem to be more observations than names given\n");
						fprintf(stdout,"found %d individuals in file, and %d names in command line\n",ni,nnames);
						fprintf(stdout,"Individual(s) not found: ");
						for (i=0; i<ni; i++){
							if ( conserve[i]==-1 ) {
								fprintf(stdout,"%s ",name[i]);
							}		
						}
						fprintf(stdout,"\n");
						nplus=1;
					}
				}
			}
			
		}
		else if (strncmp(line,"##",2)!=0) { //line contains genotype data
			sscanf(line,"%s\t%s%[^\n]",chrom,chrpos,tmpline);
			strcpy(line,tmpline);
			if(strcmp(chrom,contig)==0){ //still the same contig as before
				
				val=vcfsnp(line,nuc,&nnuc);
				if(val==-1){
					fprintf(stderr,"Error occurred while reading line %d\n",l);
					exit(1);
				}
				else if (val==1) {
					continue;
				}
				
				if(j==0){
					if(maxsize==0){
						writegenheader(outfile,contig,name,conserve,ni,nnames);
					}
					else {
						sprintf(contigpos,"%s_%s",contig,chrpos);
						writegenheader(outfile,contigpos,name,conserve,ni,nnames);
					}
				}
				else if(maxsize > 0 && j>=maxsize){
					j=0;
					sprintf(contigpos,"%s_%s",contig,chrpos);
					writegenheader(outfile,contigpos,name,conserve,ni,nnames);
				}
				
				fprintf(outfile,"%s\t",chrpos);
				//strip fields we don't use 
				sscanf(line,"\t%*s\t%*s\t%*s\t%s%[^\n]",vcfformatstring,tmpline);
				strcpy(line,tmpline);
				if(strncmp(vcfformatstring,"GT",2)!=0){
					fprintf(stderr,"Error: this doesn't seem to be a supported format (\"GT\" tag not found where expected), line %d\n",l);
					exit(1);					
				}
				
				fi=0;
				for (i=0; i<ni; i++){
					if (conserve[i]>0) {
						sscanf(line,"\t%s%[^\n]",genotypestring,tmpline);
						for(ii=0;ii<2;ii++){
							c=genotypestring[2*ii];
							if(c=='.'){
								genotypes[ii]='N';
								if(ii==0){
									genotypestring[2]='.';
								}
							}
							else if ( isdigit(c) && (c - '0') <= maxnuc ) {
								genotypes[ii]=nuc[c - '0'];
							}
							else {
								fprintf(stderr,"Error in reading genotype, line %d\n",l);
								exit(1);
							}
						}
						fprintf(outfile,"%c%c|1",genotypes[0],genotypes[1]);
						fi++;
						if(fi<nnames){
							fprintf(outfile,"\t");
						}
						else {
							fprintf(outfile,"\n");
						}	
					}
					else {
						sscanf(line,"\t%*s%[^\n]",tmpline);
					}
					strcpy(line,tmpline);
				}
				j++;
			}
			else { //new contig
				//treat last read contig
				if (firstcontig!=1){
					fprintf(stdout,"Treating contig %s...\n",contig);
					free(contig);
					k++;
				}
				
				if((contig=(char *)malloc(sizeof(char)*NAME_LEN))==NULL) { 
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}	
				strcpy(contig,chrom);
				firstcontig=0;
				j=0;				
								
				//read first line
				
				val=vcfsnp(line,nuc,&nnuc);
				if(val==-1){
					fprintf(stderr,"Error occurred while reading line %d\n",l);
					exit(1);
				}
				else if (val==1) {
					continue;
				}
				if(maxsize==0){
					writegenheader(outfile,contig,name,conserve,ni,nnames);
				}
				else {
					sprintf(contigpos,"%s_%s",contig,chrpos);
					writegenheader(outfile,contigpos,name,conserve,ni,nnames);
				}

				fprintf(outfile,"%s\t",chrpos);
				//strip fields we don't use 
				sscanf(line,"\t%*s\t%*s\t%*s\t%s%[^\n]",vcfformatstring,tmpline);
				strcpy(line,tmpline);
				if(strncmp(vcfformatstring,"GT",2)!=0){
					fprintf(stderr,"Error: this doesn't seem to be a supported format (\"GT\" tag not found where expected), line %d\n",l);
					exit(1);					
				}
				
				fi=0;
				for (i=0; i<ni; i++){
					if (conserve[i]>0) {
						sscanf(line,"\t%s%[^\n]",genotypestring,tmpline);
						for(ii=0;ii<2;ii++){
							c=genotypestring[2*ii];
							if(c=='.'){
								genotypes[ii]='N';
								if(ii==0){
									genotypestring[2]='.';
								}
							}
							else if ( isdigit(c) && (c - '0') <= maxnuc ) {
								genotypes[ii]=nuc[c - '0'];
							}
							else {
								fprintf(stderr,"Error in reading genotype, line %d\n",l);
								exit(1);
							}
						}
						fprintf(outfile,"%c%c|1",genotypes[0],genotypes[1]);
						fi++;
						if(fi<nnames){
							fprintf(outfile,"\t");
						}
						else {
							fprintf(outfile,"\n");
						}	
					}
					else {
						sscanf(line,"\t%*s%[^\n]",tmpline);
					}
					strcpy(line,tmpline);
				}
				j++;
			}
		}
	}

	fclose(fp);
	fclose(outfile);
	return 0;
	}
	
