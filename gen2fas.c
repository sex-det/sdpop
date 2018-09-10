#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "reading.h"

extern int NAME_LEN;

int main(int argc, char *argv[]) 
{
	FILE *fp,*outfile;
	int CUR_MAX = 4095;
	char *line = calloc((size_t)CUR_MAX,sizeof(char));
	char *tmpline = calloc((size_t)CUR_MAX,sizeof(char));
	char *contig;
	int i,j,l,ni,fi,firstcontig,pos,a;
	char ch;
	char **name;
	int nJ,count,length;
	char *genotypes;

	for(i=0;i<argc;i++) {
		fprintf(stdout,"%s ",argv[i]);
	}
	fprintf(stdout,"\n");

	if((fp=fopen(argv[1],"r"))==NULL){
		fprintf(stderr,"error opening input file %s\n",argv[1]);
		exit(1);
	}
	if((outfile=fopen(argv[2],"w"))==NULL){
		fprintf(stderr,"error opening output file %s\n",argv[2]);
		exit(1);
	}	
	
	//genotype file reading
	l=0; //line number
	ch='a';
		
	firstcontig=1;
	
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
				for (i=0; i<ni; i++){
					for(a=0;a<2;a++){
						fprintf(outfile,">%s_%s_Allele%d\n",contig,name[i],a+1);
						for(j=0;j<nJ;j++){
							fprintf(outfile,"%c",genotypes[j*ni*2+2*i+a]);
							if((j + 1)  % 80 == 0){
								fprintf(outfile,"\n");
							}
						}
						fprintf(outfile,"\n");
					}
				}
				free(genotypes);
				for (i=0; i < ni; i++){
					free(name[i]);
				}
				free(name);
				
				free(contig);
			}
			firstcontig=0;
			
			//Start preparing to read a new contig
			if((contig=(char *)malloc(sizeof(char)*NAME_LEN))==NULL) { 
				fprintf(stderr,"error in memory allocation\n");
				exit(15);
			}	
			
			if ( strlen(line) >= NAME_LEN ) {
				fprintf(stderr,"Error: a contig name is too long (more than %d characters) on line %d.\n",NAME_LEN-1,l);
				exit(22);
			}
			sscanf(line,">%s",contig);
			
		}
		else if ( line[0] == 'p' ){
			name=getgennames(line, &ni); //returns all names present in line, and the number of names ni
			if(ni<=0){
				fprintf(stdout, "No individual names found in line %d; format error ?\nexiting\n",l);
				exit(1);
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
				if((genotypes=(char *)malloc(sizeof(char)*(ni)*2))==NULL){
					fprintf(stderr,"error in memory allocation\n");
					exit(34);
				}
			}
			else {
				if((genotypes=realloc(genotypes,sizeof(char)*(ni)*2*(j+1)))==NULL){
					fprintf(stderr,"error in memory allocation\n");
					exit(35);
				}
			}
			fi=0;
			for (i=0; i<ni; i++){
				sscanf(line,"\t%c%c|%*f%[^\n]",&genotypes[j*ni*2+2*fi],&genotypes[j*ni*2+2*fi+1],tmpline);
				fi++;
				strcpy(line,tmpline);
			}
			j++;
		}
	}
	//	fprintf(stdout,"%d sites (individuals * positions), error rate (reads2snp) = %f\n",gensite,1.-totgenprob/gensite);
	nJ=j;
	for (i=0; i<ni; i++){
		for(a=0;a<2;a++){
			fprintf(outfile,">%s_%s_Allele%d\n",contig,name[i],a+1);
			for(j=0;j<nJ;j++){
				fprintf(outfile,"%c",genotypes[j*ni*2+2*i+a]);
				if((j +1) % 80 == 0){
					fprintf(outfile,"\n");
				}
			}
			fprintf(outfile,"\n");
		}
	}
	free(genotypes);
	for (i=0; i < ni; i++){
		free(name[i]);
	}
	free(name);	
	free(contig);

	return 0;
}