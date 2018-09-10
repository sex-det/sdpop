#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

int main(int argc, char *argv[]) 
{
	FILE *fp,*falr,*outfile;
	int CUR_MAX = 4095;
	int NAME_LEN = 100;
	int NCONTIG_BATCH = 1000;
	int NSITES_BATCH = 1000;
	int ncontigs_allocated,nsites_allocated;
	char *line = calloc((size_t)CUR_MAX,sizeof(char));
	char *tmpline = calloc((size_t)CUR_MAX,sizeof(char));
	char *contig, *alrcontig= calloc((size_t)NAME_LEN,sizeof(char));;
	int *nsites,*individuals;
	float **genprob;
	int ncontigs,totsites=0;
	int count = 0; 
	int length = 0;
	int i,j,k,l,ni,firstcontig,pos,t;
	char ch,c;
	char **names,**genotypes;
	int nJ; //number of individuals
	clock_t begin = clock();

	for(i=0;i<argc;i++) {
		fprintf(stdout,"%s ",argv[i]);
	}
	fprintf(stdout,"\n");
	
	if (argc != 4) {
		fprintf(stdout,"Usage: %s genfile-in order-file genfile-out\n",argv[0]);
		exit(1);
	}
	
	if((fp=fopen(argv[1],"r"))==NULL){
		fprintf(stderr,"error opening input file %s\n",argv[1]);
		exit(1);
	}
	if((falr=fopen(argv[2],"r"))==NULL){
		fprintf(stderr,"error opening input file %s\n",argv[2]);
		exit(1);
	}
	if((outfile=fopen(argv[3],"w"))==NULL){
		fprintf(stderr,"error opening output file %s\n",argv[3]);
		exit(1);
	}
	//genotype file reading
	l=0; //line number
	ch='a';
	if((individuals=(int *)malloc(sizeof(int **)*NCONTIG_BATCH))==NULL) {
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}		
	if((nsites=(int *)malloc(sizeof(int)*NCONTIG_BATCH))==NULL) {
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}
	if((contig=(char *)malloc(sizeof(char)*NCONTIG_BATCH*NAME_LEN))==NULL) { 
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}
	if((names=(char **)malloc(sizeof(char *)*NCONTIG_BATCH))==NULL) {
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}		
	if((genotypes=(char **)malloc(sizeof(char *)*NCONTIG_BATCH))==NULL) {
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}
	if((genprob=(float **)malloc(sizeof(float *)*NCONTIG_BATCH))==NULL) {
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
			if(count == CUR_MAX) { // time to expand (for unexpectedly large line lengths) ?
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
			if (firstcontig==0) {
				nJ=pos;
				nsites[k]=nJ;
				totsites+=nJ;
				k++;
				if(k % 5000 == 0){
					fprintf(stdout,"%d contigs, %d sites, and still reading...\n",k,totsites);
				}
			}
			else {
				firstcontig=0;
			}
			
			//Start preparing to read a new contig
			if ( k >= ncontigs_allocated ) {
				ncontigs_allocated+=NCONTIG_BATCH;
				if((contig=(char *)realloc(contig,sizeof(char)*ncontigs_allocated*NAME_LEN))==NULL) { 
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}	
				if((nsites=(int *)realloc(nsites,sizeof(int **)*ncontigs_allocated))==NULL) {
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}		
				if((individuals=(int *)realloc(individuals,sizeof(int **)*ncontigs_allocated))==NULL) {
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}		
				if((names=(char **)realloc(names,sizeof(char *)*ncontigs_allocated))==NULL) {
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}		
				if((genotypes=(char **)realloc(genotypes,sizeof(char *)*ncontigs_allocated))==NULL) {
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}		
				if((genprob=(float **)realloc(genprob,sizeof(float *)*ncontigs_allocated))==NULL) {
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
			ni=0; //we'll have "position\tname1\tname2\t...\tnameN\0 : the number of tabs is the number of individuals.
			i=0;
			while ((c=line[i]) != '\0'){
				if (c == '\t') 
					ni++;
				i++;
			}
			individuals[k]=ni;	
			if((names[k]=(char *)malloc(sizeof(char)*CUR_MAX))==NULL) {
				fprintf(stderr,"error in memory allocation\n");
				exit(1);
			}		
			strcpy(names[k],line);
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
					if((genotypes[k]=(char *)malloc(sizeof(char)*ni*2*NSITES_BATCH))==NULL){
						fprintf(stderr,"error in memory allocation\n");
						exit(1);
					}
					if((genprob[k]=(float *)malloc(sizeof(float)*ni*NSITES_BATCH))==NULL){
						fprintf(stderr,"error in memory allocation\n");
						exit(1);
					}
					nsites_allocated=NSITES_BATCH;
				}
				else if (j==nsites_allocated) {
					nsites_allocated+=NSITES_BATCH;
					if((genotypes[k]=(char *)realloc(genotypes[k],sizeof(char)*ni*2*nsites_allocated))==NULL){
						fprintf(stderr,"error in memory allocation\n");
						exit(1);
					}
					if((genprob[k]=(float *)realloc(genprob[k],sizeof(float)*ni*nsites_allocated))==NULL){
						fprintf(stderr,"error in memory allocation\n");
						exit(1);
					}
				}
				for (i=0; i<ni; i++){
					sscanf(line,"\t%c%c|%f%[^\n]",&genotypes[k][j*ni*2+2*i],&genotypes[k][j*ni*2+2*i+1],&genprob[k][j*ni+i],tmpline);
					strcpy(line,tmpline);
				}
				j++;
			
		}
//		memset(line, '\0', strlen(line)*(sizeof line));
	}
	nJ=pos;
	nsites[k]=nJ;
	totsites+=nJ;
	ncontigs=k+1;
	fprintf(stdout,"Read %d contigs and %d sites\n",ncontigs,totsites);
	//All reading has been done
	fclose(fp);

	clock_t end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	fprintf(stdout,"Reading took %f seconds\n",time_spent);
	
	//outputting in the right order
	ch='a';
	l=0;
	while (ch != EOF) { //loop through the alr file
		ch='a';
		count = 0;
		length = 0;
		while ( (ch != '\n') && (ch != EOF) ) { //loop through the line
			if(count == CUR_MAX) { // time to expand (for unexpectedly large line lengths) ?
				CUR_MAX *= 2; 
				count = 0;
				line = realloc(line, sizeof(char) * CUR_MAX); 
				tmpline = realloc(tmpline, sizeof(char) * CUR_MAX); 
			}
			ch = getc(falr); // read from stream.
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
			sscanf(line,">%s",alrcontig);
			for(k=0;k<ncontigs;k++) { //search could be optimized, e.g.: shorten the list for each found contig
				if(strcmp(&contig[k*NAME_LEN],alrcontig)==0){
					fprintf(outfile,">%s\n",&contig[k*NAME_LEN]);
					fprintf(outfile,"%s",names[k]);
					for(t=0;t<nsites[k];t++){
						fprintf(outfile,"%d\t",t+1);
						for(i=0;i<individuals[k];i++){
							if(genprob[k][t*individuals[k]+i]>0.9999999999999){
								fprintf(outfile,"%c%c|1",genotypes[k][t*individuals[k]*2+2*i],genotypes[k][t*individuals[k]*2+2*i+1]);
							}
							else if (genprob[k][t*individuals[k]+i]<0.0000000000001){
								fprintf(outfile,"%c%c|0",genotypes[k][t*individuals[k]*2+2*i],genotypes[k][t*individuals[k]*2+2*i+1]);
							}
							else {
								fprintf(outfile,"%c%c|%f",genotypes[k][t*individuals[k]*2+2*i],genotypes[k][t*individuals[k]*2+2*i+1],genprob[k][t*individuals[k]+i]);
							}
							if(i!=individuals[k]-1){
								fprintf(outfile,"\t");
							}
						}
						fprintf(outfile,"\n");
					}
					break;
				}
			}
		}
	}

	/*		
	while(fscanf(falr,">%s\n",alrcontig)!=EOF) {
		for(k=0;k<ncontigs;k++) { //search could be optimized, e.g.: shorten the list for each found contig
			if(strcmp(&contig[k*NAME_LEN],alrcontig)==0){
				fprintf(outfile,">%s\n",&contig[k*NAME_LEN]);
				fprintf(outfile,"%s",names[k]);
				for(t=0;t<nsites[k];t++){
					fprintf(outfile,"%d\t",t+1);
					for(i=0;i<individuals[k];i++){
						if(genprob[k][t*individuals[k]+i]>0.9999999999999){
							fprintf(outfile,"%c%c|1",genotypes[k][t*individuals[k]*2+2*i],genotypes[k][t*individuals[k]*2+2*i+1]);
						}
						else if (genprob[k][t*individuals[k]+i]<0.0000000000001){
							fprintf(outfile,"%c%c|0",genotypes[k][t*individuals[k]*2+2*i],genotypes[k][t*individuals[k]*2+2*i+1]);
						}
						else {
							fprintf(outfile,"%c%c|%f",genotypes[k][t*individuals[k]*2+2*i],genotypes[k][t*individuals[k]*2+2*i+1],genprob[k][t*individuals[k]+i]);
						}
						if(i!=individuals[k]-1){
							fprintf(outfile,"\t");
						}
					}
					fprintf(outfile,"\n");
				}
				break;
			}
		}
	}
	*/
	fclose(falr);
	fclose(outfile);
	free(tmpline);
	free(line);
	
	for(k=0;k<ncontigs;k++) {
		free(names[k]);
		free(genotypes[k]);
		free(genprob[k]);
	}
	free(names);
	free(genotypes);
	free(genprob);
	free(alrcontig);
	free(contig);
	free(nsites);
	free(individuals);
	
	end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	fprintf(stdout,"Total running time was %f seconds\n",time_spent);
	
	return 0;
}