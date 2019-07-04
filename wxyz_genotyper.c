#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

//1) read sdpop output and select genes considered as sex-linked
//2) read genotype file 
int DNA2int(char c) 
{
	switch (c) {
	case 'a' :
	case 'A' :
		return 1;
		break;
	case 'c' :
	case 'C' :
		return 2;
		break;
	case 'g' :
	case 'G' :
		return 3;
		break;
	case 't' :
	case 'T' :
		return 4;
		break;
	case '.' :
	case 'n' :
	case 'N' :
		return 0;
		break;
	default:
		fprintf(stderr,"DNA2int: error reading base. Allowed characters are a,t,g,c,A,T,G,C,.,n and N\n");
		exit(1);
		break;
	}
}

char int2DNA(int i) 
{
	switch (i) {
	case 1 :
		return 'A';
		break;
	case 2 :
		return 'C';
		break;
	case 3 :
		return 'G';
		break;
	case 4 :
		return 'T';
		break;
	default:
		return 'N';
		break;
	}
}

void outfunction(FILE *outfile, char *contig, int s,int t,int fixed,int ns,double pi_X,double pi_Y,double divergence,char *sequenceX,char *sequenceY){
	fprintf(outfile,">%s_X %d %d %d %f %f\n",contig,t,fixed,ns,pi_X/ns,divergence/ns);
	sequenceX[s]='\0';
	sequenceY[s]='\0';
	fprintf(outfile,"%s\n",sequenceX);
	fprintf(outfile,">%s_Y %d %d %d %f %f\n",contig,t,fixed,ns,pi_Y/ns,divergence/ns);
	fprintf(outfile,"%s\n",sequenceY);
	
	free(sequenceX);
	free(sequenceY);
}

int main(int argc, char *argv[]) 
{
	FILE *sdpfile,*genfile,*outfile;
	int CUR_MAX = 4095;
	int NAME_LEN = 9000;
	int NCONTIG_BATCH = 100;
	int ncontigs_allocated;
	char *line = calloc((size_t)CUR_MAX,sizeof(char));
	char *tmpline = calloc((size_t)CUR_MAX,sizeof(char));
	char *contig,***nuc,ch='a',gencontig[NAME_LEN],word[NAME_LEN];
	int *npolysites,toread=0,ncontigs,ngencontigs,sites_allocated,***polysite;
	double *contigpp,contigthreshold,sitethreshold1,sitethreshold2,***f;
	char *sequenceX,*sequenceY,nuc1,nuc2;
	int i,j,k,l,s,t,found,pos,nacgt[5],div,nextpos,x;
	int count,length,SNP,fixed,ns,ff[2],XY,nwords,posterior_field,mean;
	double pi_X,pi_Y,divergence;

	for(i=0;i<argc;i++) {
		fprintf(stdout,"%s ",argv[i]);
	}
	fprintf(stdout,"\n");
	
	if (argc != 10) {
		fprintf(stdout,"Usage: %s sdpopfile genfile outfile contigthreshold sitethreshold1 sitethreshold2 posterior_field system max/mean\n",argv[0]);
		exit(1);
	}
	
	if((sdpfile=fopen(argv[1],"r"))==NULL){
		fprintf(stderr,"error opening input file %s\n",argv[1]);
		exit(1);
	}
	if((genfile=fopen(argv[2],"r"))==NULL){
		fprintf(stderr,"error opening input file %s\n",argv[1]);
		exit(1);
	}
	if((outfile=fopen(argv[3],"w"))==NULL){
		fprintf(stderr,"error opening output file %s\n",argv[2]);
		exit(1);
	}	
	contigthreshold=atof(argv[4]);
	if(contigthreshold<0 || contigthreshold>1){
		fprintf(stderr,"Error: contigthreshold should be between 0 and 1; value given: %s\n",argv[4]);
		exit(1);
	}
	sitethreshold1=atof(argv[5]);
	if(sitethreshold1<0.5 || sitethreshold1>1){
		fprintf(stderr,"Error: sitethreshold1 should be between 0.5 and 1; value given: %s\n",argv[5]);
		exit(1);
	}
	sitethreshold2=atof(argv[6]);
	if(sitethreshold2<0.5 || sitethreshold2>sitethreshold1){
		fprintf(stderr,"Error: sitethreshold2 should be between 0.5 and sitethreshold1; value given: %s\n",argv[6]);
		exit(1);
	}
	if(strcmp(argv[8],"XY")==0){
		XY=1;
	}
	else if(strcmp(argv[8],"ZW")==0){
		XY=0;
	}
	else {
		fprintf(stderr,"Error: system should be either \"XY\" or \"ZW\" ; value given: \"%s\"\n",argv[8]);
		exit(1);
	}
	if(strcmp(argv[9],"mean")==0){
		mean=1;
	}
	else if(strcmp(argv[9],"max")==0){
		mean=0;
	}
	else {
		fprintf(stderr,"Error: choose either \"max\" or \"mean\" ; value given: \"%s\"\n",argv[9]);
		exit(1);
	}
		
		
	k=-1;
	if((contig=(char *)malloc(sizeof(char)*NCONTIG_BATCH*NAME_LEN))==NULL) { 
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}			
	if((npolysites=(int *)malloc(sizeof(int)*NCONTIG_BATCH))==NULL) {
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}
	if((contigpp=(double *)malloc(sizeof(double)*NCONTIG_BATCH))==NULL) {
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}
	if((polysite=(int ***)malloc(sizeof(int **)*NCONTIG_BATCH))==NULL) {
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}
	if((f=(double ***)malloc(sizeof(double **)*NCONTIG_BATCH))==NULL) {
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}
	if((nuc=(char ***)malloc(sizeof(char **)*NCONTIG_BATCH))==NULL) {
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}
	ncontigs_allocated=NCONTIG_BATCH;
	
	fprintf(stdout,"Reading data from file %s...\n",argv[1]);
	//sdpop file reading
	l=0; //line number
	ch='a';
	
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
			ch = getc(sdpfile); // read from stream.
			line[length] = ch;
			length++;
			count++;
		}
		line[length] = '\0';
//		if (length <= 1) { //empty line : suppose it's the end of the file
//			break;
//		}
		//We've read one line :
		l++;
		
		//parsing header
		if ( line[0] == '#' ){
			if ( line[1] == '>' ){
//				printf("%s",line);
				//finding contig posteriors
				nwords=0;
				posterior_field=-1;
				while ( sscanf(line,"%[^\t ]%*[\t ]%[^\n]",word,tmpline)==2)	{
					if(strcmp(word,argv[7])==0){
						posterior_field=nwords;
						break;
					}
					nwords++;
					strcpy(line,tmpline);
				}
//				printf("%d\n",posterior_field);
				if(posterior_field<0){
					printf("Error: word \"%s\" not found in line %d\n",argv[7],l);
					exit(1);
				}
			}
		}
		if ( line[0] == '#' ){
			if ( line[1] == 'p' ){
//				printf("%s",line);
				//finding allele frequency fields
				//X and Y (or Z and W) frequencies are separated by a space (not a tab)
				for(i=0;i<2;i++){
					ff[i]=-1;
				}
				nwords=0;
				while ( sscanf(line,"%[^\t ]%*[\t ]%[^\n]",word,tmpline)==2)	{
//					printf("%d \"%s\"\t",nwords,word);
					strcpy(line,tmpline);
					if(XY){
						if (mean==0){
							if(strcmp(word,"fx_max")==0) {
								ff[0]=nwords;
							}
							if(strcmp(word,"fy_max")==0) {
								ff[1]=nwords;
							}
						}
						else {
							if(strcmp(word,"fx_mean")==0) {
								ff[0]=nwords;
							}
							if(strcmp(word,"fy_mean")==0) {
								ff[1]=nwords;
							}
						}
					}
					else {
						if (mean==0){
							if(strcmp(word,"fz_max")==0) {
								ff[0]=nwords;
							}
							if(strcmp(word,"fw_max")==0) {
								ff[1]=nwords;
							}
						}
						else {
							if(strcmp(word,"fz_mean")==0) {
								ff[0]=nwords;
							}
							if(strcmp(word,"fw_mean")==0) {
								ff[1]=nwords;
							}
						}
					}
					nwords++;
				}			

				fprintf(outfile,"#contig_name N_poly_sites N_fixed_diff N_total pi divergence\n");
			}
		}
		
		if ( line[0] == '>') {
			k++;
			//Start preparing to read a new contig
			if ( k >= ncontigs_allocated ) {
				ncontigs_allocated+=NCONTIG_BATCH;
				if((contig=(char *)realloc(contig,sizeof(char)*ncontigs_allocated*NAME_LEN))==NULL) { 
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}
				if((npolysites=(int *)realloc(npolysites,sizeof(int)*ncontigs_allocated))==NULL) {
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}		
				if((contigpp=(double *)realloc(contigpp,sizeof(double)*ncontigs_allocated))==NULL) {
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}		
				if((polysite=(int ***)realloc(polysite,sizeof(int **)*ncontigs_allocated))==NULL) {
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}		
				if((f=(double ***)realloc(f,sizeof(double **)*ncontigs_allocated))==NULL) {
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}		
				if((nuc=(char ***)realloc(nuc,sizeof(char **)*ncontigs_allocated))==NULL) {
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}		
			}
			if ( strlen(line) >= NAME_LEN ) {
				fprintf(stderr,"Error: a line name is too long (more than %d characters), line %d.\n",NAME_LEN-1,l);
				exit(1);
			}
//			printf("%s",line);
			sscanf(line,">%s\t%d\t%[^\n]",&contig[k*NAME_LEN],&npolysites[k],tmpline);
//			sscanf(line,">%s\t%d\t%*f\t%*d\t%*f\t%*f\t%*d\t%*f\t%*f\t%*d\t%*f\t%*f\t%*d\t%*f\t%*f\t%*d\t%*f\t%lf\t%*d\n",
//				&contig[k*NAME_LEN],&npolysites[k],&contigpp[k]);
			nwords=2;
			strcpy(line,tmpline);
			while ( sscanf(line,"%[^\t ]%*[\t ]%[^\n]",word,tmpline)==2)	{
				strcpy(line,tmpline);
				if(nwords==posterior_field){
					contigpp[k]=atof(word);
//					printf("%d %s %f\n",nwords,word,contigpp[k]);
					break;
				}
				nwords++;
			}	

			if(contigpp[k]>=contigthreshold){
//				printf("%f\n",contigpp[k]);
				toread=1;
				if((polysite[k]=(int **)malloc(sizeof(int *)*npolysites[k]))==NULL) {
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}		
				if((f[k]=(double **)malloc(sizeof(double *)*npolysites[k]))==NULL) {
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}		
				if((nuc[k]=(char **)malloc(sizeof(char *)*npolysites[k]))==NULL) {
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}		
				for(t=0;t<npolysites[k];t++){
					if((polysite[k][t]=(int *)malloc(sizeof(int)*7))==NULL) {
						fprintf(stderr,"error in memory allocation\n");
						exit(1);
					}		
					if((f[k][t]=(double *)malloc(sizeof(double)*2))==NULL) {
						fprintf(stderr,"error in memory allocation\n");
						exit(1);
					}		
					if((nuc[k][t]=(char *)malloc(sizeof(char)*2))==NULL) {
						fprintf(stderr,"error in memory allocation\n");
						exit(1);
					}		
					
				}
				t=0;
			}
			else {
				k--;
				toread=0;
			}

		}
		else if(isdigit(line[0]) && toread==1){
			if(t>=npolysites[k]){
				fprintf(stderr,"Error: contig %s does seem to have more polymorphic sites than announced\n",&contig[k*NAME_LEN]);
				exit(1);				
			}
//			printf("%s",line);
			sscanf(line,"%d\t%c%c\t%d\t%d\t%d\t%d\t%d\t%d\t%[^\n]",
				&polysite[k][t][0],&nuc[k][t][0],&nuc[k][t][1],&polysite[k][t][1],&polysite[k][t][2],&polysite[k][t][3],&polysite[k][t][4],&polysite[k][t][5],&polysite[k][t][6],tmpline);
			strcpy(line,tmpline);
			nwords=8;
			j=0;
			while ( sscanf(line,"%[^\t ]%*[\t ]%[^\n]",word,tmpline)==2)	{
				strcpy(line,tmpline);
				if(nwords==ff[j]){
					f[k][t][j]=atof(word);
					j++;
					if(j>1){
						break;
					}
				}
				nwords++;
			}
			t++;
		}
	}
	ncontigs=k+1;

	fprintf(stdout,"%d contigs selected.\n",ncontigs);
	fclose(sdpfile);

	if(ncontigs>0){
	fprintf(stdout,"Reading data from file %s...\n",argv[2]);
	ch='a';
	l=0;
	found=-1;
	k=0;
	ngencontigs=0;
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
			ch = getc(genfile); // read from stream.
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
		if ( line[0] == '>'){
			if(found>=0){ // print last sequences
				outfunction(outfile,&contig[k*NAME_LEN],s,t,fixed,ns,pi_X,pi_Y,divergence,sequenceX,sequenceY);
			}
			sscanf(line,">%s",gencontig);
			found=-1;
			for(k=0;k<ncontigs;k++){
				if(strcmp(gencontig,&contig[k*NAME_LEN])==0){
					found=k;
					t=0;
					s=0;
					ns=0;
					fixed=0;
					ngencontigs+=1;
					sites_allocated=2*polysite[k][(npolysites[k]-1)][0];
					if((sequenceX=(char *)malloc(sizeof(char)*(sites_allocated+1)))==NULL){
						fprintf(stderr,"error in memory allocation\n");
						exit(1);
					}				
					if((sequenceY=(char *)malloc(sizeof(char)*(sites_allocated+1)))==NULL){
						fprintf(stderr,"error in memory allocation\n");
						exit(1);
					}
					pi_X=0;
					pi_Y=0;
					divergence=0;
					break;
				}
			}
		}
		else if ( line[0] != 'p' && found>=0 ){
			if(s==sites_allocated){ //time to expand
				sites_allocated+=npolysites[found];
				if((sequenceX=realloc(sequenceX,sizeof(char)*(sites_allocated+1)))==NULL){
					fprintf(stderr,"error in memory allocation (realloc)\n");
					exit(1);
				}
				if((sequenceY=realloc(sequenceY,sizeof(char)*(sites_allocated+1)))==NULL){
					fprintf(stderr,"error in memory allocation (realloc)\n");
					exit(1);
				}				
			}
			sscanf(line,"%d%[^\n]",&pos,tmpline);
			strcpy(line,tmpline);
//			if(pos!=s+1){
//				fprintf(stderr,"Error in input file line number %d: expecting to read \"%d ...\"\n",l,s+1);
//				fprintf(stderr,"read \"%s\" instead\n",line);
//				exit(1);
//			}
			if(t<npolysites[k]) {
				if(pos==polysite[k][t][0]) { //a SNP, and we know the genotype
					SNP=1;
				}
				else {
					SNP=0;
				}
			}
			else {
				SNP=0;
			}
			
			if(SNP==1){
				ns++;
				if(f[k][t][0]>=sitethreshold1){
					sequenceX[s]=toupper(nuc[k][t][0]);
				}
				else if(f[k][t][0]>=sitethreshold2){
					sequenceX[s]=tolower(nuc[k][t][0]);
				}
				else if(f[k][t][0]<=1-sitethreshold1){
					sequenceX[s]=toupper(nuc[k][t][1]);
				}
				else if(f[k][t][0]<=1-sitethreshold2){
					sequenceX[s]=tolower(nuc[k][t][1]);
				}
				else {
					sequenceX[s]='n';
				}
				if(f[k][t][1]>=sitethreshold1){
					sequenceY[s]=toupper(nuc[k][t][0]);
				}
				else if(f[k][t][1]>=sitethreshold2){
					sequenceY[s]=tolower(nuc[k][t][0]);
				}
				else if(f[k][t][1]<=1-sitethreshold1){
					sequenceY[s]=toupper(nuc[k][t][1]);
				}
				else if(f[k][t][1]<=1-sitethreshold2){
					sequenceY[s]=tolower(nuc[k][t][1]);
				}
				else {
					sequenceY[s]='n';
				}
				pi_X+=2*f[k][t][0]*(1.-f[k][t][0]);
				pi_Y+=2*f[k][t][1]*(1.-f[k][t][1]);
				divergence+=f[k][t][0]*(1.-f[k][t][1])+f[k][t][1]*(1.-f[k][t][0]);
//				if( (f[k][t][0]>=sitethreshold1 || f[k][t][0]<=1-sitethreshold1) && (f[k][t][1]>=sitethreshold1 || f[k][t][1]<=1-sitethreshold1) ){
				if( (f[k][t][0]>=sitethreshold1 && f[k][t][1]<=1-sitethreshold1) || (f[k][t][1]>=sitethreshold1 && f[k][t][0]<=1-sitethreshold1) ){
					fixed++;
				}
				t++;
			}
			else{ //not a SNP : find nucleotide in genotype file
				for(i=0;i<5;i++){
					nacgt[i]=0;
				}
				j=0;
//				printf("%s\n",line);
				nextpos=0;
				while(strlen(line)-nextpos>3){
					sscanf(line + nextpos,"\t%c%c|%*f%n",&nuc1,&nuc2,&x);
 //					strcpy(line,tmpline);
 //					printf("%d\t%s\t%c%c\n",nextpos,&line[nextpos],nuc1,nuc2);
 					nacgt[DNA2int(nuc1)]++;
 					nacgt[DNA2int(nuc2)]++;
					j++;
					nextpos+=x;
				}
				if(nacgt[0]==2*j){ //N
					sequenceX[s]='N';
					sequenceY[s]='N';
				}
				else {
					ns++;
					div=0;
					for(i=1;i<5;i++){
						if(nacgt[i]>0){
							div++;
						}
					}
					if(div==1){
						for(i=1;i<5;i++){
							if(nacgt[i]>0){
								sequenceX[s]=int2DNA(i);
								sequenceY[s]=int2DNA(i);
								break;
							}
						}
					}
					else { //if there is some polymorphism that has passed unnoticed...
						sequenceX[s]='X';
						sequenceY[s]='X';
					}
				}
			}
			s++;
		}
	}
	if(found>=0){ // print last sequences
		outfunction(outfile,&contig[k*NAME_LEN],s,t,fixed,ns,pi_X,pi_Y,divergence,sequenceX,sequenceY);
	}
	}
	
	fclose(genfile);
	if(ngencontigs!=ncontigs){
		fprintf(stderr,"Warning: %d contigs selected from %s, but %d found in %s\n",ncontigs,argv[1],ngencontigs,argv[2]);
	}
	else {
		fprintf(stdout,"%d contigs outputted.\n",ngencontigs);
	}
	fclose(outfile);
	
	
	return 0;
}