#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

int main(int argc, char *argv[]) {
	
	int namelen=1000;
	char tcont[namelen],word[namelen];
	char line[100000],tmpline[100000];
	int c,l,j,i,f[7],nwords,length,nsites,n;
	double L[7],logL[7],sumL,sumlogL,m;
	FILE *fp;

	if((fp=fopen(argv[1],"r"))==NULL){
		fprintf(stderr,"error opening input file %s\n",argv[1]);
		exit(1);
	}

	i=0;
	l=0; //line number
	c=0;
		
	while ((c = fgetc(fp)) != EOF) { //loop through the file
		length=0;
		while ( (c != '\n') && (c != EOF) ) { //loop through the line
			line[length]=c;
			c = fgetc(fp);
			length++;
		}
		
		line[length] = '\0';
		if (length <= 1) { //empty line : suppose it's the end of the file
			break;
		}
		//We've read one line :
		l++;
//		printf("%s\n",line);

		if ( line[0] == '#' ){
			if ( line[1] == 'p' ){
				for(j=0;j<7;j++){
					f[j]=-1;
				}
				nwords=0;
				while ( sscanf(line,"%[^\t ]%*[\t ]%[^\n]",word,tmpline)==2)	{
//					printf("%d \"%s\"\t",nwords,word);
					strcpy(line,tmpline);
					if(strcmp(word,"logL_autosomal")==0) {
						f[0]=nwords;
					}
					if(strcmp(word,"logL_haploid")==0) {
						f[1]=nwords;
					}
					if(strcmp(word,"logL_paralog")==0) {
						f[2]=nwords;
					}
					if(strcmp(word,"logL_xhemizygote")==0) {
						f[3]=nwords;
					}
					if(strcmp(word,"logL_xy")==0) {
						f[4]=nwords;     
					}                    
					if(strcmp(word,"logL_zhemizygote")==0) {
						f[5]=nwords;     
					}                    
					if(strcmp(word,"logL_zw")==0) {
						f[6]=nwords;     
					}
					nwords++;
				}
				printf("contig_name\tN_sites\tmean_coverage");
				printf("\tmean_autosomal");
				printf(" geom_autosomal");
				printf("\tmean_haploid");
				printf(" geom_haploid");
				printf("\tmean_paralog");
				printf(" geom_paralog");
				printf("\tmean_xhemizygote");
				printf(" geom_xhemizygote");
				printf("\tmean_xy");
				printf(" geom_xy");
				printf("\tmean_zhemizygote");
				printf(" geom_zhemizygote");
				printf("\tmean_zw");
				printf(" geom_zw");
				printf("\n");
				
			}
		}
		else if ( line[0] == '>' ){
			if (i > 0){
				printf("%s\t%d\t%f",tcont,n,m);
				sumL=0.;
				sumlogL=0.;
				for(j=0;j<7;j++){
					if(f[j]>-1){
						logL[j]=exp(logL[j]/nsites);
						sumlogL+=logL[j];
						sumL+=L[j];
					}
				}
				for(j=0;j<7;j++){
					printf("\t%e %e",L[j]/sumL,logL[j]/sumlogL);
				}
				printf("\n");
			}
			sscanf(line,">%s\t%d\t%lf%*[^\n]",tcont,&n,&m);
			i++;
			for(j=0;j<7;j++){
				L[j]=0;
				logL[j]=0;
				nsites=0;
			}
		}
		else if (line[0] != '#' && i > 0){
			nsites++;
			nwords=0;
			j=0;
			while ( sscanf(line,"%[^\t ]%*[\t ]%[^\n]",word,tmpline)==2)	{
				strcpy(line,tmpline);
				if(nwords==f[j]){
					L[j]+=exp(atof(word));
					logL[j]+=atof(word);
					j++;
					while (f[j]<0){
						j++;
					}
				}
				nwords++;
			}
		}
	}
	if (i > 0){
		printf("%s\t%d\t%f",tcont,n,m);
		sumL=0.;
		sumlogL=0.;
		for(j=0;j<7;j++){
			if(f[j]>-1){
				logL[j]=exp(logL[j]/nsites);
				sumlogL+=logL[j];
				sumL+=L[j];
			}
		}
		for(j=0;j<7;j++){
			printf("\t%e %e",L[j]/sumL,logL[j]/sumlogL);
		}
		printf("\n");
	}
	fclose(fp);			
	return 0;
}

