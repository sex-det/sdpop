#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "reading.h"
#include "types.h"

int main(int argc, char *argv[]) 
{
	FILE *fp,*outfile;
	std::string line;
	char **names;
	int t,c,i,ii,l,j,nJ,pos,ipos,npos,lastpos,nind,k,chrpos,ni,firstcontig,ncontigs;
	int *sex,*het;
	double theta,a;
	int val,nnuc,maxnuc=5; //A, T, G, C, or N
	char nuc[maxnuc];
	char chrom[1000];
	ContigGenotypes contiggenotypes;
	int nA,nT,nC,nG,nN,div;
	int window=1000;
	
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
	
	k=0;
	l=0; //line number
	firstcontig=0;
	
	while ((c = std::fgetc(fp)) != EOF) { //loop through the file
		line="";
		while ( (c != '\n') && (c != EOF) ) { //loop through the line
			line+=char(c);
			c = std::fgetc(fp);
		}
		if (line.size() < 1) { //empty line : suppose it's the end of the file
			break;
		}
		//We've read one line :
		l++;
		if( strncmp(line.data(),"#CHROM",6)==0 ) { //only one line with individual names in vcf file				
			names=getvcfnames(line.data(),&ni);
			if((sex=(int *)calloc((size_t)(ni),sizeof(int)))==NULL) { 
				fprintf(stderr,"error in memory allocation\n");
				exit(1);
			}
			if((het=(int *)calloc((size_t)(ni),sizeof(int)))==NULL) { 
				fprintf(stderr,"error in memory allocation\n");
				exit(1);
			}
			for(i=0;i<ni;i++){
				sex[i]=0;
				Individual tmpindividual;
				tmpindividual.name=names[i];
				tmpindividual.sex=0;
				contiggenotypes.individuals.push_back(tmpindividual);
			}
		}
		else if (strncmp(line.data(),"##",2)!=0) { //line contains genotype data
			sscanf(line.data(),"%s\t%d",chrom,&chrpos);
			if(strcmp(chrom,contiggenotypes.name.data())==0){ //still the same contig as before
				
				val=vcfsnp(line.data(),nuc,&nnuc);
				if(val==-1){
					fprintf(stderr,"Error occurred while reading line %d\n",l);
					exit(1);
				}
				else if (val==1) {
					continue;
				}
				
				//here, we are sure it's a SNP. nnuc is the number of observed nucleotides
				Genotypes tempgenotypes;                         
				tempgenotypes.position=chrpos;
				tempgenotypes.individualgenotypes=vcfgenotypes(ni,sex,line.data(),nuc,maxnuc);
				contiggenotypes.genotypes.push_back(tempgenotypes);
				j++;				
			}
			else { //new contig
				//treat last read contig
				if (firstcontig!=1){
					fprintf(stdout,"Treating contig %s...\n",contiggenotypes.name.data());
					nJ=j;
					if(nJ>0){
						ipos=0;
						npos=contiggenotypes.genotypes.size();
						nind=contiggenotypes.individuals.size();
						pos=contiggenotypes.genotypes[ipos].position;
						lastpos=contiggenotypes.genotypes[npos-1].position;
						theta=0;
						for(t=1;t<=lastpos;t++){
							if(t % window == 0){ //output
								theta/=(double)window;
								fprintf(outfile,"%s\t%d\t%f",contiggenotypes.name.data(),t,theta);
								for (i=0; i<nind; i++){
									if (het[i]>0){
//										fprintf(outfile,"\t%f",(double)het[i]/double(window));
										fprintf(outfile,"\t%d",het[i]);
										het[i]=0;
									}
									else {
										fprintf(outfile,"\t%d",0);
									}
								}
								fprintf(outfile,"\n");
								theta=0;
							}
							if(t==pos){ //we have the genotypes
								nA=nT=nG=nC=nN=0;
								for (i=0; i<nind; i++){ 
									if (contiggenotypes.genotypes[ipos].individualgenotypes[i].nucleotides[0]!=contiggenotypes.genotypes[ipos].individualgenotypes[i].nucleotides[1]){
										het[i]++;
									}
									for (ii=0; ii<2; ii++){
										switch (contiggenotypes.genotypes[ipos].individualgenotypes[i].nucleotides[ii]) {
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
											fprintf(stderr,"Error in parsing genotypes at position %d of contig %s\n",contiggenotypes.genotypes[j].position,contiggenotypes.name.data());
											fprintf(stderr,"allowed are only A, T, G, C, and N (or n).\n");
											fprintf(stderr,"Found %c instead\n",contiggenotypes.genotypes[j].individualgenotypes[i].nucleotides[ii]);
											exit(1);
										}
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
								if (div > 1) {
									a=1.;
									for (i=2; i<2*nind; i++){
										a+=1./(double)i;
									}
									theta+=1./a;
								}
								
								ipos++;
								pos=contiggenotypes.genotypes[ipos].position;
							}
							else { //we don't have the genotypes; suppose the site is monoallelic
							}
						}
						
						
						contiggenotypes.genotypes.clear();
					}
					k++;
				}
				
				contiggenotypes.name=chrom;
				firstcontig=0;
				j=0;				
				
				//read first line
				
				val=vcfsnp(line.data(),nuc,&nnuc);
				if(val==-1){
					fprintf(stderr,"Error occurred while reading line %d\n",l);
					exit(52);
				}
				else if (val==1) {
					continue;
				}
				//here, we are sure it's a SNP. nnuc is the number of observed nucleotides
				
				Genotypes tempgenotypes;
				tempgenotypes.position=chrpos;
				tempgenotypes.individualgenotypes=vcfgenotypes(ni,sex,line.data(),nuc,maxnuc);
				contiggenotypes.genotypes.push_back(tempgenotypes);
				j++;
			}
		}
	}
	nJ=j;
	
	//Filtering polymorphisms for the last contig	
	if(nJ>0){
	}
	ncontigs=k+1;
	fprintf(stdout,"Read %d contigs\n",ncontigs);	

	for(i=0; i<ni; i++){
		free(names[i]);
	}
	free(names);
	free(sex);
	free(het);

	fclose(fp);
	fclose(outfile);
	
	return 0;
	
}
