#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "reading.h"
#include "types.h"

struct GenotypeInformation {
	int heterozygoussites;
	int homozygoussites;
	int nongenotypessites;
};

void countandoutput(FILE *outfile,ContigGenotypes contiggenotypes,int window)
{
	int t,i,ii,pos,ipos,npos,lastpos,nind,nobs;
	double theta,a,pi,pis;
	char nuc1, nuc2;
	int nA,nT,nC,nG,nN,div,missingsites=0,nsites=0;
	std::vector<GenotypeInformation> genotypeinformation;
	static int times=0;
	
	nind=contiggenotypes.individuals.size();
	if(times==0){
		fprintf(outfile,"Contig\tlastsite\tn_sites\tmissing_sites\ttheta\tpi");
		for (i=0; i<nind; i++){
			fprintf(outfile,"\t%s_het",contiggenotypes.individuals[i].name.data());
			fprintf(outfile,"\t%s_hom",contiggenotypes.individuals[i].name.data());
			fprintf(outfile,"\t%s_N",contiggenotypes.individuals[i].name.data());
		}
		fprintf(outfile,"\n");		
	}
	npos=contiggenotypes.genotypes.size();
	if(npos>0){
		for (i=0; i<nind; i++){
			GenotypeInformation tmpgenotypeinfo;
			tmpgenotypeinfo.heterozygoussites=0;
			tmpgenotypeinfo.homozygoussites=0;
			tmpgenotypeinfo.nongenotypessites=0;
			genotypeinformation.push_back(tmpgenotypeinfo);
		}
		ipos=0;
		pos=contiggenotypes.genotypes[ipos].position;
		lastpos=contiggenotypes.genotypes[npos-1].position;
		theta=0.;
		pi=0.;
		for(t=1;t<=lastpos;t++){
			if(t==pos){ //we have the genotypes
				nA=nT=nG=nC=nN=0;
				nobs=0;
				for (i=0; i<nind; i++){
					nuc1=contiggenotypes.genotypes[ipos].individualgenotypes[i].nucleotides[0];
					nuc2=contiggenotypes.genotypes[ipos].individualgenotypes[i].nucleotides[1];
					if (nuc1!=nuc2){
						genotypeinformation[i].heterozygoussites++;
						nobs++;
					}
					else if (nuc1=='N' || nuc1=='n'){
						genotypeinformation[i].nongenotypessites++;
					}
					else {
						genotypeinformation[i].homozygoussites++;
						nobs++;
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
							fprintf(stderr,"Error in parsing genotypes at position %d of contig %s\n",contiggenotypes.genotypes[ipos].position,contiggenotypes.name.data());
							fprintf(stderr,"allowed are only A, T, G, C, and N (or n).\n");
							fprintf(stderr,"Found %c instead\n",contiggenotypes.genotypes[ipos].individualgenotypes[i].nucleotides[ii]);
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
					for (i=2; i<2*nobs; i++){
						a+=1./(double)i;
					}
					theta+=1./a;
					pis=nA*nT;
					pis+=nA*nC;
					pis+=nA*nG;
					pis+=nT*nC;
					pis+=nT*nG;
					pis+=nG*nC;
					pis*=(2./(double)(2*nobs*(2*nobs-1)));
					pi+=pis;
				}
				
				ipos++;
				pos=contiggenotypes.genotypes[ipos].position;
			}
			else { //we don't have the genotypes
				missingsites++;
			}
			nsites++;
			if(t==lastpos || (window > 0 && t % window == 0)){ //output
				//								theta/=(double)window;
				fprintf(outfile,"%s\t%d\t%d\t%d\t%f\t%f",contiggenotypes.name.data(),t,nsites,missingsites,theta,pi);
				for (i=0; i<nind; i++){
					fprintf(outfile,"\t%d",genotypeinformation[i].heterozygoussites);
					fprintf(outfile,"\t%d",genotypeinformation[i].homozygoussites);
					fprintf(outfile,"\t%d",genotypeinformation[i].nongenotypessites);
					genotypeinformation[i].heterozygoussites=0;
					genotypeinformation[i].homozygoussites=0;
					genotypeinformation[i].nongenotypessites=0;
				}
				fprintf(outfile,"\n");
				theta=0;
				pi=0;
				missingsites=0;
				nsites=0;
			}							
		}
	}
	times++;
}


int main(int argc, char *argv[]) 
{
	FILE *fp,*outfile;
	std::string line;
	char **names;
	int c,i,l,chrpos,ni,firstcontig;
	int *sex;
	int val,nnuc,maxnuc=5; //A, T, G, C, or N
	char nuc[maxnuc];
	char chrom[1000];
	ContigGenotypes contiggenotypes;
	int window=1000;
	
	for(i=0;i<argc;i++) {
		fprintf(stdout,"%s ",argv[i]);
	}
	fprintf(stdout,"\n");
	
	fprintf(stdout,"This program gives nucleotide diversity and heterozygosity counts for all the samples in a vcf file.\n");
	fprintf(stdout,"Note that theta a pi values are given as the sums of sites in the window.\n");
	fprintf(stdout,"It's up to the user to decide if these should be divided by the window size, or only by the sites for which genotypes were available.\n");
		
	if (argc != 4) {
		fprintf(stdout,"Usage: %s infile outfile windowsize\n",argv[0]);
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
	window=atoi(argv[3]);
	if(window==0){
		printf("calculating statistics per contig\n");
	}
	else {
		printf("window size: %d\n",window);
	}

	l=0; //line number
	firstcontig=1;
	
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
			}
			else { //new contig
				//treat last read contig
				if (firstcontig!=1){
					fprintf(stdout,"Treating contig %s...\n",contiggenotypes.name.data());
					countandoutput(outfile,contiggenotypes,window);
					contiggenotypes.genotypes.clear();
				}
				
				contiggenotypes.name=chrom;
				firstcontig=0;
				
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
			}
		}
	}
	
	//Filtering polymorphisms for the last contig	
	fprintf(stdout,"Treating contig %s...\n",contiggenotypes.name.data());
	countandoutput(outfile,contiggenotypes,window);

	for(i=0; i<ni; i++){
		free(names[i]);
	}
	free(names);
	free(sex);

	fclose(fp);
	fclose(outfile);
	
	return 0;
	
}
