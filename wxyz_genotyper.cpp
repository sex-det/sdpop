/*
	This file is part of SDpop.	

	SDpop is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SDpop is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "types.h"
#include "reading.h"
#include "calc.h"
#include <cerrno>
#include <string>
#include <vector>
#include <algorithm>

extern int NAME_LEN;


//1) read sdpop output and select genes considered as sex-linked
//2) read genotype file 
std::string StringToUpper(std::string strToConvert)
{
    std::transform(strToConvert.begin(), strToConvert.end(), strToConvert.begin(), ::toupper);

    return strToConvert;
}
std::string StringToLower(std::string strToConvert)
{
    std::transform(strToConvert.begin(), strToConvert.end(), strToConvert.begin(), ::tolower);

    return strToConvert;
}

struct uppercase {
   void operator()(char& c) { c = toupper((unsigned char)c); }
};

struct lowercase {
   void operator()(char& c) { c = tolower((unsigned char)c); }
};

//void outfunction(FILE *outfile, char *contig, int s,int t,int fixed,int ns,double pi_X,double pi_Y,double divergence,char *sequenceX,char *sequenceY){
//	fprintf(outfile,">%s_X %d %d %d %f %f\n",contig,t,fixed,ns,pi_X/ns,divergence/ns);
//	sequenceX[s]='\0';
//	sequenceY[s]='\0';
//	fprintf(outfile,"%s\n",sequenceX);
//	fprintf(outfile,">%s_Y %d %d %d %f %f\n",contig,t,fixed,ns,pi_Y/ns,divergence/ns);
//	fprintf(outfile,"%s\n",sequenceY);
//	
//	free(sequenceX);
//	free(sequenceY);
//}
//
//void outfunction_ref(FILE *outfile, char *contig, int s,int t,int fixed,int ns,double pi_X,double pi_Y,double divergence,int Xref,int Yref,char *sequenceX,char *sequenceY){
//	fprintf(outfile,">%s_X %d %d %d %f %f %d\n",contig,t,fixed,ns,pi_X/ns,divergence/ns,Xref);
//	sequenceX[s]='\0';
//	sequenceY[s]='\0';
//	fprintf(outfile,"%s\n",sequenceX);
//	fprintf(outfile,">%s_Y %d %d %d %f %f %d\n",contig,t,fixed,ns,pi_Y/ns,divergence/ns,Yref);
//	fprintf(outfile,"%s\n",sequenceY);
//	
//	free(sequenceX);
//	free(sequenceY);
//}


void treat_and_output(ContigGenotypesA contiggenotypes, ContigA contigsites, std::vector<std::vector<double>> contigf, const double sitethreshold1, const double sitethreshold2, FILE *outfile){
	int npos=contiggenotypes.genotypes.size();
	int nind=contiggenotypes.individuals.size();
	int amax,i,ii,j,div,t;
	std::string sequenceX;
	std::string sequenceY;
	int nsites=0,nmsites=0,nmlen=0,fixed=0,totreflen=0,refX=0,refY=0;
	double pi_X=0,pi_Y=0,divergence=0;
	
	for (j=0; j<npos; j++){
		if(contiggenotypes.genotypes[j].alleles.size()==0){
				sequenceX+="N";
				sequenceY+="N";
				totreflen+=1;
				continue;
		}
		int reflen=contiggenotypes.genotypes[j].alleles[0].length();
		totreflen+=reflen;
		amax=-1;
		for (i=0; i<nind; i++){ //find maximum allele number
			for (ii=0; ii<2; ii++){
				if(contiggenotypes.genotypes[j].individualgenotypes[i].allele[ii]>amax){
					amax=contiggenotypes.genotypes[j].individualgenotypes[i].allele[ii];
				}
			}
		}
		int allelecount[amax+1] = {0};
		//count the number of alleles
		for (i=0; i<nind; i++){ //loop through all chromosomes
			for (ii=0; ii<2; ii++){
				if(contiggenotypes.genotypes[j].individualgenotypes[i].allele[ii]>=0){
					(allelecount[contiggenotypes.genotypes[j].individualgenotypes[i].allele[ii]])++;
				}
			}
		}
		div=0;
		for(i=0;i<=amax;i++){
			if(allelecount[i]>0){
				div++;
			}
		}
		if(div==0){
			for(i=0;i<reflen;i++){
				sequenceX+="N";
				sequenceY+="N";
			}
		}
		else if(div==1){
			//all individuals have the same allele ; suppose that X and Y are the same
			for(i=0;i<=amax;i++){
				if(allelecount[i]>0){
					sequenceX+=contiggenotypes.genotypes[j].alleles[i];
					sequenceY+=contiggenotypes.genotypes[j].alleles[i];
					int allelelen=contiggenotypes.genotypes[j].alleles[i].length();
					if(allelelen == 1 ){
						nsites++;
					}
					else {
						nmsites++;
						nmlen+=allelelen;
					}
					break; //as there is only one allele
				}
			}
		}
		else {
			Varsite tempvarsite;
			for(t=0; t<contigsites.varsites.size(); t++){
				if(contigsites.varsites[t].position==contiggenotypes.genotypes[j].position){
					tempvarsite=contigsites.varsites[t];
					break;
				}
			}
			if (t < contigsites.varsites.size() ){ // we have found the site
				std::string allele1=tempvarsite.alleles[0];
				std::string allele2=tempvarsite.alleles[1];
				int allele1len=allele1.length();
				int allele2len=allele2.length();
				//adjust lengths if the two alleles are included ; if none of the alleles are to be included, use the reference allele length
				int maxlen=0;
				if ( contigf[t][0]>=sitethreshold2 ){
					maxlen=allele1len;
				}
				if ( contigf[t][1]>=sitethreshold2 ){
					maxlen=allele1len;
				}
				if ( 1-contigf[t][0]>=sitethreshold2 ){
					maxlen=(maxlen > allele2len ? maxlen : allele2len);
				}
				if ( 1-contigf[t][1]>=sitethreshold2 ){
					maxlen=(maxlen > allele2len ? maxlen : allele2len);
				}
				for(i=allele1len;i<maxlen;i++){
					allele1+="-";
				}
				for(i=allele2len;i<maxlen;i++){
					allele2+="-";
				}
				if(maxlen == 0){
					maxlen=reflen;
				}
				if(contigf[t][0]>=sitethreshold1){
					sequenceX+=StringToUpper(allele1);
				}
				else if(contigf[t][0]>=sitethreshold2){
					sequenceX+=StringToLower(allele1);
				}
				else if(contigf[t][0]<=1-sitethreshold1){
					sequenceX+=StringToUpper(allele2);
				}
				else if(contigf[t][0]<=1-sitethreshold2){
					sequenceX+=StringToLower(allele2);
				}
				else {
					for(i=0;i<maxlen;i++){
						sequenceX+="n";
					}
				}
				if(contigf[t][1]>=sitethreshold1){
					sequenceY+=StringToUpper(allele1);
				}
				else if(contigf[t][1]>=sitethreshold2){
					sequenceY+=StringToLower(allele1);
				}
				else if(contigf[t][1]<=1-sitethreshold1){
					sequenceY+=StringToUpper(allele2);
				}
				else if(contigf[t][1]<=1-sitethreshold2){
					sequenceY+=StringToLower(allele2);
				}
				else {
					for(i=0;i<maxlen;i++){
						sequenceY+="n";
					}
				}
				if(allele1len==1 && allele2len==1){ //only for SNPs
					nsites++;
					pi_X+=2*contigf[t][0]*(1.-contigf[t][0]);
					pi_Y+=2*contigf[t][1]*(1.-contigf[t][1]);
					divergence+=contigf[t][0]*(1.-contigf[t][1])+contigf[t][1]*(1.-contigf[t][0]);
					if( contigf[t][0]>=sitethreshold1 && contigf[t][1]<=1-sitethreshold1 ) {
						fixed++;
						if(contiggenotypes.genotypes[j].alleles[0]==StringToUpper(allele1)){
							refX++;
						}
						else if(contiggenotypes.genotypes[j].alleles[0]==StringToUpper(allele2)){
							refY++;
						}
					}
					else if(contigf[t][1]>=sitethreshold1 && contigf[t][0]<=1-sitethreshold1 ) {
						fixed++;
						if(contiggenotypes.genotypes[j].alleles[0]==StringToUpper(allele1)){
							refY++;
						}
						else if(contiggenotypes.genotypes[j].alleles[0]==StringToUpper(allele2)){
							refX++;
						}
					}
				}
				else {
					nmsites++;
					nmlen+=maxlen;
				}

				////				if( (contigf[t][0]>=sitethreshold1 || contigf[t][0]<=1-sitethreshold1) && (contigf[t][1]>=sitethreshold1 || contigf[t][1]<=1-sitethreshold1) ){
				//if( contigf[t][0]>=sitethreshold1 && contigf[t][1]<=1-sitethreshold1) { // nuc0 is fixed on X, nuc1 on Y
				//	(*fixed)++;
				//	if(toupper(nuc[k][t][0])==toupper(alleles[0])){
				//		(*Xref)++;
				//	}
				//	else if(toupper(nuc[k][t][1])==toupper(alleles[0])){
				//		(*Yref)++;
				//		printf("Yref: %d\n",s);
				//	}
				//}
				//else if (contigf[t][1]>=sitethreshold1 && contigf[t][0]<=1-sitethreshold1) { //nuc1 is fixed on X, nuc0 on Y
				//	(*fixed)++;
				//	if(toupper(nuc[k][t][1])==toupper(alleles[0])){
				//		(*Xref)++;
				//	}
				//	else if(toupper(nuc[k][t][0])==toupper(alleles[0])){
				//		(*Yref)++;
				//		printf("Yref: %d\n",s);
				//	}
				//}
				
			}
			else { //there is some polymorphism that has passed unnoticed 
				sequenceX+="X";
				sequenceY+="X";
			}
			
		}
		j=j+reflen-1;
	}
	//fprintf(outfile,">%s_X %d %d %d %f %f %d\n",contig,t,fixed,ns,pi_X/ns,divergence/ns,Xref);
	
	fprintf(outfile,">%s_X %d %d %d %d %d %d %d %d %f %f\n",contigsites.name.data(),contigsites.varsites.size(),sequenceX.length(),totreflen,nsites,nmsites,nmlen,fixed,refX,pi_X/nsites,divergence/nsites);
	fprintf(outfile,"%s\n",sequenceX.data());
	fprintf(outfile,">%s_Y %d %d %d %d %d %d %d %d %f %f\n",contigsites.name.data(),contigsites.varsites.size(),sequenceY.length(),totreflen,nsites,nmsites,nmlen,fixed,refY,pi_Y/nsites,divergence/nsites);
	fprintf(outfile,"%s\n",sequenceY.data());
	//outfunction_ref(outfile,&contig[k*NAME_LEN],s,t,fixed,ns,pi_X,pi_Y,divergence,Xref,Yref,sequenceX,sequenceY);
}

int main(int argc, char *argv[]) 
{
	FILE *sdpfile,*genfile,*outfile;
	int CUR_MAX = 4095;
	char *line = (char *)calloc((size_t)CUR_MAX,sizeof(char));
	char *tmpline = (char *)calloc((size_t)CUR_MAX,sizeof(char));
	char ch='a',gencontig[NAME_LEN],tcont[NAME_LEN],word[NAME_LEN],chrom[NAME_LEN],oldchrom[NAME_LEN],**name;
	int npolysites,toread=0,ncontigs,ngencontigs;
	double pp,contigthreshold,sitethreshold1,sitethreshold2;
	char nuc1[NAME_LEN],nuc2[NAME_LEN];
	int i,j,k,l,t,found,pos,ni,chrpos;
	int count,length,ff[2],xy,nwords,posterior_field,mean;
	int *sex;
	int datafmt=1;
	int READS2SNP=0;
	int VCF=1;
	int FASTA=2;
	std::vector<ContigA> contigs;
	std::vector<std::vector<std::vector<double>>> f;
	std::vector<double> contigpp;
	std::vector<std::string> contig;
	std::vector<std::vector<std::vector<std::string>>> nuc;

	
	for(i=0;i<argc;i++) {
		fprintf(stdout,"%s ",argv[i]);
	}
	fprintf(stdout,"\n");
	
	if (argc != 11) {
		fprintf(stdout,"Usage: %s sdpopfile reffile reftype outfile contigthreshold sitethreshold1 sitethreshold2 posterior_field system max/mean\n",argv[0]);
		exit(1);
	}
	
	if((sdpfile=fopen(argv[1],"r"))==NULL){
		fprintf(stderr,"error opening input file %s\n",argv[1]);
		exit(1);
	}
	if((genfile=fopen(argv[2],"r"))==NULL){
		fprintf(stderr,"error opening input file %s\n",argv[2]);
		exit(1);
	}
	if(strcmp(argv[3],"r")==0 || strcmp(argv[3],"1")==0) {
		datafmt=READS2SNP;
	}
	else if (strcmp(argv[3],"v")==0 || strcmp(argv[3],"3")==0) {
		datafmt=VCF;
	}
	else if (strcmp(argv[3],"f")==0 || strcmp(argv[3],"2")==0) {
		datafmt=FASTA;
	}
	else {
		fprintf(stdout,"Usage: %s infile outfile type namelist1 namelist2 (r)\n",argv[0]);
		fprintf(stderr,"type should be either \"r\" or \"1\" for reads2snp genotype file, or \"f\" or \"2\" for fasta file, or \"v\" or \"3\" for VCF.\n");
		exit(1);	
	}

	if((outfile=fopen(argv[4],"w"))==NULL){
		fprintf(stderr,"error opening output file %s\n",argv[4]);
		exit(1);
	}	
	contigthreshold=atof(argv[5]);
	if(contigthreshold<0 || contigthreshold>1){
		fprintf(stderr,"Error: contigthreshold should be between 0 and 1; value given: %s\n",argv[5]);
		exit(1);
	}
	sitethreshold1=atof(argv[6]);
	if(sitethreshold1<0.5 || sitethreshold1>1){
		fprintf(stderr,"Error: sitethreshold1 should be between 0.5 and 1; value given: %s\n",argv[6]);
		exit(1);
	}
	sitethreshold2=atof(argv[7]);
	if(sitethreshold2<0.5 || sitethreshold2>sitethreshold1){
		fprintf(stderr,"Error: sitethreshold2 should be between 0.5 and sitethreshold1; value given: %s\n",argv[7]);
		exit(1);
	}
	if(strcmp(argv[9],"XY")==0){
		xy=1;
	}
	else if(strcmp(argv[9],"ZW")==0){
		xy=0;
	}
	else {
		fprintf(stderr,"Error: system should be either \"XY\" or \"ZW\" ; value given: \"%s\"\n",argv[9]);
		exit(1);
	}
	if(strcmp(argv[10],"mean")==0){
		mean=1;
	}
	else if(strcmp(argv[10],"max")==0){
		mean=0;
	}
	else {
		fprintf(stderr,"Error: choose either \"max\" or \"mean\" ; value given: \"%s\"\n",argv[10]);
		exit(1);
	}
		
		
	k=-1;
	
	fprintf(stdout,"Reading data from file %s...\n",argv[1]);
	//sdpop file reading
	l=0; //line number
	ch='a';
	std::vector<std::vector<double>> tempfcontig;
	while (ch != EOF) { //loop through the file
		ch='a';
		count = 0;
		length = 0;
		while ( (ch != '\n') && (ch != EOF) ) { //loop through the line
			if(count == CUR_MAX) { // time to expand (for unexepectedly large line lengths) ?
				CUR_MAX *= 2; 
				count = 0;
				line = (char *)realloc(line, sizeof(char) * CUR_MAX); 
				tmpline = (char *)realloc(tmpline, sizeof(char) * CUR_MAX); 
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
					if(strcmp(word,argv[8])==0){
						posterior_field=nwords;
						break;
					}
					nwords++;
					strcpy(line,tmpline);
				}
//				printf("%d\n",posterior_field);
				if(posterior_field<0){
					printf("Error: word \"%s\" not found in line %d\n",argv[8],l);
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
					if(xy){
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

				fprintf(outfile,"#contig_name N_poly_sites length reference_length snp_genotyped_length non_snps non_snp_length N_fixed_diff N_refalleles pi divergence\n");
			}
		}
		
		if ( line[0] == '>') {
			if(toread==1){
				f.push_back(tempfcontig);
				tempfcontig.clear();
			}
			ContigA tempcontig;
			k++;
			//Start preparing to read a new contig
			if ( strlen(line) >= NAME_LEN ) {
				fprintf(stderr,"Error: a line name is too long (more than %d characters), line %d.\n",NAME_LEN-1,l);
				exit(1);
			}
//			printf("%s",line);
			sscanf(line,">%s\t%d\t%[^\n]",tcont,&npolysites,tmpline);
			nwords=2;
			strcpy(line,tmpline);
			while ( sscanf(line,"%[^\t ]%*[\t ]%[^\n]",word,tmpline)==2)	{
				strcpy(line,tmpline);
				if(nwords==posterior_field){
					pp=atof(word);
					break;
				}
				nwords++;
			}	
			if(pp>=contigthreshold){
//				printf("%f\n",contigpp[k]);
				toread=1;
				t=0;
				contigpp.push_back(pp);
				tempcontig.name=tcont;
				contigs.push_back(tempcontig);
			}
			else {
				k--;
				toread=0;
			}

		}
		else if(isdigit(line[0]) && toread==1){
			ContigA & current_contig = contigs.back();
			if(t>=npolysites){
				fprintf(stderr,"Error: contig %s does seem to have more polymorphic sites than announced\n",current_contig.name.data());
				exit(1);				
			}
//			printf("%s",line);
			Varsite tempvarsite;
			std::vector<double> tempfsite;
			if(sscanf(line,"%d\t%[^,],%s\t%d\t%d\t%d\t%d\t%d\t%d\t%[^\n]",&tempvarsite.position,nuc1,nuc2,
						&tempvarsite.genotypes_by_sex[N11F],&tempvarsite.genotypes_by_sex[N12F],
					&tempvarsite.genotypes_by_sex[N22F],&tempvarsite.genotypes_by_sex[N11M],&tempvarsite.genotypes_by_sex[N12M],
					&tempvarsite.genotypes_by_sex[N22M],tmpline)!=10){
					fprintf(stderr,"Error reading line %d\n",l);
					fprintf(stderr,"Line: %s\n",line);
					exit(1);		
			}
			tempvarsite.alleles.push_back(nuc1);				
			tempvarsite.alleles.push_back(nuc2);				
			current_contig.varsites.push_back(tempvarsite);
			strcpy(line,tmpline);
			nwords=8;
			j=0;
			while ( sscanf(line,"%[^\t ]%*[\t ]%[^\n]",word,tmpline)==2)	{
				strcpy(line,tmpline);
				if(nwords==ff[j]){
					tempfsite.push_back(atof(word));
					j++;
					if(j>1){
						break;
					}
				}
				nwords++;
			}
			tempfcontig.push_back(tempfsite);
			t++;
		}
	}
	if(toread==1){
		f.push_back(tempfcontig);
		tempfcontig.clear();
	}
	ncontigs=k+1;

	fprintf(stdout,"%d contigs selected.\n",ncontigs);
	fclose(sdpfile);

	if(ncontigs>0){
		fprintf(stdout,"Reading data from file %s...\n",argv[2]);
		if(datafmt==READS2SNP){
			ContigGenotypesA contiggenotypes;
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
						line = (char *)realloc(line, sizeof(char) * CUR_MAX); 
						tmpline = (char *)realloc(tmpline, sizeof(char) * CUR_MAX); 
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
						treat_and_output(contiggenotypes,contigs[found],f[found],sitethreshold1,sitethreshold2,outfile);
						free(name);
						free(sex);
					}
					sscanf(line,">%s",gencontig);
					found=-1;
					for(k=0;k<ncontigs;k++){
						if(strcmp(gencontig,contigs[k].name.data())==0){
							found=k;
							ngencontigs+=1;
							contiggenotypes.individuals.clear();								
							contiggenotypes.genotypes.clear();								
							break;
						}
					}
				}
				else if ( line[0] == 'p' && found>=0){
					name=getgennames(line, &ni); //returns all names present in line, and the number of names ni
					if((sex=(int *)calloc((size_t)(ni),sizeof(int)))==NULL) { 
						fprintf(stderr,"error in memory allocation\n");
						exit(1);
					}
					for(i=0;i<ni;i++){
						sex[i]=1; //use all individuals
						Individual tmpindividual;
						tmpindividual.name=name[i];
						tmpindividual.sex=sex[i];
						contiggenotypes.individuals.push_back(tmpindividual);
					}
					
				}

				else if ( line[0] != 'p' && found>=0 ){
					sscanf(line,"%d",&pos);
					contiggenotypes.genotypes.push_back(gengenotypes(pos,ni,sex,line));
				}
			}
			if(found>=0){ // print last sequences
				treat_and_output(contiggenotypes,contigs[found],f[found],sitethreshold1,sitethreshold2,outfile);
				free(name);
				free(sex);
			}
		}
/*		else if (datafmt==FASTA){
			ch='a';
			l=0;
			found=-1;
			k=0;
			ngencontigs=0;
			while (ch != EOF) { //loop through the fasta file
				ch='a';
				count = 0;
				length = 0;
				while ( (ch != '\n') && (ch != EOF) ) { //loop through the line
					if(count == CUR_MAX) { // time to expand (for unexepectedly large line lengths) ?
						CUR_MAX *= 2; 
						count = 0;
						line = (char *)realloc(line, sizeof(char) * CUR_MAX); 
						tmpline = (char *)realloc(tmpline, sizeof(char) * CUR_MAX); 
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
				if (strncmp(line,">",1)==0) { //line with fasta name
//							printf("%s\n",line);
					sscanf(line,">%s",chrom);
					//fprintf(stdout,"%s\n",chrom);						
					if(found>=0){ // do the actual work and print last sequences
						ns=strlen(sequenceX);
						for(t=0;t<npolysites[k];t++){
							s=polysite[k][t][0];
							if(s>=strlen(sequenceX)){
								fprintf(stderr,"SNP position (%d) is outside the reference for sequence %s (which has length %d in fasta file)\n",s,&contig[k*NAME_LEN],(int)strlen(sequenceX));
								exit(1);				
							}
							refallele=toupper(sequenceX[s]);
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
							if( f[k][t][0]>=sitethreshold1 && f[k][t][1]<=1-sitethreshold1) { // nuc0 is fixed on X, nuc1 on Y
								fixed++;
								if(toupper(nuc[k][t][0])==refallele){
									Xref++;
								}
								else if(toupper(nuc[k][t][1])==refallele){
									Yref++;
									printf("Yref: %d\n",s);
								}
							}
							else if (f[k][t][1]>=sitethreshold1 && f[k][t][0]<=1-sitethreshold1) { //nuc1 is fixed on X, nuc0 on Y
								fixed++;
								if(toupper(nuc[k][t][1])==refallele){
									Xref++;
								}
								else if(toupper(nuc[k][t][0])==refallele){
									Yref++;
									printf("Yref: %d\n",s);
								}
							}
							for(s=0;s<strlen(sequenceX);s++){
								if(toupper(sequenceX[s])=='N'){
									ns--;
								}
							}
						}
						outfunction_ref(outfile,&contig[k*NAME_LEN],s,t,fixed,ns,pi_X,pi_Y,divergence,Xref,Yref,sequenceX,sequenceY);
					}
					
					found=-1;
					for(k=0;k<ncontigs;k++){
//							printf("%s\n",chrom);
						if(strcmp(chrom,&contig[k*NAME_LEN])==0){
							printf("%s\n",chrom);
							found=k;
							t=0;
							s=0;
							ns=0;
							fixed=0;
							Xref=0;
							Yref=0;
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
				else { // line contains nucleotides
					if(found>-1){
						if(strlen(sequenceX)+strlen(line)>=sites_allocated){//time to expand
							sites_allocated+=strlen(line);
							if((sequenceX=(char *)realloc(sequenceX,sizeof(char)*(sites_allocated+1)))==NULL){
								fprintf(stderr,"error in memory allocation (realloc)\n");
								exit(1);
							}
							if((sequenceY=(char *)realloc(sequenceY,sizeof(char)*(sites_allocated+1)))==NULL){
								fprintf(stderr,"error in memory allocation (realloc)\n");
								exit(1);
							}				
						}
						sprintf(sequenceX,"%s",line);
						sprintf(sequenceY,"%s",line);
					}
				}
			}
		}*/		
		else if (datafmt==VCF){
			ContigGenotypesA contiggenotypes;
			ch='a';
			l=0;
			found=-1;
			k=0;
			ngencontigs=0;
			strcpy(oldchrom," ");
			while (ch != EOF) { //loop through the vcf file
				
				ch='a';
				count = 0;
				length = 0;
				while ( (ch != '\n') && (ch != EOF) ) { //loop through the line
					if(count == CUR_MAX) { // time to expand (for unexepectedly large line lengths) ?
						CUR_MAX *= 2; 
						count = 0;
						line = (char *)realloc(line, sizeof(char) * CUR_MAX); 
						tmpline = (char *)realloc(tmpline, sizeof(char) * CUR_MAX); 
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
				if( strncmp(line,"#CHROM",6)==0 ) { //only one line with individual names in vcf file				
					name=getvcfnames(line,&ni);
					if((sex=(int *)calloc((size_t)(ni),sizeof(int)))==NULL) { 
						fprintf(stderr,"error in memory allocation\n");
						exit(1);
					}
					for(i=0;i<ni;i++){
						sex[i]=1; //use all individuals
						Individual tmpindividual;
						tmpindividual.name=name[i];
						tmpindividual.sex=sex[i];
						contiggenotypes.individuals.push_back(tmpindividual);
					}
								
				}
				else if (strncmp(line,"#",1)!=0) { //line contains genotype data
					sscanf(line,"%s\t%d",chrom,&chrpos);
//					fprintf(stdout,"%s\n",chrom);
					if(strcmp(chrom,oldchrom)==0){ //still the same contig as before
						if(found>-1){
//							fprintf(stdout,"%d %d %d\n",s,t,chrpos);
							GenotypesA tempgenotypes;
							tempgenotypes.position=chrpos;
							tempgenotypes.alleles=vcfsnpindel(line);
							tempgenotypes.individualgenotypes=vcfalleles(ni,sex,line);
							contiggenotypes.genotypes.push_back(tempgenotypes);
						}
					}
					else {
						strcpy(oldchrom,chrom);
						if(found>=0){ // print last sequences
							treat_and_output(contiggenotypes,contigs[found],f[found],sitethreshold1,sitethreshold2,outfile);
						}
						
						found=-1;
						for(k=0;k<ncontigs;k++){
							if(strcmp(chrom,contigs[k].name.data())==0){
								printf("%s\n",chrom);
								found=k;
								ngencontigs+=1;
								contiggenotypes.genotypes.clear();								
								break;
							}
						}
						if(found>=0){
							GenotypesA tempgenotypes;
							tempgenotypes.position=chrpos;
							tempgenotypes.alleles=vcfsnpindel(line);
							tempgenotypes.individualgenotypes=vcfalleles(ni,sex,line);
							contiggenotypes.genotypes.push_back(tempgenotypes);
						}
					}
				}
			}
			if(found>=0){ // print last sequences
				treat_and_output(contiggenotypes,contigs[found],f[found],sitethreshold1,sitethreshold2,outfile);
			}			
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