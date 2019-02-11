#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "reading.h"
#include "types.h"

#define N11 1
#define N12 2
#define N22 3
#define READS2SNP 1
#define GBS_FILTERED 2
#define VCF 3

extern int NAME_LEN;
			
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

Contig polyfilter(ContigGenotypes contiggenotypes, int *n3, int *n4, double *theta) {//number of positions and individuals,
	// genotype matrix, pointer to the number of polymorphic sites
	// function treats one contig
	char nucvec[4] = {0};
	int i,j,ii,randbit;
	int div,n11f,n12f,n22f,n11m,n12m,n22m,nA,nT,nG,nC,nN;
	char pos1,pos2,nuc1,nuc2;
	int n,n2,nt,nf,npos,nind;
	double a;
	Contig contig;
	
	contig.name=contiggenotypes.name;
	
	nt=0;
	nf=0;
	*theta=0.;
	
	npos=contiggenotypes.genotypes.size();
	nind=contiggenotypes.individuals.size();
	
	n=0;
	for (j=0; j<npos; j++){
		//count the number of alleles : A, T, G, C, N
		nA=nT=nG=nC=nN=0;
		for (i=0; i<nind; i++){ //loop through all chromosomes
			for (ii=0; ii<2; ii++){
				switch (contiggenotypes.genotypes[j].individualgenotypes[i].nucleotides[ii]) {
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
		if (div == 1) { //no polymorphism
			nuc1=nucvec[0];
			for (i=0; i<nind; i++){ //loop through all chromosomes
				pos1=contiggenotypes.genotypes[j].individualgenotypes[i].nucleotides[0];
				pos2=contiggenotypes.genotypes[j].individualgenotypes[i].nucleotides[1];
				if (contiggenotypes.individuals[i].sex==FEMALE) {
					if (pos1 == pos2 && pos1 == nucvec[0]) { //homozygous
						n11f++;
					}
				}
				else if (contiggenotypes.individuals[i].sex==MALE) {
					if (pos1 == pos2 && pos1 == nucvec[0]) { //homozygous
						n11m++;
					}
				}
			}
		}
		else if (div == 2) { //simple polymorphism
			//randomize 
			nuc1=nucvec[0];
			nuc2=nucvec[1];
			randbit=rand()%2;
			nucvec[randbit]=nuc1;
			nucvec[!randbit]=nuc2;
			for (i=0; i<nind; i++){ //loop through all chromosomes
				pos1=contiggenotypes.genotypes[j].individualgenotypes[i].nucleotides[0];
				pos2=contiggenotypes.genotypes[j].individualgenotypes[i].nucleotides[1];
				if (contiggenotypes.individuals[i].sex==FEMALE) {
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
				else if (contiggenotypes.individuals[i].sex==MALE) {
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
		else if (div ==3) {
			nt++;
		}
		else if (div ==4) {
			nf++;
		}
//		if (div == 1 || div == 2 ) {
		if (div == 2 ) {
			SNP tempsnp;
			tempsnp.position=contiggenotypes.genotypes[j].position; //position
			tempsnp.genotypes_by_sex[N11F]=n11f; //female counts
			tempsnp.genotypes_by_sex[N12F]=n12f; 
			tempsnp.genotypes_by_sex[N22F]=n22f; 			
			tempsnp.genotypes_by_sex[N11M]=n11m; //male counts
			tempsnp.genotypes_by_sex[N12M]=n12m; 
			tempsnp.genotypes_by_sex[N22M]=n22m; 	
			tempsnp.alleles[0]=DNA2int(nucvec[0]);
			if (div == 2 ) {
				tempsnp.alleles[1]=DNA2int(nucvec[1]);
			}
			else if (div == 1 ) {
				tempsnp.alleles[1]=DNA2int(nucvec[0]);
			}
			contig.snps.push_back(tempsnp);
			n++;
		}
	}
	*n3=nt;
	*n4=nf;
	*theta/=(double)npos;
	return contig;
}

void write_contig(FILE *outfile, Contig contig, int n3, int n4, double theta)
{
	int t,i;
	int npolysites;
	npolysites=contig.snps.size();
	fprintf(outfile,">%s\t%d\t%d\t%d\t%f\n",contig.name.data(),npolysites,n3,n4,theta);
	if(npolysites>0) {
		for (t=0; t<npolysites; t++){
			fprintf(outfile,"%d\t",contig.snps[t].position);
			fprintf(outfile,"%c%c\t",int2DNA(contig.snps[t].alleles[0]),int2DNA(contig.snps[t].alleles[1]));
			for (i=0; i<6; i++) {
				fprintf(outfile,"%d\t",contig.snps[t].genotypes_by_sex[i]);
			}
			fprintf(outfile,"\n");
		}
	}
}

int main(int argc, char *argv[]) 
{
	FILE *fp,*outfile;
	std::string line;
	int n3,n4;
	double theta;
	int ncontigs,totsites=0;
	int i,j,k,l,ni,firstcontig,pos;
	int nfem,nmal,ifem,imal,nplus=0,chrpos;
	int *sex,*foundsex;
	int c;
	char **name,**femname,**malname;
	int *ffound,*mfound;
	int nJ,nf,nm,nfound; //number of individuals
	int randomise=0,ri,tempsex,datafmt;
	int val,nnuc,maxnuc=5; //A, T, G, C, or N
	char nuc[maxnuc];
	char chrom[NAME_LEN];

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
	else if (strcmp(argv[3],"v")==0 || strcmp(argv[3],"3")==0) {
		datafmt=VCF;
	}
	else {
		fprintf(stdout,"Usage: %s infile outfile type namelist1 namelist2 (r)\n",argv[0]);
		fprintf(stderr,"type should be either \"r\" or \"1\" for reads2snp genotype file, or \"g\" or \"2\" for filtered GBS data, or \"v\" or \"3\" for VCF.\n");
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
		exit(6);
	}
	if((mfound=(int *)calloc((size_t)(nmal),sizeof(int)))==NULL) {
		fprintf(stderr,"error in memory allocation\n");
		exit(7);
	}
	for (ifem=0;ifem<nfem;ifem++) {
		ffound[ifem]=0;
	}
	for (imal=0;imal<nmal;imal++) {
		mfound[imal]=0;
	}

	//genotype file reading
	l=0; //line number
		
	firstcontig=1;
	k=0;
	
	fprintf(stdout,"Reading data...\n");
	
	if(datafmt==READS2SNP) {
		ContigGenotypes contiggenotypes;
		
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
			
			if ( strncmp(line.data(),">",1)==0){
				if (firstcontig==0) {	//Filtering polymorphisms for the last read contig
					nJ=pos;
					Contig contig=polyfilter(contiggenotypes,&n3,&n4,&theta);
					contiggenotypes.individuals.clear();
					contiggenotypes.genotypes.clear();
					write_contig(outfile,contig,n3,n4,theta);
					totsites+=contig.snps.size();
					free(sex);
					free(foundsex);
					for (i=0; i < ni; i++){
						free(name[i]);
					}
					free(name);
					k++;
					if(k % 5000 == 0){
						fprintf(stdout,"%d contigs, %d polymorphic sites, and still reading...\n",k,totsites);
					}
				}
				firstcontig=0;
				
				//Start preparing to read a new contig
				sscanf(line.data(),">%s",chrom);
				contiggenotypes.name=chrom;
			}
			else if ( strncmp(line.data(),"p",1)==0){
				name=getgennames(line.data(), &ni); //returns all names present in line, and the number of names ni
				if(ni<=0){
					fprintf(stdout, "No individual names found in line %d; format error ?\nexiting\n",l);
					exit(1);
				}
				
				//Find out the sex of the individuals.
				if((sex=(int *)calloc((size_t)(ni),sizeof(int)))==NULL) { 
					fprintf(stderr,"error in memory allocation\n");
					exit(29);
				}
				if((foundsex=(int *)calloc((size_t)(ni),sizeof(int)))==NULL) { 
					fprintf(stderr,"error in memory allocation\n");
					exit(30);
				}
				
				nfound=findsex(name, ni, femname, nfem, malname, nmal, sex, foundsex, ffound, &nf, mfound, &nm);
				
				if(randomise){ //Fisher-Yates shuffling
					for (i=0;i<ni-1;i++) {
						if(sex[i]>=0){
							while(sex[ri=i+rand()/(RAND_MAX/(ni-i)+1)]<0){
							}
							tempsex=sex[i];
							sex[i]=sex[ri];
							sex[ri]=tempsex;
						}
					}
				}
				for(i=0;i<ni;i++){
					if ( sex[i]>=0 ) {
						Individual tmpindividual;
						tmpindividual.name=name[i];
						tmpindividual.sex=sex[i];
						contiggenotypes.individuals.push_back(tmpindividual);
					}
				}
				
				
				if(ni>nfem+nmal){
					if(nplus==0) {
						fprintf(stdout,"Contig %s seems to have more observations than names given:\n",contiggenotypes.name.data());
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
				j=0;
			}
			else { //line contains genotype data
				sscanf(line.data(),"%d",&pos);
//				if(pos != j+1){
//					fprintf(stderr,"Error in input file line number %d: expecting to read \"%d ...\"\n",l,j+1);
//					fprintf(stderr,"read \"%s\" instead\n",line.data());
//					exit(1);
//				}
				Genotypes tempgenotypes;
				tempgenotypes.position=pos;
				tempgenotypes.individualgenotypes=gengenotypes(ni,sex,line.data());
				contiggenotypes.genotypes.push_back(tempgenotypes);
				
				j++;
				
			}
		}
		nJ=j;
		Contig contig=polyfilter(contiggenotypes,&n3,&n4,&theta);
		write_contig(outfile,contig,n3,n4,theta);
		totsites+=contig.snps.size();
		free(sex);
		free(foundsex);
		for (i=0; i < ni; i++){
			free(name[i]);
		}
		free(name);
		k++;
		fprintf(stdout,"Read %d contigs and %d polymorphic sites\n",k,totsites);
		
		
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
	}
	else if (datafmt==VCF){
		
		ContigGenotypes contiggenotypes;
		k=0;
		l=0; //line number
		
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
				name=getvcfnames(line.data(),&ni);
				//Find out the sex of the individuals.
				if((sex=(int *)calloc((size_t)(ni),sizeof(int)))==NULL) { 
					fprintf(stderr,"error in memory allocation\n");
					exit(38);
				}
				if((foundsex=(int *)calloc((size_t)(ni),sizeof(int)))==NULL) { 
					fprintf(stderr,"error in memory allocation\n");
					exit(39);
				}
				nfound=findsex(name, ni, femname, nfem, malname, nmal, sex, foundsex, ffound, &nf, mfound, &nm);
				if(nfound==0){
					fprintf(stderr,"Nothing to do; exiting\n");
					exit(0);
				}
					if(randomise){ //Fisher-Yates shuffling
						for (i=0;i<ni-1;i++) {
							if(sex[i]>=0){
							while(sex[ri=i+rand()/(RAND_MAX/(ni-i)+1)]<0){
							}
							tempsex=sex[i];
							sex[i]=sex[ri];
							sex[ri]=tempsex;
							}
						}
					}
				for(i=0;i<ni;i++){
					if ( sex[i]>=0 ) {
						Individual tmpindividual;
						tmpindividual.name=name[i];
						tmpindividual.sex=sex[i];
						contiggenotypes.individuals.push_back(tmpindividual);
					}
				}
				if(ni>nfem+nmal){
					if(nplus==0) {
						fprintf(stdout,"There seem to be more observations than names given\n");
						fprintf(stdout,"found %d individuals in file, and %d names in command line\n",ni,nfem+nmal);
						fprintf(stdout,"No sex found for individual(s): ");
						for (i=0; i<ni; i++){
							if ( sex[i]==-1 ) {
								fprintf(stdout,"%s ",name[i]);
							}		
						}
						fprintf(stdout,"\n");
						nplus=1;
					}
				}
				if(nfound==0){
					fprintf(stderr,"Individual names given are not found in the file. Nothing to do; exiting\n");
					exit(0);
				}
				
//				strcpy(contname,"\0");
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
					if(randomise){ //Fisher-Yates shuffling
						for (i=0;i<ni-1;i++) {
							if(sex[i]>=0){
							while(sex[ri=i+rand()/(RAND_MAX/(ni-i)+1)]<0){
							}
							tempsex=sex[i];
							sex[i]=sex[ri];
							sex[ri]=tempsex;
							}
						}
					}
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
						Contig contig=polyfilter(contiggenotypes,&n3,&n4,&theta);
							contiggenotypes.genotypes.clear();
							write_contig(outfile,contig,n3,n4,theta);
							totsites+=contig.snps.size();
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
			Contig contig=polyfilter(contiggenotypes,&n3,&n4,&theta);
			write_contig(outfile,contig,n3,n4,theta);
			totsites+=contig.snps.size();
		}
		ncontigs=k+1;
		fprintf(stdout,"Read %d contigs and %d polymorphic sites\n",ncontigs,totsites);	
		for(i=0; i<ni; i++){
			free(name[i]);
		}
		free(name);
	}

	//All reading has been done : clear memory.
	fclose(fp);

	
	//outputting
//	for (k=0;k<ncontigs;k++) {
//		write_contig(outfile,datafmt,&contig[k*NAME_LEN],npolysites[k],n3[k],n4[k],theta[k],sitename[k],polysite[k]);
//	}
	
	free(ffound);
	free(mfound);
	fclose(outfile);
	
	return 0;

}
