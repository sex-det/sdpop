#include "types.h"
#include "reading.h"

int NAME_LEN = 9000;

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

//int read_cnt2(FILE *fp, int namelen, int chromosomes, char **contig_p, int **npolysites_p, int ****polysite_p) 
int read_cnt2(FILE *fp, int namelen, int chromosomes, std::vector<Contig>& contigs) 
//for reading cnt files that have the identity (A,T,C or G) of the alleles
{
	char tcont[namelen];
	std::string line;
	char nuc1,nuc2;
	int c,l;

	//reading counts from file 
	l=0; //line number
	c=0;
		
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
		
		
		if ( line[0] == '>' ){
			Contig tempcontig;

			if ( line.size() >= namelen ) {
				fprintf(stderr,"Error: a contig name is too long (more than %d characters) on line %d.\n",namelen-1,l);
				exit(1);
			}
			sscanf(line.data(),">%s",tcont);
			tempcontig.name=tcont;
			contigs.push_back(tempcontig);
			printf("%s\n",tempcontig.name.data());
//			sscanf(line.data(),">%s",&contig[k*namelen]);
			//			printf("%s\n",&contig[k*namelen]);
		}
		else { //line contains counts
				
				SNP tempsnp;
				if(chromosomes==XY){
					if(sscanf(line.data(),"%d\t%c%c\t%d\t%d\t%d\t%d\t%d\t%d\t",&tempsnp.position,&nuc1,&nuc2,
						&tempsnp.genotypes_by_sex[N11F],&tempsnp.genotypes_by_sex[N12F],
					&tempsnp.genotypes_by_sex[N22F],&tempsnp.genotypes_by_sex[N11M],&tempsnp.genotypes_by_sex[N12M],
					&tempsnp.genotypes_by_sex[N22M])!=9){
						fprintf(stderr,"In readcnt2: Error reading line %d\n",l);
						fprintf(stderr,"Line: %s\n",line.data());
						exit(1);
					}
				}
				else {
					if(sscanf(line.data(),"%d\t%c%c\t%d\t%d\t%d\t%d\t%d\t%d\t",&tempsnp.position,&nuc1,&nuc2,
						&tempsnp.genotypes_by_sex[N11M],&tempsnp.genotypes_by_sex[N12M],
					&tempsnp.genotypes_by_sex[N22M],&tempsnp.genotypes_by_sex[N11F],&tempsnp.genotypes_by_sex[N12F],
					&tempsnp.genotypes_by_sex[N22F])!=9){
						fprintf(stderr,"In readcnt2: Error reading line %d\n",l);
						fprintf(stderr,"Line: %s\n",line.data());
						exit(1);
					}
				}
				tempsnp.alleles[0]=DNA2int(nuc1);				
				tempsnp.alleles[1]=DNA2int(nuc2);				
				Contig & current_contig = contigs.back();
				current_contig.snps.push_back(tempsnp);
		}
	}
	
	return contigs.size();
}

int read_cnt(FILE *fp, 	int namelen, int chromosomes, char **contig_p, int **npolysites_p, int ****polysite_p) 
{
	char *contig;
	int *npolysites,***polysite;
	int CUR_MAX = 4095;
	int NCONTIG_BATCH = 100;
	int ncontigs_allocated,ncontigs;
	char *line = (char*)calloc((size_t)CUR_MAX,sizeof(char));
	char *tmpline = (char*)calloc((size_t)CUR_MAX,sizeof(char));
	int count = 0; 
	int length = 0;
	char ch;
	int i,k,l,firstcontig,t;
	int ni,nI=0;

	//reading counts from file 
	l=0; //line number
	ch='a';
	if((npolysites=(int *)malloc(sizeof(int)*NCONTIG_BATCH))==NULL) {
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}
	if((polysite=(int ***)malloc(sizeof(int **)*NCONTIG_BATCH))==NULL) {
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}
	if((contig=(char *)malloc(sizeof(char)*NCONTIG_BATCH*namelen))==NULL) { 
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}
	ncontigs_allocated=NCONTIG_BATCH;
	
	firstcontig=1;
	k=0;
		
	while (ch != EOF) { //loop through the file
		ch='a';
		count = 0;
		length = 0;
		while ( (ch != '\n') && (ch != EOF) ) { //loop through the line
			if(count == CUR_MAX) { // time to expand (for unexepectedly large line lengths) ?
				CUR_MAX *= 2; 
				count = 0;
				line = (char*)realloc(line, sizeof(char) * CUR_MAX); 
				tmpline = (char*)realloc(tmpline, sizeof(char) * CUR_MAX); 
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
				npolysites[k]=t;
				t=0;
				k++;
			}
			
			if (firstcontig==1) {
				firstcontig=0;
				t=0;
			}
			//Start preparing to read a new contig
			if ( k >= ncontigs_allocated ) {
				ncontigs_allocated+=NCONTIG_BATCH;
				if((contig=(char *)realloc(contig,sizeof(char)*ncontigs_allocated*namelen))==NULL) { 
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}	
				if((polysite=(int ***)realloc(polysite,sizeof(int **)*ncontigs_allocated))==NULL) {
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}		
				if((npolysites=(int *)realloc(npolysites,sizeof(int **)*ncontigs_allocated))==NULL) {
					fprintf(stderr,"error in memory allocation\n");
					exit(1);
				}		
			}
			if ( strlen(line) >= namelen ) {
				fprintf(stderr,"Error: a contig name is too long (more than %d characters) on line %d.\n",namelen-1,l);
				exit(1);
			}
			sscanf(line,">%s",&contig[k*namelen]);
			//			printf("%s\n",&contig[k*namelen]);
		}
		else { //line contains counts
				if (t==0) {
					if((polysite[k]=(int **)malloc(sizeof(int *)))==NULL){
						fprintf(stderr,"error in memory allocation\n");
						exit(1);
					}
					if((polysite[k][t]=(int *)malloc(sizeof(int)*7))==NULL){
						fprintf(stderr,"error in memory allocation\n");
						exit(1);
					}
				}
				else {
					if((polysite[k]=(int**)realloc(polysite[k],sizeof(int *)*(t+1)))==NULL){
						fprintf(stderr,"error in memory allocation\n");
						exit(1);
					}
					if((polysite[k][t]=(int *)malloc(sizeof(int)*7))==NULL){
						fprintf(stderr,"error in memory allocation\n");
						exit(1);
					}
				}
//				polysite[k][t][0]=t;
//				sscanf(line,"%*d\t%*f\t%d\t%d\t%d\t%d\t%d\t%d\t",&polysite[k][t][1],&polysite[k][t][2],&polysite[k][t][3],&polysite[k][t][4],&polysite[k][t][5],&polysite[k][t][6]);
//				if(chromosomes==XY){
//					sscanf(line,"%*d\t%*d\t%*f\t%d\t%d\t%d\t%d\t%d\t%d\t",&polysite[k][t][1],&polysite[k][t][2],&polysite[k][t][3],&polysite[k][t][4],&polysite[k][t][5],&polysite[k][t][6]);
//				}
//				else {
//					sscanf(line,"%*d\t%*d\t%*f\t%d\t%d\t%d\t%d\t%d\t%d\t",&polysite[k][t][4],&polysite[k][t][5],&polysite[k][t][6],&polysite[k][t][1],&polysite[k][t][2],&polysite[k][t][3]);
//				}
				if(chromosomes==XY){
					sscanf(line,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t",&polysite[k][t][0],&polysite[k][t][1],&polysite[k][t][2],&polysite[k][t][3],&polysite[k][t][4],&polysite[k][t][5],&polysite[k][t][6]);
				}
				else {
					sscanf(line,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t",&polysite[k][t][0],&polysite[k][t][4],&polysite[k][t][5],&polysite[k][t][6],&polysite[k][t][1],&polysite[k][t][2],&polysite[k][t][3]);
				}
				ni=0;
				for(i=1;i<=6;i++){
					ni+=polysite[k][t][i];
				}
				if(ni>nI){
					nI=ni;
				}
				t++;
		}
		memset(line, '\0', strlen(line)*(sizeof line));
	}
//	fprintf(stdout,"%d sites (individuals * positions), error rate (reads2snp) = %f\n",gensite,1.-totgenprob/gensite);
	
	npolysites[k]=t;
	ncontigs=k+1;
		
	//All reading has been done : clear memory.
	free(tmpline);
	free(line);
	
	*npolysites_p=npolysites;
	*polysite_p=polysite;
	*contig_p=contig;
	
	return ncontigs;
}

char **getgennames(const char *linein, int *nind)
{
	int ni=0,i=0,ichar,pipepos;
	char c;
	char *line = (char*)calloc(strlen(linein)+1,sizeof(char));
	char *tmpline = (char*)calloc(strlen(linein)+1,sizeof(char));
	char **name;

	strcpy(line,linein);
	
	//we'll have "position\tname1\tname2\t...\tnameN\0 : the number of tabs is the number of individuals.
	
	while ((c=line[i]) != '\0'){
		if (c == '\t') 
			ni++;
		i++;
	}
	if((name=(char **)calloc((size_t)(ni),sizeof(char *)))==NULL) { 
		fprintf(stderr,"getgennames: error in memory allocation\n");
		exit(1);
	}
	//find names iteratively ; shorten the line for each name found :
	sscanf(line,"position%[^\n]",tmpline);
	strcpy(line,tmpline);
	for (i=0; i<ni; i++){
		if((name[i]= (char *)malloc(sizeof(char) * NAME_LEN))==NULL) { 
			fprintf(stderr,"getgennames: error in memory allocation\n");
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
				fprintf(stderr,"getgennames: Error: one of the individual names seems to have reached the maximum length of %d characters (including prefixes).\n",NAME_LEN);
				exit(1);
			}
		}
		if(pipepos>0) {
			for(ichar=pipepos+1;ichar<NAME_LEN;ichar++) {
				name[i][ichar-pipepos-1]=name[i][ichar];
			}
		}
		
	}
	free(line);
	free(tmpline);
	*nind=ni;
	return name;
}

char **getvcfnames(const char *linein, int *nind)
{
	int ni=0,i=0,c,oldc;
	char *line = (char*)calloc(strlen(linein)+1,sizeof(char));
	char *tmpline = (char*)calloc(strlen(linein)+1,sizeof(char));
	char **name;
	
	strcpy(line,linein);
	
	oldc='a';
	while ((c=line[i]) != '\0'){
		if ( c == '\t' || c == ' ' ){
			if ( oldc == '\t' || oldc == ' '){
				i++;							
				continue;
			}						
			ni++;
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
	for (i=0; i<ni; i++){
		if((name[i]= (char *)malloc(sizeof(char) * NAME_LEN))==NULL) { 
			fprintf(stderr,"error in memory allocation\n");
			exit(37);
		}
		sscanf(line,"\t%s%[^\n]",name[i],tmpline);
		strcpy(line,tmpline);
	}
	free(line);
	free(tmpline);
	*nind=ni;
	return name;
}

int vcfsnp(const char *line,char *nuc, int *nnuc) //test if position is a snp or monomorphic; skip indels
{
	char allele1[NAME_LEN],allele2[NAME_LEN];
	int maxnuc=5; //A, T, G, C, or N
	int i,len;
	char c;
	int mnuc;
	
	//find out if this is a SNP, and what the observed alleles are
	sscanf(line,"%*s\t%*s\t%*s\t%s\t%s",allele1,allele2);
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
					return -1;
				}
				nuc[mnuc]=c;
				mnuc++;
			}
			else {
				return 1; //not a SNP
			}
			i=i+2;
		}
	}
	else { //strlen(allele2)==1
		if(strncmp(allele2,".",1)==0){
			mnuc=1;
		}
		else {
			nuc[1]=allele2[0];
			mnuc=2;
		}
	}
//	fprintf(stdout,"--%d--",mnuc);
	*nnuc=mnuc;
	return 0;
}

int findsex(char **name, int ni, char **femname, int nfem, char **malname, int nmal, int *sex, int *foundsex, int *ffound, int *nfgen, int *mfound, int *nmgen)
{
	//Find out the sex of the individuals.
	//sex has the sex of all individuals in the *gen or *vcf file (sex can be male, female, or absent)
	//foundsex has the sex of the individuals we want to study, in the order of the *gen or *vcf file
	
	int i,ifem,imal,nf,nm,nfound;
	nf=0;
	nm=0;
	for (i=0; i<ni; i++){			
		sex[i]=-1;
		for (ifem=0;ifem<nfem;ifem++) {
			if (strcmp(name[i],femname[ifem])==0) {
				sex[i]=FEMALE;
				foundsex[nm+nf]=FEMALE;
				ffound[ifem]=1;
				nf++;
			}
		}
		for (imal=0;imal<nmal;imal++) {
			if (strcmp(name[i],malname[imal])==0) {
				sex[i]=MALE;
				foundsex[nm+nf]=MALE;
				mfound[imal]=1;
				nm++;
			}
		}
	}
	nfound=nf+nm;
	*nfgen=nf;
	*nmgen=nm;
	return nfound;
}

std::vector<Genotype> vcfgenotypes(int ni, int *sex, const char *linein, char *nuc, int maxnuc)
{
	int i,c,ii;
	char vcfformatstring[NAME_LEN],genotypestring[NAME_LEN];
	char *line = (char*)calloc(strlen(linein)+1,sizeof(char));
	char *tmpline = (char*)calloc(strlen(linein)+1,sizeof(char));
	
	std::vector<Genotype> genotypes;
	strcpy(line,linein);
	
	//strip fields we don't use 
	sscanf(line,"%*s\t%*s\t%*s\t%*s\t%*s\t%*s\t%*s\t%*s\t%s%[^\n]",vcfformatstring,tmpline);
	strcpy(line,tmpline);
	if(strncmp(vcfformatstring,"GT",2)!=0){
		fprintf(stderr,"Error: this doesn't seem to be a supported format (\"GT\" tag not found where expected)\n");
		exit(1);					
	}
	
	for (i=0; i<ni; i++){
		if (sex[i]>=0) {
			Genotype genotype;
			sscanf(line,"\t%s%[^\n]",genotypestring,tmpline);
			for(ii=0;ii<2;ii++){
				c=genotypestring[2*ii];
				if(c=='.'){
					genotype.nucleotides[ii]='N';
					//									genotypes[j*nfound*2+2*fi+ii]='N';
					if(ii==0){
						genotypestring[2]='.';
					}
				}
				else if ( isdigit(c) && (c - '0') <= maxnuc ) {
					genotype.nucleotides[ii]=nuc[c - '0'];
					//									genotypes[j*nfound*2+2*fi+ii]=nuc[c - '0'];
				}
				else {
					fprintf(stderr,"Error in reading genotype\n");
					exit(1);
				}
			}
			genotypes.push_back(genotype);
		}
		else {
			sscanf(line,"\t%*s%[^\n]",tmpline);
		}
		strcpy(line,tmpline);
	}
	free(line);
	free(tmpline);
	return genotypes;
}

std::vector<Genotype> gengenotypes(int ni,int *sex,const char *linein)
{
	int i;
	char *line = (char*)calloc(strlen(linein)+1,sizeof(char));
	char *tmpline = (char*)calloc(strlen(linein)+1,sizeof(char));
	char nuc1,nuc2;
	
	std::vector<Genotype> genotypes;
	strcpy(line,linein);

	sscanf(line,"%*d%[^\n]",tmpline);
	strcpy(line,tmpline);
	for (i=0; i<ni; i++){
		if (sex[i]>=0) {
			Genotype genotype;
			sscanf(line,"\t%c%c|%*f%[^\n]",&nuc1,&nuc2,tmpline);
			genotype.nucleotides[0]=nuc1;
			genotype.nucleotides[1]=nuc2;
			genotypes.push_back(genotype);
		}
		else {
			sscanf(line,"\t%*c%*c|%*f%[^\n]",tmpline);
		}
		strcpy(line,tmpline);
	}
	free(line);
	free(tmpline);
	return genotypes;
}
