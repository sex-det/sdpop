#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int npoly=0;
int e0t=0,e1t=0,e2t=0;
int genformat=0;

void countgenotypes(FILE *fp,int **sequence, int npos, int nind, int contiglength, double e[3]) {
	int n11m,n12m,n22m,n11f,n12f,n22f;
	int t,i,nxdiff,nydiff,nucfix,nfix,ny0,nx0;
	double a,fa,fx,fy,pi_a=0,pi_x=0,pi_y=0,D_xy=0,pif_a=0,pif_x=0,pif_y=0,Df_xy=0;
	int e0,e1,e2;
	
	//	Watterson's theta: autosomes
	a=0;
	for(i=1;i<4*nind;i++){
		a+=1./(double)i;
	}
	fprintf(fp,"th_WA = %f ",((double)npos/(double)contiglength)/a);
	//	Watterson's theta: count X and Y polymorphisms and differences
	nydiff=0;	
	nxdiff=0;
	nfix=0;
	for(t=0;t<npos;t++){
		ny0=0;
		nx0=0;
		//Y
		for(i=0;i<nind;i++){
			if(sequence[i][t]==0){
				ny0++;
			}
		}
		if(ny0 > 0 && ny0 < nind){
			nydiff++;
		}
		 //X
		for(i=nind;i<4*nind;i++){
			if(sequence[i][t]==0){
				nx0++;
			}
		}
		if(nx0 > 0 && nx0 < 3*nind){
			nxdiff++;
		}
		if( (nx0==0 && ny0==nind) || (ny0==0 && nx0==3*nind) ){
			nfix++;
		}
		pi_a+=(double)((nx0+ny0)*(3*nind-nx0+nind-ny0))/(0.5*(4*nind)*(4*nind-1));
		pi_y+=(double)((ny0)*(nind-ny0))/(0.5*(nind)*(nind-1));
		pi_x+=(double)((nx0)*(3*nind-nx0))/(0.5*(3*nind)*(3*nind-1));
		D_xy+=(double)(nx0*(nind-ny0)+ny0*(3*nind-nx0))/(double)(nind*3*nind);

		fx=(double)nx0/(3.*nind);
		fy=(double)ny0/(1.*nind);
		fa=(double)(nx0+ny0)/(4.*nind);
		pif_a+=2.*fa*(1-fa);
		pif_y+=2.*fy*(1-fy);
		pif_x+=2.*fx*(1-fx);
		Df_xy+=fx*(1-fy)+fy*(1-fx);
	}
 	// output Watterson's theta for X
	a=0;
	for(i=1;i<3*nind;i++){
		a+=1./(double)i;
	}
	fprintf(fp,"th_WX = %f ",((double)nxdiff/(double)contiglength)/a);
 	// output Watterson's theta for Y
	a=0;
	for(i=1;i<nind;i++){
		a+=1./(double)i;
	}
	fprintf(fp,"th_WY = %f\t",((double)nydiff/(double)contiglength)/a);
 	// output fixed differences 
 	fprintf(fp,"fixed = %f\t",(double)nfix/(double)contiglength);
 	fprintf(fp,"pi_a = %f pi_x = %f pi_y = %f D_xy = %f\t",pi_a/contiglength,pi_x/contiglength,pi_y/contiglength,D_xy/contiglength);
 	fprintf(fp,"pif_a = %f pif_x = %f pif_y = %f Df_xy = %f\n",pif_a/contiglength,pif_x/contiglength,pif_y/contiglength,Df_xy/contiglength);
	
	if(genformat){
		fprintf(fp,"position");
        for(i=0;i<nind;i++){ //males
        	fprintf(fp,"\tsp|male%d",i);
        }
        for(i=0;i<nind;i++){ //females
        	fprintf(fp,"\tsp|female%d",i);
        }
        fprintf(fp,"\n");
        
        //output gen format
        for(t=0;t<contiglength;t++){
        	fprintf(fp,"%d",t+1);
        	for(i=0;i<nind;i++){ //males
        		if(t<npos){ //"real" polymorphic positions 
        			if(sequence[i][t]!=sequence[i+2*nind][t]){ //heterozygous
        				if((double)rand()/RAND_MAX > 2.*e[1]) { //no errors
        					fprintf(fp,"\tAT|1");
        				}
        				else if ((double)rand()/RAND_MAX < 0.5) { //errors
        					fprintf(fp,"\tAA|1");
        				}
        				else {
        					fprintf(fp,"\tTT|1");
        				}
        			}
        			else if(sequence[i][t]==0){ //homozygous for allele 1
        				if((double)rand()/RAND_MAX > e[0]+e[2]) { //no errors
        					fprintf(fp,"\tAA|1");
        				}
        				else if ((double)rand()/RAND_MAX < e[0]/(e[0]+e[2])) { //errors
        					fprintf(fp,"\tAT|1");
        				}
        				else {
       					fprintf(fp,"\tTT|1");
         				}
        			}
        			else { //homozygous for allele 2
        				if((double)rand()/RAND_MAX > e[0]+e[2]) { //no errors
       					fprintf(fp,"\tTT|1");
         				}
        				else if ((double)rand()/RAND_MAX < e[0]/(e[0]+e[2])) { //errors
        					fprintf(fp,"\tAT|1");
        				}
        				else {
        					fprintf(fp,"\tAA|1");
        				}
        			}
        		}
        		else { //no polymorphism observed: add errors
        			if((double)rand()/RAND_MAX > e[0]+e[2]){ //no errors
        					fprintf(fp,"\tAA|1");
        			}
        			else if ((double)rand()/RAND_MAX < e[0]/(e[0]+e[2])) { //errors
        					fprintf(fp,"\tAT|1");
        			}
        			else {
       					fprintf(fp,"\tTT|1");
         			}
        		}
        	}
        	for(i=nind;i<2*nind;i++){ //females
        		if(t<npos){ //"real" polymorphic positions 
        			if(sequence[i][t]!=sequence[i+2*nind][t]){
        				if((double)rand()/RAND_MAX > 2.*e[1]) { //no errors
        					fprintf(fp,"\tAT|1");
        				}
        				else if ((double)rand()/RAND_MAX < 0.5) { //errors
        					fprintf(fp,"\tAA|1");
        				}
        				else {
       					fprintf(fp,"\tTT|1");
         				}
        			}
        			else if(sequence[i][t]==0){
        				if((double)rand()/RAND_MAX > e[0]+e[2]) { //no errors
        					fprintf(fp,"\tAA|1");
        				}
        				else if ((double)rand()/RAND_MAX < e[0]/(e[0]+e[2])) { //errors
        					fprintf(fp,"\tAT|1");
        				}
        				else {
       					fprintf(fp,"\tTT|1");
         				}
        			}
        			else {
        				if((double)rand()/RAND_MAX > e[0]+e[2]) { //no errors
        					fprintf(fp,"\tTT|1");
        				}
        				else if ((double)rand()/RAND_MAX < e[0]/(e[0]+e[2])) { //errors
        					fprintf(fp,"\tAT|1");
        				}
        				else {
        					fprintf(fp,"\tAA|1");
        				}
        			}
        		}
        		else { //no polymorphism observed: add errors
        			if((double)rand()/RAND_MAX > e[0]+e[2]){ //no errors
        					fprintf(fp,"\tAA|1");
        			}
        			else if ((double)rand()/RAND_MAX < e[0]/(e[0]+e[2])) { //errors
        					fprintf(fp,"\tAT|1");
        			}
        			else {
       					fprintf(fp,"\tTT|1");
         			}
        		}
        	}
        		fprintf(fp,"\n");
        }
    }
	else {
	//output cnt format
	for(t=0;t<contiglength;t++){
		n11m=n12m=n22m=0;
		n11f=n12f=n22f=0;
		e0=e1=e2=0;
		for(i=0;i<nind;i++){ //males
			if(t<npos){ //"real" polymorphic positions 
				if(sequence[i][t]!=sequence[i+2*nind][t]){ //heterozygous
					if((double)rand()/RAND_MAX > 2.*e[1]) { //no errors
						n12m++;
					}
					else if ((double)rand()/RAND_MAX < 0.5) { //errors
						n11m++;
						e1++;
					}
					else {
						n22m++;
						e1++;
					}
				}
				else if(sequence[i][t]==0){ //homozygous for allele 1
					if((double)rand()/RAND_MAX > e[0]+e[2]) { //no errors
						n11m++;
					}
					else if ((double)rand()/RAND_MAX < e[0]/(e[0]+e[2])) { //errors
						n12m++;
						e0++;
					}
					else {
						n22m++;
						e2++;
					}
				}
				else { //homozygous for allele 2
					if((double)rand()/RAND_MAX > e[0]+e[2]) { //no errors
						n22m++;
					}
					else if ((double)rand()/RAND_MAX < e[0]/(e[0]+e[2])) { //errors
						n12m++;
						e0++;
					}
					else {
						n11m++;
						e2++;
					}
				}
			}
			else { //no polymorphism observed: add errors
				if((double)rand()/RAND_MAX > e[0]+e[2]){ //no errors
					n11m++;
				}
				else if ((double)rand()/RAND_MAX < e[0]/(e[0]+e[2])) { //errors
					n12m++;
					e0++;
				}
				else {
					n22m++;
					e2++;
				}
			}
		}
		for(i=nind;i<2*nind;i++){ //females
			if(t<npos){ //"real" polymorphic positions 
				if(sequence[i][t]!=sequence[i+2*nind][t]){
					if((double)rand()/RAND_MAX > 2.*e[1]) { //no errors
						n12f++;
					}
					else if ((double)rand()/RAND_MAX < 0.5) { //errors
						n11f++;
						e1++;
					}
					else {
						n22f++;
						e1++;
					}
				}
				else if(sequence[i][t]==0){
					if((double)rand()/RAND_MAX > e[0]+e[2]) { //no errors
						n11f++;
					}
					else if ((double)rand()/RAND_MAX < e[0]/(e[0]+e[2])) { //errors
						n12f++;
						e0++;
					}
					else {
						n22f++;
						e2++;
					}
				}
				else {
					if((double)rand()/RAND_MAX > e[0]+e[2]) { //no errors
						n22f++;
					}
					else if ((double)rand()/RAND_MAX < e[0]/(e[0]+e[2])) { //errors
						n12f++;
						e0++;
					}
					else {
						n11f++;
						e2++;
					}
				}
			}
			else { //no polymorphism observed: add errors
				if((double)rand()/RAND_MAX > e[0]+e[2]){ //no errors
					n11f++;
				}
				else if ((double)rand()/RAND_MAX < e[0]/(e[0]+e[2])) { //errors
					n12f++;
					e0++;
				}
				else {
					n22f++;
					e2++;
				}
			}
		}
		if((n11m==nind && n11f==nind) || (n22m==nind && n22f==nind)) { //check for polymorphism
			continue;
		}
		else {
			e0t+=e0;
			e1t+=e1;
			e2t+=e2;
			npoly++;
			fprintf(fp,"%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n",t,"AT",n11f,n12f,n22f,n11m,n12m,n22m);
		}
	}
	}
}

int main (int argc, char *argv[]) {
	int contiglength=100,linelength;
	int ncontigs=1000;
	int nind=7;
	int nchar,npos,firstcontig;
	int i,ii,k,t;
	FILE *inputf,*outputf;	
	char *line=NULL;
	char contigname[1000];
	char c;
	int **sequence;
	double e[3],xyfrac,xhemifrac;
	
	for(i=0;i<argc;i++) {
		fprintf(stdout,"%s ",argv[i]);
	}
	fprintf(stdout,"\n");
	if (argc != 8) {
		fprintf(stdout,"Usage: %s infile outfile e0 e1 e2 x-hemi genformat\n",argv[0]);
		exit(1);
	}

	if((inputf=fopen(argv[1],"r"))==NULL){
		fprintf(stderr,"error opening input file %s\n",argv[1]);
		exit(1);
	}
	fscanf(inputf,"nsex=%d ncontigs=%d xyfrac=%lf pi=%*f length=%d time=%*f",&nind,&ncontigs,&xyfrac,&contiglength);
	if((outputf=fopen(argv[2],"w"))==NULL){
		fprintf(stderr,"error opening output file %s\n",argv[2]);
		exit(1);
	}
	
	e[0]=atof(argv[3]);
	e[1]=atof(argv[4]);
	e[2]=atof(argv[5]);
	xhemifrac=atof(argv[6]);
	genformat=atoi(argv[7]);
	
	if((sequence=(int **)malloc(sizeof(int *)*nind*2*2))==NULL) {
		fprintf(stderr,"error in memory allocation\n");
		exit(1);
	}
	for(ii=0;ii<nind*2*2;ii++){
		if((sequence[ii]=(int *)malloc(sizeof(int)*contiglength))==NULL) {
			fprintf(stderr,"error in memory allocation\n");
			exit(1);
		}	
	}

	linelength=(1000>contiglength) ? 1000 : contiglength;
	line=(char *)malloc(sizeof(char)*linelength);
	firstcontig=1;
	k=0;
	ii=0;

	while((nchar=getline(&line, (size_t *) &linelength, inputf))>0){
		if(line[0] == '>' ){
			if(!firstcontig && ii>0){ // treat last read contig		
				if(k>ncontigs*(1.-xyfrac) && k<ncontigs*(1.-(1.-xhemifrac)*xyfrac)){ //x-hemizygotes
					fprintf(outputf,">h%s\t",contigname);
				}
				else {
					fprintf(outputf,">%s\t",contigname);
				}
				if(ii!=2*2*nind){
					fprintf(stderr,"error: expected %d sequences, read %d\n",2*2*nind,ii);
					exit(1);
				}
				if(k>ncontigs*(1.-xyfrac) && k<ncontigs*(1.-(1.-xhemifrac)*xyfrac)){ //x-hemizygotes
					for(i=0;i<nind;i++){ //males
						for(t=0;t<contiglength;t++){
							sequence[i][t]=sequence[i+2*nind][t];
						}
					}
				}	
				countgenotypes(outputf,sequence,npos,nind,contiglength,e);
			}
			sscanf(line,">%s",contigname);
			firstcontig=0;
//			if(k<ncontigs*(1.-xyfrac)){
//				fprintf(outputf,">%s\t%d\t%f\t0\n",contigname,k,(double)ncontigs*(1.-xyfrac));
//			}
//			else if(k<ncontigs*(1.-(1.-xhemifrac)*xyfrac)){ //x-hemizygotes
//				fprintf(outputf,">%s\t%d\t1\n",contigname,k);
//			}
//			else {
//				fprintf(outputf,">%s\t%d\t2\n",contigname,k);
//			}
			k++;
			ii=0;				
		}
		else {
			npos=nchar-1; //number of polymorphic sites
			for(t=0;t<npos;t++){
				c=line[t];
				sequence[ii][t]=atoi(&c);
			}
			ii++;
		}
	}
	if(ii>0){ // treat last read contig		
		if(k>ncontigs*(1.-xyfrac) && k<ncontigs*(1.-(1.-xhemifrac)*xyfrac)){ //x-hemizygotes
			fprintf(outputf,">h%s\t",contigname);
		}
		else {
			fprintf(outputf,">%s\t",contigname);
		}
		if(ii!=2*2*nind){
			fprintf(stderr,"error: expected %d sequences, read %d\n",2*2*nind,ii);
			exit(1);
		}
		if(k>ncontigs*(1.-xyfrac) && k<ncontigs*(1.-(1.-xhemifrac)*xyfrac)){ //x-hemizygotes
			for(i=0;i<nind;i++){ //males
				for(t=0;t<contiglength;t++){
					sequence[i][t]=sequence[i+2*nind][t];
				}
			}
		}	
		countgenotypes(outputf,sequence,npos,nind,contiglength,e);
	}
	fprintf(stdout,"npoly=%d (%d genotypes), e0=%f, e1=%f, e2=%f\n",npoly,npoly*nind*2,(double)e0t/(npoly*nind*2),(double)e1t/(npoly*nind*2),(double)e2t/(npoly*nind*2));
	
	free(line);	
	for(i=0;i<nind;i++){
	free(sequence[i]);
	}
	free(sequence);
	
	return 0;
}
