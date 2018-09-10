#include <stdlib.h>
#include <stdio.h>

int main (int argc, char *argv[]) {

	int MAXLINE=300;
	char line[MAXLINE],contigname[100];
	int i,k,t,j,error,pos, n11f, n12f, n22f,n11m,n12m,n22m,S,type,npos,kerror;
	float E[4],f[4],f_sim,totprob;
	float pi[4];
	int falsepos[4][100],truepos[4][100],ntype[4],ktype[4],falseclass[100],trueclass[100];
	int kfp[4],ktp[4],kfn[4];
	FILE *fp,*simfile,*outfile;
	int simerror,totalsimerrors=0;
	
	if((fp=fopen(argv[1],"r"))==NULL){
		fprintf(stderr,"error opening input file %s\n",argv[1]);
		exit(1);
	}
	if((simfile=fopen(argv[2],"r"))==NULL){
		fprintf(stderr,"error opening input file %s\n",argv[1]);
		exit(1);
	}
	if((outfile=fopen(argv[3],"w"))==NULL){
		fprintf(stderr,"error opening output file %s\n",argv[2]);
		exit(1);
	}

	for(i=0;i<4;i++) {
		ntype[i]=0;
		ktype[i]=0;
		kfp[i]=0;
		kfn[i]=0;
		ktp[i]=0;
		for (j=0;j<100;j++) {
			falsepos[i][j]=0;
			truepos[i][j]=0;
			falseclass[j]=0;
			trueclass[j]=0;
		}
	}
	
	i=0;
	k=0;
	t=0;
	error=0;
	kerror=0;
	totprob=0;
	while (fgets(line,MAXLINE,fp)!=NULL) {
//		if (i==0) { //treat the first line of the EM output file
//			sscanf(line,"log-likelihood: %*f, pi:\t%f\t%f\t%f\t%f",&pi[0],&pi[1],&pi[2],&pi[4]);
//			i++;
//		}
//		else
		if (line[0] == '>' ) {
//			sscanf(line,">%s\t%d\tj : %*d; pp : %f\tj : %*d; pp : %f\tj : %*d; pp : %f\tj : %*d; pp : %f\t%d",contigname,&npos,&E[0],&E[1],&E[2],&E[3],&S);
			sscanf(line,">%s\t%d\tl : %*d; pp : %f\tl : %*d; pp : %f\tl : %*d; pp : %f\t%d",contigname,&npos,&E[0],&E[1],&E[2],&S);
//			printf("%s",line);
//			printf("%s\t%d\t%f\t%f\t%f\t%d\n",contigname,npos,E[0],E[1],E[2],S);
			fgets(line,MAXLINE,simfile);
			sscanf(line,">%*s\t%d",&type);
//			printf("%s\n",line);
			k++;
			ktype[type]++;
			if (type != S) {
				kerror++;
				kfp[S]++;
				kfn[type]++;
//				printf("%d\t%d\n",S,type);
			}
			else {
				ktp[S]++;
				//no count for true negatives
			}
		}
		else {
			sscanf(line,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t(%f)\t%f\t(%f)\t%f\t(%f)\t%f\t(%f)\t%d", &pos, &n11f, &n12f, &n22f,&n11m,&n12m,&n22m,&E[0],&f[0],&E[1],&f[1],&E[2],&f[2],&E[3],&f[3],&S);
			fgets(line,MAXLINE,simfile);
			sscanf(line,"%d %f %d %d %d %d %d %d %d",&type,&f_sim,&n11f,&n12f,&n22f,&n11m,&n12m,&n22m,&simerror);
			t++;
			totalsimerrors+=simerror;
//			ntype[type]++;
//			j=(int)(E[S]*100);
//			if(j-1<=0) {
//				fprintf(stderr,"j too small : %d\n",j);
//				exit(1);
//			}
//			if (type != S) {
//				falsepos[S][j-1]++;
//				falseclass[j-1]++;
//			}
//			else {
//				truepos[S][j-1]++;
//				trueclass[j-1]++;
//			}
//			if (type != S) {
//				//fprintf(outfile,"%d %d %d %d %d %d %d %f %f %f %f %f %f %f %f %d %d %f\n",pos,n11f,n12f,n22f,n11m,n12m,n22m,E[0],f[0],E[1],f[1],E[2],f[2],E[3],f[3],S,type,f_sim);
//				error++;
//			}
//			totprob+=E[S];
		}
	}
//	printf("%s %f %f %f %f %d %d %d %d %f\n",argv[1],pi[0],pi[1],pi[2],pi[4],error,t,kerror,k,totprob/t);
//	printf("%s %f %f %f %f %d %d %d\n",argv[1],pi[0],pi[1],pi[2],pi[4],kerror,k,totalsimerrors);
//	printf("%s %f %f %f %f ",argv[1],pi[0],pi[1],pi[2],pi[4]);
	for (i=0;i<3;i++) {
	printf("%d %f %f ",ktype[i],(double)kfp[i]/(double)(k-ktype[i]),(double)ktp[i]/(double)ktype[i]);
	}
	printf("%d %d %d\n",kerror,k,totalsimerrors);
	//cumulative counts
//	for (j=98;j>=0;j--) {
//		falseclass[j]+=falseclass[j+1];
//		trueclass[j]+=trueclass[j+1];
//	}
//	for(i=0;i<4;i++) {
//		for (j=98;j>=0;j--) {
//			falsepos[i][j]+=falsepos[i][j+1];
//			truepos[i][j]+=truepos[i][j+1];
//		}
//	//		printf("%d : count = %d, true + false positives = %d\n",i,ntype[i],falsepos[i][0]+truepos[i][0]);
//	}
	//false positive rate : #false positives / #negatives
	//true positive rate : #true positives / #positives
//	for (j=0;j<100;j++) {
//		fprintf(outfile,"%f %f %f ",(float)j/100.,(float)falseclass[j]/(float)t,(float)trueclass[j]/(float)t);
//		for(i=0;i<4;i++) {
//			fprintf(outfile,"%f %f ",(float)falsepos[i][j]/(float)(t-ntype[i]),(float)truepos[i][j]/(float)ntype[i]);
//		}
//		fprintf(outfile,"\n");
//	}

	fclose(fp);
	fclose(simfile);
	fclose(outfile);
	return 0;
}