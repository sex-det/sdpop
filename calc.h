#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_machine.h>

#define FEMALE 0
#define MALE 1
#define SEXES 2
//number of biological segregation types
enum SegregationType {
	J_AUTO, J_HAPLOID, J_PARA, J_HEMI, J_SEX, JTYPES
};
#define JLTYPES 9 //number of detailed segregation types
#define JL_AUTO 0
#define JL_HAPLOID 1
#define JL_PARA1 2
#define JL_PARA2 3
#define JL_HEMI 4
#define JL_SEX1 5
#define JL_SEX2 6
#define JL_SEX3 7
#define JL_SEX4 8
#define LTYPES 4 //number of XY segregation types
#define SITE 0
#define CONTIG 1
#define NONE 0
#define ERRORS 1
#define FIXED 2
#define N11 1
#define N12 2
#define N22 3
#define N11F 1
#define N12F 2
#define N22F 3
#define N11M 4
#define N12M 5
#define N22M 6
#define XY 1
#define ZW 2
#define DIP 0 //diploid only
#define HAPDIP 1 //haploid or diploid
#define NOPARA 0
#define PARALOGS 1

double intpow(double x,int y);
double* vectorize_d(double x0, double x1, double x2);
unsigned int* vectorize_ui(unsigned int x0, unsigned int x1, unsigned int x2);
double* errormult(double ematrix[3][3], double x[3]);
void CondSiteProbs(int ncontigs, int *npolysites, int ***polysite, double Q[3][3], double *****P, long double ***condsiteprob);
void CondSegProbs(int ncontigs, int *npolysites, double *rho, long double ***condsiteprob, long double ***condsegprob);
void horner(int length, int lmax, long double* array, double* param, double* result);
void loghorner(int length, int lmax, long double* logarray, double* param, double* result);
void calcQ(double m33[3][3], double e0, double e1, double e2);
double zscore(int k, int npolysites, int **polysite, double *rho, double Q[3][3], int lmax, int iterations, double ****P, long double **condsegprob);
