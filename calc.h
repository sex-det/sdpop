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

#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_machine.h>
#include "types.h"

#define FEMALE 0
#define MALE 1
#define SEXES 2
//number of biological segregation types
enum SegregationType {
	J_AUTO, J_HAPLOID, J_PARA, J_HEMI, J_SEX, J_ZHEMI, J_ZW, JTYPES
};
enum DetailedSegregationTypes {
	JL_AUTO, JL_HAPLOID, JL_PARA1, JL_PARA2, JL_HEMI, JL_SEX1, JL_SEX2, JL_SEX3, JL_SEX4,
	JL_ZHEMI, JL_ZW1, JL_ZW2, JL_ZW3, JL_ZW4, JLTYPES
};
#define LTYPES 8 //number of XY segregation types
#define SITE 0
#define CONTIG 1
#define NONE 0
#define ERRORS 1
#define FIXED 2
#define READ_FROM_CNT_FILE 3
#define N11 1
#define N12 2
#define N22 3
#define XY 1
#define ZW 2
#define DIP 0 //diploid only
#define HAPDIP 1 //haploid or diploid
#define NOPARA 0
#define PARALOGS 1

template <typename F> void foreach_j (const Model model, const F & f) {
	f(J_AUTO);
	if(model.haploid){
		f(J_HAPLOID);
	}
	if(model.paralogs){
		f(J_PARA);
	}
	if(model.xy){
		f(J_HEMI);
		f(J_SEX);
	}
	if(model.zw){
		f(J_ZHEMI);
		f(J_ZW);
	}
}

template <typename F> void foreach_jl (const Model model, const F & f) {
	f(JL_AUTO);
	if(model.haploid){
		f(JL_HAPLOID);
	}
	if(model.paralogs){
		f(JL_PARA1);
		f(JL_PARA2);
	}
	if(model.xy){
		f(JL_HEMI);
		f(JL_SEX1);
		f(JL_SEX2);
		f(JL_SEX3);
		f(JL_SEX4);
	}
	if(model.zw){
		f(JL_ZHEMI);
		f(JL_ZW1);
		f(JL_ZW2);
		f(JL_ZW3);
		f(JL_ZW4);
	}
}

template <typename F> void foreach_l_para (const Model model, const F & f) {
	if(model.paralogs){
		f(0,JL_PARA1,J_PARA);
		f(1,JL_PARA2,J_PARA);		
	}
}
template <typename F> void foreach_l_xy (const Model model, const F & f) {
	if(model.xy){
		f(0,JL_SEX1,J_SEX);
		f(1,JL_SEX2,J_SEX);
		f(2,JL_SEX3,J_SEX);
		f(3,JL_SEX4,J_SEX);
	}
}
template <typename F> void foreach_l_zw (const Model model, const F & f) {
	if(model.zw){
		f(0,JL_ZW1,J_ZW);
		f(1,JL_ZW2,J_ZW);
		f(2,JL_ZW3,J_ZW);
		f(3,JL_ZW4,J_ZW);
	}
}

double intpow(double x,int y);
double* vectorize_d(double x0, double x1, double x2);
unsigned int* vectorize_ui(unsigned int x0, unsigned int x1, unsigned int x2);
double* errormult(double ematrix[3][3], double x[3]);
//void CondSiteProbs(int ncontigs, int *npolysites, int ***polysite, double Q[3][3], double *****P, long double ***condsiteprob);
//void CondSegProbs(int ncontigs, int *npolysites, double *rho, long double ***condsiteprob, long double ***condsegprob);
void CondSiteProbs(std::vector<ContigA>& contigs, const Model model, double Q[3][3], double *****P, long double ***condsiteprob);
void CondSegProbs(std::vector<ContigA>& contigs, const Model model, double **rho, long double ***condsiteprob, long double ***condsegprob);
void horner(int length, int lmax, long double* array, double* param, double* result);
void loghorner(int length, int lmax, long double* logarray, double* param, double* result);
void calcQ(double m33[3][3], double e0, double e1, double e2);
double zscore(int k, int npolysites, int **polysite, double *rho, double Q[3][3], int lmax, int iterations, double ****P, long double **condsegprob);
