#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include <vector>
#include <cstdio>
#include "types.h"

#define NONE 0
#define XY 1
#define ZW 2
#define FEMALE 0
#define MALE 1
#define SEXES 2

int DNA2int(char c);
char int2DNA(int i);
int read_cnt(FILE *fp, int namelen, int chromosomes, char **contig_p, int **npolysites_p, int ****polysite_p);
//int read_cnt2(FILE *fp, int namelen, int chromosomes, char **contig_p, int **npolysites_p, int ****polysite_p);
//int read_cnt2(FILE *fp, int namelen, int chromosomes, std::vector<std::string>& contig_p, int **npolysites_p, int ****polysite_p);
int read_cnt_model(FILE *fp, int namelen, const Model model, std::vector<ContigA>& contigs);
int read_cnt_model_error(FILE *fp, int namelen, const Model model, std::vector<ContigA>& contigs, double *e);
int read_cnt2(FILE *fp, int namelen, int chromosomes, std::vector<ContigA>& contigs);
int read_cnt3(FILE *fp, int namelen, int chromosomes, std::vector<ContigA>& contigs, double *e);
char **getgennames(const char *linein, int *ni);
char **getvcfnames(const char *linein, int *nind);
int vcfsnp(const char *line,char *nuc, int *nnuc);
std::vector<std::string> vcfsnpindel(const char *line);
int findsex(char **name, int ni, char **femname, int nfem, char **malname, int nmal, int *sex, int *foundsex, int *ffound, int *nfgen, int *mfound, int *nmgen);
std::vector<Genotype> vcfgenotypes(int ni, int *sex, const char *linein, char *nuc, int maxnuc);
std::vector<GenotypeA> vcfalleles(int ni, int *sex, const char *linein);
GenotypesA gengenotypes(int pos, int ni,int *sex, const char *linein);
