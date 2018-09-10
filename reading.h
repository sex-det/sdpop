#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define XY 1
#define ZW 2
#define FEMALE 0
#define MALE 1
#define SEXES 2

int DNA2int(char c);
char int2DNA(int i);
int read_cnt(FILE *fp, int namelen, int chromosomes, char **contig_p, int **npolysites_p, int ****polysite_p);
int read_cnt2(FILE *fp, int namelen, int chromosomes, char **contig_p, int **npolysites_p, int ****polysite_p);
char **getgennames(char *line, int *ni);
int findsex(char **name, int ni, char **femname, int nfem, char **malname, int nmal, int *sex, int *foundsex, int *ffound, int *nfgen, int *mfound, int *nmgen);
