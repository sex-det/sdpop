#pragma once

#include <string>
#include <vector>


enum GenotypeCounts {
	N11F, N12F, N22F, N11M, N12M, N22M, _GenotypeCountsSize
};

struct SNP {
	int position;
	int genotypes_by_sex[_GenotypeCountsSize]; //size of enum GenotypeCounts
	int alleles[2];
};

struct Contig {
	std::string name;
	std::vector<SNP> snps; //number of snps
};

struct Individual {
	std::string name;
	int sex;
};

struct Genotype {
	char nucleotides[2];
};

struct Genotypes {
	int position;
	std::vector<Genotype> individualgenotypes;
};

struct ContigGenotypes {
	std::string name;
	std::vector<Individual> individuals;
	std::vector<Genotypes> genotypes;
};

struct Model {
	int haploid;
	int paralogs;
	int xy;
	int zw;
};
