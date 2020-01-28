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
	double coverage;
	std::vector<SNP> snps; //number of snps
};

struct Individual {
	std::string name;
	int sex;
};

struct Genotype {
	char nucleotides[2];
	double probability;
};

struct Genotypes {
	int position;
	std::vector<Genotype> individualgenotypes;
};

struct GenotypeA {
	std::vector<int> allele;
	double probability;
};

struct GenotypesA {
	int position;
	std::vector<std::string> alleles;
	std::vector<GenotypeA> individualgenotypes;
};

struct ContigGenotypes {
	std::string name;
	std::vector<Individual> individuals;
	std::vector<Genotypes> genotypes;
};

struct ContigGenotypesA {
	std::string name;
	std::vector<Individual> individuals;
	std::vector<GenotypesA> genotypes;
};

struct Model {
	int haploid;
	int paralogs;
	int xy;
	int zw;
};

struct Varsite {
	int position;
	int genotypes_by_sex[_GenotypeCountsSize]; //size of enum GenotypeCounts
	std::vector<std::string> alleles;
};

struct ContigA {
	std::string name;
	double coverage;
	std::vector<Varsite> varsites; //number of snps
};
