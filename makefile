

sdpop: sdpop.c
	$(CC) -Wall -o sdpop -lm sdpop.c
sdpop_2018: sdpop_2018.cpp reading.o calc.o
	$(CXX) -Wall -o sdpop_2018 -lm reading.o calc.o sdpop_2018.cpp
sdpop_hap: sdpop_hap.c reading.o calc.o
	$(CC) -Wall -o sdpop_hap -lm reading.o calc.o sdpop_hap.c
sdpop_perm: sdpop_perm.c
	$(CC) -Wall -o sdpop_perm -lm sdpop_perm.c
sdpop_nonpar: sdpop_nonpar.c
	$(CC) -Wall -lm -o sdpop_nonpar sdpop_nonpar.c
popsum: popsum.cpp reading.o 
	$(CXX) -Wall -o popsum -lm reading.o popsum.cpp
readgen: readgen.c
	$(CC) -Wall -o readgen -lgsl -lgslcblas -lm readgen.c
ms_simul: ms_simul.c
	$(CC) -Wall -o ms_simul -lm ms_simul.c
simul: simul.c
	$(CC) -Wall -o simul -lgsl -lgslcblas -lm simul.c
eer: eer.c
	$(CC) -Wall -o eer eer.c
heterozygosity: heterozygosity.c
	$(CC) -Wall -o heterozygosity -lm heterozygosity.c
postrecalc: postrecalc.c
	$(CC) -Wall -o postrecalc -lm postrecalc.c
vcf2gen: vcf2gen.c
	$(CC) -Wall -o vcf2gen vcf2gen.c

reading.o: reading.cpp reading.h
	$(CXX) -Wall -c reading.cpp
calc.o: calc.cpp calc.h
	$(CXX) -Wall -c calc.cpp

wxyz_genotyper: wxyz_genotyper.c reading.o
	$(CC) -Wall -o wxyz_genotyper -lm reading.o wxyz_genotyper.c
fisher: fisher.c reading.o calc.o
	$(CC) -Wall -o fisher -lm reading.o calc.o fisher.c
gen2fas: gen2fas.c reading.o
	$(CC) -Wall -o gen2fas reading.o gen2fas.c
vcfhet: reading.o vcfhet.cpp
	$(CXX) -Wall -o vcfhet -lm reading.o vcfhet.cpp
