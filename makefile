sdpop: sdpop.c
	gcc -Wall -o sdpop -lm sdpop.c
sdpop_2018: sdpop_2018.c reading.o calc.o
	gcc -Wall -o sdpop_2018 -lm reading.o calc.o sdpop_2018.c
sdpop_hap: sdpop_hap.c reading.o calc.o
	gcc -Wall -o sdpop_hap -lm reading.o calc.o sdpop_hap.c
sdpop_perm: sdpop_perm.c
	gcc -Wall -o sdpop_perm -lm sdpop_perm.c
sdpop_nonpar: sdpop_nonpar.c
	gcc -Wall -lm -o sdpop_nonpar sdpop_nonpar.c
popsum: popsum.c reading.o 
	gcc -Wall -o popsum -lm reading.o popsum.c
readgen: readgen.c
	gcc -Wall -o readgen -lgsl -lgslcblas -lm readgen.c
ms_simul: ms_simul.c
	gcc -Wall -o ms_simul -lm ms_simul.c
simul: simul.c
	gcc -Wall -o simul -lgsl -lgslcblas -lm simul.c
eer: eer.c
	gcc -Wall -o eer eer.c
heterozygosity: heterozygosity.c
	gcc -Wall -o heterozygosity -lm heterozygosity.c
vcf2gen: vcf2gen.c
	gcc -Wall -o vcf2gen vcf2gen.c

reading.o: reading.c reading.h
	gcc -Wall -c reading.c
calc.o: calc.c calc.h
	gcc -Wall -c calc.c

wxyz_genotyper: wxyz_genotyper.c reading.o
	gcc -Wall -o wxyz_genotyper -lm reading.o wxyz_genotyper.c
fisher: fisher.c reading.o calc.o
	gcc -Wall -o fisher -lm reading.o calc.o fisher.c
gen2fas: gen2fas.c reading.o
	gcc -Wall -o gen2fas reading.o gen2fas.c

