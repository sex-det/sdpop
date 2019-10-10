reading.o: reading.cpp reading.h types.h
	$(CXX) -Wall -c reading.cpp -o reading.o
calc.o: calc.cpp calc.h types.h
	$(CXX) -Wall -c calc.cpp

sdpop_2018: sdpop_2018.cpp reading.o calc.o types.h
	$(CXX) -Wall -o sdpop_2018 -lm reading.o calc.o sdpop_2018.cpp
popsum: popsum.cpp reading.o types.h
	$(CXX) -Wall -o popsum -lm reading.o popsum.cpp
wxyz_genotyper: wxyz_genotyper.c
	$(CC) -Wall -o wxyz_genotyper -lm wxyz_genotyper.c

