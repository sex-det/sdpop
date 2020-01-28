reading.o: reading.cpp reading.h types.h
	$(CXX) -Wall -c reading.cpp -o reading.o
calc.o: calc.cpp calc.h types.h
	$(CXX) -Wall -c calc.cpp

sdpop: sdpop.cpp reading.o calc.o types.h
	$(CXX) -Wall -o sdpop -lm reading.o calc.o sdpop.cpp
popsum: popsum.cpp reading.o types.h
	$(CXX) -Wall -o popsum -lm reading.o popsum.cpp
wxyz_genotyper: wxyz_genotyper.cpp reading.o calc.o types.h
	$(CXX) -Wall -o wxyz_genotyper -lm reading.o wxyz_genotyper.cpp

