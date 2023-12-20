#for static compilation, on computer clusters for instance, try:
#make CXXFLAGS="--static" all

#uncomment next line for debugging
#CXXFLAGS += -ggdb

all: popsum sdpop wxyz_genotyper

reading.o: reading.cpp reading.h types.h
	$(CXX) $(CXXFLAGS) -Wall -c reading.cpp -o reading.o
calc.o: calc.cpp calc.h types.h
	$(CXX) $(CXXFLAGS) -Wall -c calc.cpp

sdpop: sdpop.cpp reading.o calc.o types.h
	$(CXX) $(CXXFLAGS) -Wall -o sdpop -lm reading.o calc.o sdpop.cpp
popsum: popsum.cpp reading.o types.h
	$(CXX) $(CXXFLAGS) -Wall -o popsum -lm reading.o popsum.cpp
wxyz_genotyper: wxyz_genotyper.cpp reading.o calc.o types.h
	$(CXX) $(CXXFLAGS) -Wall -o wxyz_genotyper -lm reading.o wxyz_genotyper.cpp

clean:
	rm -f reading.o calc.o popsum sdpop wxyz_genotyper

