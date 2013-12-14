CC = g++
CPPFLAGS = -std=c++11 -Wall -Werror -pedantic -Wno-sign-compare -g -O3
INCLUDES = -Isrc
SOURCES = $(wildcard src/*.cpp)
HEADERS = $(wildcard src/*.h)

bin/Test : src/Test.cpp bin/MLFMM.o bin/BHNode.o $(HEADERS)
	$(CC) $(INCLUDES) $(CPPFLAGS) -o $@ $< bin/MLFMM.o bin/BHNode.o

bin/MLFMM.o : src/MLFMM.cpp $(HEADERS)
	$(CC) -c $(INCLUDES) $(CPPFLAGS) -o $@ $<

bin/BHNode.o : src/BHNode.cpp $(HEADERS)
	$(CC) -c $(INCLUDES) $(CPPFLAGS) -o $@ $<

documentation : 
	doxygen Doxyfile && cd doc/latex && make && open doc/latex/reman.pdf

clean :
	rm -rf bin/*