all: extraction-test

INCLUDES=-I../libgeom -I/usr/include/eigen3
CXXFLAGS=-std=c++17 -Wall -pedantic $(INCLUDES)
LIBS=-L../libgeom/release -lgeom

extraction-test: extraction-test.o bezier-extractions.o
	$(CXX) -o $@ $^ $(LIBS)
