CC=g++
CFLAGS= -Wall -std=c++11 -pedantic
#CFLAGS= -fpermissive
LDFLAGS=
SOURCES=waveguide.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=waveguide
COMMON=

all: clean waveguide

waveguide:
	$(CC) $(CFLAGS) $(SOURCES) -o waveguide
	
test: clean
	./waveguide 0.01 0.01

clean:
	
.PHONY : all clean
