#
EXT=.exe
NAME=DCG
DIR=../bin/
NAMEFUL=$(DIR)/$(NAME)$(EXT)

FC = gfortran
FFLAGS = -c -u -O -fcray-pointer
CC = gcc
CPP = g++
CFLAGS = -c -O
CPPFLAGS = -c -O
LDFLAGS = -O
INCLUDEDIRS =
LIBRARIES = -lpthread -lgfortran -lm -lstdc++
LIBS = -L/usr/local/gfortran/lib/

.c.o :
	$(CC) $(CFLAGS) $(INCLUDEDIRS) $<

.cpp.o :
	$(CPP)  $(CPPFLAGS) $(INCLUDEDIRS) $<

.f90 .f.o :
	$(FC) $(FFLAGS) $<

OBJECTS = \
$(NAME).o \
ProcessInput.o \
FindTemp.o ProcessTemp.o \
RandomWalk.o mt.o \
SpectralClustering.o kmns.o \
EnsembleEigen.o basic_arpack.o basic_lapack.o \
ManipDCG.o

$(NAMEFUL) : $(OBJECTS)
	$(CC) -o $(NAMEFUL) $(LDFLAGS) $(OBJECTS) $(INCLUDEDIRS) $(LIBS) $(LIBRARIES) 

all: $(OBJECTS)
	$(CC) -o $(NAMEFUL) $(LDFLAGS) $(OBJECTS) $(INCLUDEDIRS) $(LIBS) $(LIBRARIES) 

clean:
	touch junk.o; rm -f *.o $(NAMEFUL)

$(OBJECTS) :
