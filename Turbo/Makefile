# The program name can be overwritten from command line as follows:
#'make PROGRAM=pccc' or 'make clean PROGRAM=pccc' 

#program name
PROGRAM			=		pccc sccc equalizer IDMA BICM STBICM EXIT_pccc TC_descrambler TC_mud

#mex program name
MEX_PROGRAM		=		C_SISOrsc C_SISOnsc C_SISOequalizer C_SISOmud C_SISOdescrambler C_SISOdemapper

#SISO class path
VPATH			=		../SISO

#compiler options
CXX 			=		gcc
CXXFLAGS 		= 		`itpp-config --cflags-opt` -I$(VPATH)
LDFLAGS			=		`itpp-config --libs-opt`
#DBGFLAGS 		= 		`itpp-config --cflags-debug` -I$(VPATH)
DBGFLAGS 		= 		-I/usr/local/include -Wall -ggdb -O0 -pipe -I$(VPATH)
DBGLIBS 		= 		`itpp-config --libs-debug`
MEX_OPTS		=		-litpp -L/opt/acml4.2.0/gfortran64_mp/lib -lacml_mp -I$(VPATH)

#additional file names
PROGRAM_d       =       $(addsuffix _d, ${PROGRAM})
MEX_PROGRAM.mexa64		=		$(addsuffix .mexa64, ${MEX_PROGRAM})
MEX_PROGRAM_d.mexa64	=		$(addsuffix _d.mexa64, ${MEX_PROGRAM})

#target definition
.PHONY : all debug mex mex_debug clean

% : %.cpp
	$(CXX) $(CXXFLAGS) $< $(LDFLAGS) -o $@

all : $(PROGRAM)

%_d : %.cpp
	$(CXX) $(DBGFLAGS) $< $(DBGLIBS) -o $@

debug : $(PROGRAM_d)

%.mexa64 : %.cpp
	mex $< $(MEX_OPTS) -output $@ 

mex : $(MEX_PROGRAM.mexa64)

%_d.mexa64 : %.cpp
	mex -g $< $(MEX_OPTS) -output $@ 

mex_debug : $(MEX_PROGRAM_d.mexa64)

clean :
	$(RM) $(PROGRAM) $(PROGRAM_d) $(MEX_PROGRAM.mexa64) $(MEX_PROGRAM_d.mexa64)
