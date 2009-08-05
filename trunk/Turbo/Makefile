# The program name can be overwritten from command line as follows:
#'make PROGRAM=pccc' or 'make clean PROGRAM=pccc' 

#program name
PROGRAM		=		pccc sccc equalizer IDMA BICM STBICM EXIT_pccc TC_descrambler TC_mud

#SISO class path
VPATH			=		../SISO

#dependencies
DEPS			=		SISO

#compiler options
CXX 			=		gcc
CXXFLAGS 		= 		`itpp-config --cflags-opt` -I$(VPATH)
LDFLAGS			=		`itpp-config --libs-opt`
#DBGFLAGS 		= 		`itpp-config --cflags-debug` -I$(VPATH)
DBGFLAGS 		= 		-I/usr/local/include -Wall -ggdb -O0 -pipe -I$(VPATH)
DBGLIBS 		= 		`itpp-config --libs-debug`

#additional file names
PROGRAM.o       =       $(addsuffix .o, ${PROGRAM})
PROGRAM_d.o     =       $(addsuffix _d.o, ${PROGRAM})
PROGRAM_d       =       $(addsuffix _d, ${PROGRAM})
DEPS.o			= 	 	$(addsuffix .o, $(DEPS))
DEPS_d.o		=		$(addsuffix _d.o, $(DEPS))

#target definition
.PHONY : all debug clean

%.o : %.cpp
	$(CXX) -c $(CXXFLAGS) $< -o $@

$(PROGRAM) : $(PROGRAM.o) $(DEPS.o)
	$(CXX) $(PROGRAM.o) $(DEPS.o) $(LDFLAGS) -o $(PROGRAM)

all : $(PROGRAM)

%_d.o : %.cpp
	$(CXX) -c $(DBGFLAGS) $< -o $@

$(PROGRAM_d) : $(PROGRAM_d.o) $(DEPS_d.o)
	$(CXX) $(PROGRAM_d.o) $(DEPS_d.o) $(DBGLIBS) -o $(PROGRAM_d)

debug : $(PROGRAM_d)

clean :
	$(RM) $(PROGRAM) $(PROGRAM_d) $(PROGRAM_d.o) $(PROGRAM.o) $(DEPS.o) $(DEPS_d.o)
