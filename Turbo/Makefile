# Use 'make program_name' for single program compilation 

#test program name
TESTS			=		pccc sccc equalizer IDMA BICM STBICM EXIT_pccc\
						sova_example pccc_punctured pccc_hybrid\
						TC_RSC

#SISO class path
VPATH			=		/home/bogdan/C++/SISO

#compiler options
CXX 			=		gcc
CXXFLAGS 		= 		`itpp-config --cflags-opt` -I$(VPATH)
LDFLAGS			=		`itpp-config --libs-opt` -L$(VPATH) -lsiso
CXXFLAGS_DEBUG 	= 		-I/usr/local/include -Wall -ggdb -O0 -pipe -I$(VPATH)
LDFLAGS_DEBUG 	= 		`itpp-config --libs-debug` -L$(VPATH) -lsiso_debug

#additional file names
TESTS.o			= 	 	$(addsuffix .o, $(TESTS))
TESTS_d			= 	 	$(addsuffix _d, $(TESTS))

#target definition
.PHONY : all clean debug

#multiple program compilation (implicit rules are used)
all : $(TESTS)

#debug program compilation
debug : $(TESTS_d)

%_d : %.cpp
	$(CXX) $(CXXFLAGS_DEBUG) $< $(LDFLAGS_DEBUG) -o $@

#delete all
clean :
	$(RM) *.o $(TESTS) $(TESTS_d)
