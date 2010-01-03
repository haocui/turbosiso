PROGRAMS					=		C_SISOmud C_SISOdescrambler C_SISOnsc C_SISOrsc\
									C_SISOequalizer C_SISOdemapper C_SISOsova

#mex extension on this machine
MEX_EXT						:=		$(shell mexext)

#SISO library path
VPATH						=		../SISO

#compiler flags
CXX							=		mex
CXXFLAGS					=		-I$(VPATH)
#LDFLAGS						=		-litpp -L/opt/acml4.3.0/gfortran64_mp/lib -lacml_mp -L$(VPATH) -lsiso
LDFLAGS                    =       -litpp -L$(VPATH) -lsiso
#LDFLAGS_DEBUG				=		-litpp -L/opt/acml4.3.0/gfortran64_mp/lib -lacml_mp -L$(VPATH) -lsiso_debug
LDFLAGS_DEBUG               =       -litpp -L$(VPATH) -lsiso_debug

#additional file names
PROGRAMS.MEX_EXT			=		$(addsuffix .${MEX_EXT}, ${PROGRAMS})
PROGRAMS_D.MEX_EXT			=		$(addsuffix _d.${MEX_EXT}, ${PROGRAMS})

#target definitions
.PHONY : all debug clean

all : $(PROGRAMS.MEX_EXT)

%.${MEX_EXT} : %.cpp
	$(CXX) $(CXXFLAGS) $< $(LDFLAGS) -output $@

debug : $(PROGRAMS_D.MEX_EXT)

%_d.${MEX_EXT} : %.cpp
	$(CXX) -g $(CXXFLAGS) $< $(LDFLAGS_DEBUG) -output $@

clean :
	$(RM) $(PROGRAMS.MEX_EXT) $(PROGRAMS_D.MEX_EXT)
