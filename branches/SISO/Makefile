#Makefile for creating the SISO shared library and its documentation

COMPONENTS		=		SISO_rsc SISO_nsc SISO_eq SISO_mud SISO_dem EXIT STC
LIBRARY_NAME	=		libsiso
MAJOR_NO		=		1
MINOR_NO		=		1
RELEASE_NO		=		1
LIBRARY_VERSION =		$(MAJOR_NO).$(MINOR_NO).$(RELEASE_NO)

#additional file names
COMPONENTS.o       =       $(addsuffix .o, ${COMPONENTS})
COMPONENTS_d.o     =       $(addsuffix _d.o, ${COMPONENTS})
COMPONENTS_s.o     =       $(addsuffix _s.o, ${COMPONENTS})

#compiler options
CXX 			=		gcc
CXXFLAGS 		= 		`itpp-config --cflags-opt`
LDFLAGS			=		`itpp-config --libs-opt`
CXXFLAGS_DEBUG	= 		-I/usr/local/include -Wall -ggdb -O0 -pipe
LDFLAGS_DEBUG	= 		`itpp-config --libs-debug`

#target definition
.PHONY : all lib lib_debug static help wiki

all : lib lib_debug

lib : $(COMPONENTS.o)
	$(CXX) -shared -Wl,-soname,$(LIBRARY_NAME).so.$(MAJOR_NO) \
	-o $(LIBRARY_NAME).so.$(LIBRARY_VERSION) $(COMPONENTS.o) $(LDFLAGS)
	ln -sf $(LIBRARY_NAME).so.$(LIBRARY_VERSION) $(LIBRARY_NAME).so.$(MAJOR_NO)
	ln -sf $(LIBRARY_NAME).so.$(MAJOR_NO) $(LIBRARY_NAME).so

%.o : %.cpp
	$(CXX) -fPIC -c $(CXXFLAGS) $< -o $@

static : $(COMPONENTS_s.o)
	ar -cr $(LIBRARY_NAME).a $(COMPONENTS_s.o)
	
%_s.o : %.cpp
	$(CXX) -c $(CXXFLAGS) $< -o $@	

lib_debug : $(COMPONENTS_d.o)
	$(CXX) -shared -Wl,-soname,$(LIBRARY_NAME)_debug.so.$(MAJOR_NO) \
	-o $(LIBRARY_NAME)_debug.so.$(LIBRARY_VERSION) $(COMPONENTS_d.o) $(LDFLAGS_DEBUG)
	ln -sf $(LIBRARY_NAME)_debug.so.$(LIBRARY_VERSION) $(LIBRARY_NAME)_debug.so.$(MAJOR_NO)
	ln -sf $(LIBRARY_NAME)_debug.so.$(MAJOR_NO) $(LIBRARY_NAME)_debug.so

%_d.o : %.cpp
	$(CXX) -fPIC -c $(CXXFLAGS_DEBUG) $< -o $@

help :
	doxygen Doxyfile
	
wiki : help
	mkdir -p wiki
	cd wiki && html2wiki \
	 --dialect GoogleCode \
	 --base-uri=http://turbosiso.googlecode.com/svn/wiki/\
	 ../html/index.html > SISOindex.wiki
	 cd wiki && html2wiki \
	 --dialect GoogleCode \
	 --base-uri=http://turbosiso.googlecode.com/svn/wiki/\
	 ../html/pccc.html > pccc.wiki
	 cd wiki && html2wiki \
	 --dialect GoogleCode \
	 --base-uri=http://turbosiso.googlecode.com/svn/wiki/\
	 ../html/sccc.html > sccc.wiki
	 cd wiki && html2wiki \
	 --dialect GoogleCode \
	 --base-uri=http://turbosiso.googlecode.com/svn/wiki/\
	 ../html/idma.html > idma.wiki
	 cd wiki && html2wiki \
	 --dialect GoogleCode \
	 --base-uri=http://turbosiso.googlecode.com/svn/wiki/\
	 ../html/equalizer.html > equalizer.wiki
	 cd wiki && html2wiki \
	 --dialect GoogleCode \
	 --base-uri=http://turbosiso.googlecode.com/svn/wiki/\
	 ../html/bicm.html > bicm.wiki
	 cd wiki && html2wiki \
	 --dialect GoogleCode \
	 --base-uri=http://turbosiso.googlecode.com/svn/wiki/\
	 ../html/stbicm.html > stbicm.wiki

clean :
	$(RM) -r html wiki $(COMPONENTS.o) $(COMPONENTS_d.o)\
			$(LIBRARY_NAME).* $(LIBRARY_NAME)_debug.*\
			$(LIBRARY_NAME).a $(COMPONENTS_s.o)
