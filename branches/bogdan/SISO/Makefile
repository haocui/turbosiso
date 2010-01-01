#Makefile for creating the SISO shared library and its documentation

COMPONENTS		=		SISO_rsc SISO_nsc SISO_eq SISO_mud SISO_dem EXIT STC
LIBRARY_NAME	=		libsiso
MAJOR_NO		=		2
MINOR_NO		=		0
RELEASE_NO		=		0
LIBRARY_VERSION =		$(MAJOR_NO).$(MINOR_NO).$(RELEASE_NO)
INSTALL_PATH	=		/usr/local/lib

#additional file names
COMPONENTS.o       =       $(addsuffix .o, ${COMPONENTS})
COMPONENTS.cpp     =       $(addsuffix .cpp, ${COMPONENTS})
COMPONENTS_d.o     =       $(addsuffix _d.o, ${COMPONENTS})

#compiler options
CXX 			=		gcc
CXXFLAGS 		= 		`itpp-config --cflags-opt`
LDFLAGS			=		`itpp-config --libs-opt`
CXXFLAGS_DEBUG	= 		-I/usr/local/include -Wall -ggdb -O0 -pipe
LDFLAGS_DEBUG	= 		`itpp-config --libs-debug`

#target definition
.PHONY : all lib lib_debug help wiki install clean

all : lib lib_debug

lib : $(COMPONENTS.o)
	$(CXX) -shared -Wl,-soname,$(LIBRARY_NAME).so.$(MAJOR_NO) \
	-o $(LIBRARY_NAME).so.$(LIBRARY_VERSION) $(COMPONENTS.o) $(LDFLAGS)

$(COMPONENTS.o) : $(COMPONENTS.cpp)

%.o : %.cpp
	$(CXX) -fPIC -c $(CXXFLAGS) $< -o $@	

lib_debug : $(COMPONENTS_d.o)
	$(CXX) -shared -Wl,-soname,$(LIBRARY_NAME)_debug.so.$(MAJOR_NO) \
	-o $(LIBRARY_NAME)_debug.so.$(LIBRARY_VERSION) $(COMPONENTS_d.o) $(LDFLAGS_DEBUG)

$(COMPONENTS_d.o) : $(COMPONENTS.cpp)

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

install :
	cp -f $(LIBRARY_NAME).so.$(LIBRARY_VERSION) $(INSTALL_PATH)
	cp -f $(LIBRARY_NAME)_debug.so.$(LIBRARY_VERSION) $(INSTALL_PATH)
	cd $(INSTALL_PATH)
	ln -sf $(LIBRARY_NAME).so.$(LIBRARY_VERSION) $(LIBRARY_NAME).so.$(MAJOR_NO)
	ln -sf $(LIBRARY_NAME).so.$(MAJOR_NO) $(LIBRARY_NAME).so
	ln -sf $(LIBRARY_NAME)_debug.so.$(LIBRARY_VERSION) $(LIBRARY_NAME)_debug.so.$(MAJOR_NO)
	ln -sf $(LIBRARY_NAME)_debug.so.$(MAJOR_NO) $(LIBRARY_NAME)_debug.so
	/sbin/ldconfig $(INSTALL_PATH)

clean :
	$(RM) -r html wiki $(COMPONENTS.o) $(COMPONENTS_d.o)\
			$(LIBRARY_NAME).* $(LIBRARY_NAME)_debug.*
