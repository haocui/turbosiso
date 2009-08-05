all : help

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
	$(RM) -r html wiki
