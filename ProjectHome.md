This project aims at providing a library with signal processing algorithms used in turbo receivers. The main application is in the simulation of digital communication systems.

The algorithms are implemented as stand alone Soft Input Soft Output (SISO) modules, from which turbo receivers can be constructed. The programming language is C++. These SISO modules are implemented using [IT++ library](http://sourceforge.net/apps/wordpress/itpp/), however an interface is provided in order to call these modules from [MATLAB](http://www.mathworks.com/). The main development platform is Linux ([openSUSE](http://www.opensuse.org/en/)).

For further details see [Detailed description of SISO modules](SISOindex.md).

An offline HTML help can be generated using [Doxygen](http://www.stack.nl/~dimitri/doxygen/). Go to [Downloads](http://code.google.com/p/turbosiso/downloads/list) tab in order to download the stable version. Alternatively, you could browse the sources using the [Source](http://code.google.com/p/turbosiso/source/browse/#svn/trunk) tab.

**This project has been merged into [release 4.2.0](http://sourceforge.net/projects/itpp/files/itpp/4.2.0/) of IT++ library**