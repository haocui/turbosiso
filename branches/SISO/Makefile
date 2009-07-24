CC = gcc
FLAGS = `itpp-config --cflags`
DBGFLAGS = `itpp-config --cflags-debug`
#FLAGS = `itpp-config --cflags-opt`
LIBS = `itpp-config --libs`
DBGLIBS = `itpp-config --libs-debug`
#LIBS = `itpp-config --libs-opt`

DBGOPTS = -ggdb
PROFOPTS = -pg

all : IDMA IDMA_eq pg progress_timer

debug : IDMA_eq_d progress_timer_d

IDMA : IDMA.cpp Progress.cpp SISO.cpp
	$(CC) $(FLAGS) -o IDMA IDMA.cpp $(LIBS)

IDMA_eq : IDMA_eq.cpp Progress.cpp SISO.cpp
	$(CC) $(FLAGS) -o IDMA_eq IDMA_eq.cpp $(LIBS)

IDMA_eq_d : IDMA_eq.cpp Progress.cpp SISO.cpp
	$(CC) $(DBGOPTS) $(DBGFLAGS) -o IDMA_eq_d IDMA_eq.cpp $(DBGLIBS)

pg : pg.cpp Progress.cpp
	g++ -o pg pg.cpp
	
pg_d : pg.cpp Progress.cpp
	g++ $(DBGOPTS) -o pg_d pg.cpp

progress_timer : progress_timer.cpp
	$(CC) $(FLAGS) -o progress_timer progress_timer.cpp $(LIBS)
	
progress_timer_d : progress_timer.cpp
	$(CC) $(DBGFLAGS) -o progress_timer_d progress_timer.cpp $(DBGLIBS)

help :
	doxygen Doxyfile

clean :
	rm -f IDMA_eq IDMA_eq_d
	rm -f progress_timer progress_timer_d
	rm -f pg pg_d
	rm -f *~
