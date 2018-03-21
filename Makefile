CXX           = g++
CXXFLAGS      = -std=c++11 $(DEFINES) 
SOURCE=./src/

all: $(SOURCE)main.o
	g++ -Ofast -march=native -std=c++11  -o giscup $(SOURCE)main.o -Wall

main.o: $(SOURCE)main.cpp $(SOURCE)freespace.h \
		$(SOURCE)trajectory.h \
		$(SOURCE)point.h \
		$(SOURCE)cell.h \
		$(SOURCE)input.h \
		$(SOURCE)RTree.h \
		$(SOURCE)query.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o main.o main.cpp
clean:
	$(RM) *.o $(objects) giscup
	(cd src ; $(RM) *.o)
