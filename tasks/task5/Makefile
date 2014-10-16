all: task5
OBJECTS = parparser.o main.o

task5: $(OBJECTS)
	mpicxx -std=c++11 -L../../lib $^ -o $@

main.o: main.cpp
	mpicxx -std=c++11 -I../../lib -o $@

../../lib/parparser/parparser.o: ../../lib/parparser/parparser.cpp
	mpicxx -std=c++11 $^ -o $@
	
clean:
	rm -rf *.o task5