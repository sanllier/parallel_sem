all: task7_5
OBJECTS = parparser.o main.o

task7_5: $(OBJECTS)
	mpicxx $^ -o $@

main.o: main.cpp
	mpicxx -c -O3 -I../../lib $^ -o $@

parparser.o: ../../lib/parparser/parparser.cpp
	mpicxx -c -O3 $^ -o $@
	
clean:
	rm -rf *.o task7_5
