all: task7_6
OBJECTS = parparser.o main.o

task7_6: $(OBJECTS)
	mpicxx $^ -o $@

main.o: main.cpp
	mpicxx -c -O3 -I../../lib $^ -o $@

parparser.o: ../../lib/parparser/parparser.cpp
	mpicxx -c -O3 $^ -o $@
	
clean:
	rm -rf *.o task7_6
