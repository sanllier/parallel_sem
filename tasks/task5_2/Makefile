all: task5_2
OBJECTS = parparser.o main.o

task5_2: $(OBJECTS)
	mpicxx $^ -o $@

main.o: main.cpp
	mpicxx -c -O3 -I../../lib $^ -o $@

parparser.o: ../../lib/parparser/parparser.cpp
	mpicxx -c -O3 $^ -o $@
	
clean:
	rm -rf *.o task5_2
