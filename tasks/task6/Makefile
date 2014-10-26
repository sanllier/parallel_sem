all: task6
OBJECTS = parparser.o main.o

task6: $(OBJECTS)
	mpicxx $^ -o $@

main.o: main.cpp
	mpicxx -c -O3 -I../../lib $^ -o $@

parparser.o: ../../lib/parparser/parparser.cpp
	mpicxx -c -O3 $^ -o $@
	
clean:
	rm -rf *.o task6
