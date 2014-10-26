all: rand_matr
OBJECTS = parparser.o main.o

rand_matr: $(OBJECTS)
	mpicxx $^ -o $@

main.o: main.cpp
	mpicxx -c -I../../lib $^ -o $@

parparser.o: ../../lib/parparser/parparser.cpp
	mpicxx -c $^ -o $@
	
clean:
	rm -rf *.o rand_matr
