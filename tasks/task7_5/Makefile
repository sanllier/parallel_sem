all: task7_5
OBJECTS = main.o parparser.o

task7_5: $(OBJECTS)
	g++ -std=c++11 -L../../lib/parparser $^ -o $@

main.o: main.cpp
	g++ -std=c++11 -I../../lib/parparser -I../../lib/matrix -o $@

parparser.o: parparser.cpp
	g++ -std=c++11 $^ -o $@
	
clean:
	rm -rf *.o task7_5