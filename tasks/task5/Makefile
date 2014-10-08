all: task5
OBJECTS = main.o parparser.o
OBJECTSDEBUG = main_debug.o helpers_debug.o kernel_debug.o

task5: $(OBJECTS)
	g++ -std=c++11 -L../../lib/parparser $^ -o $@

main.o: main.cpp
	g++ -std=c++11 -I../../lib/parparser -I../../lib/matrix -o $@

parparser.o: parparser.cpp
	g++ -std=c++11 $^ -o $@
	
clean:
	rm -rf *.o task5