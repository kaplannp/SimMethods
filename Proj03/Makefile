CC = g++
CFLAGS = -ggdb -fmax-errors=1

sanitized: sim.o main.o
	$(CC) $^ -o $@ $(CFLAGS) -fsanitize=address
main: sim.o main.o 
	$(CC) $^ -o $@ $(CFLAGS)

sim.o: sim.cpp sim.h
	$(CC) -c sim.cpp $(CFLAGS)

main.o: main.cpp
	$(CC) -c main.cpp $(CFLAGS)
clean:
	rm -r *.o
