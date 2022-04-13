CC = g++
CFLAGS = -g -fmax-errors=4

main: main.o
	$(CC) $< -o $@ 

%.o: %.cpp
	$(CC) -c -o $@ $< $(CFLAGS)

scratch: scratch.cpp
	$(CC) scratch.cpp -o scratch && ./scratch
