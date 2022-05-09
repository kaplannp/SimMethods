CC = g++
CFLAGS = -ggdb -fmax-errors=1

main: sim.o main.o 
	$(CC) $^ -o $@ $(CFLAGS)

clean:
	rm -r *.o
