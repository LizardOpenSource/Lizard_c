CC = gcc

CFLAGS=-O3 -fomit-frame-pointer -mavx2 -march=native -std=c99

all :
	$(CC) $(CFLAGS) -c Lizard.c 
	$(CC) $(CFLAGS) -o Lizard Lizard.o
	
run : all
	./Lizard

clean :
	rm -f *.o
	rm -f Lizard

new :
	make clean
	make all
	./Lizard