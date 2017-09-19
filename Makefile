CC = gcc

CFLAGS=-O3 -fomit-frame-pointer -msse2avx -mavx2 -march=native -std=c99

all :
	$(CC) $(CFLAGS) -c Lizard.c fips202.c 
	$(CC) $(CFLAGS) -o Lizard Lizard.o fips202.o 
	
run : all
	./Lizard

clean :
	rm -f *.o
	rm -f Lizard

new :
	make clean
	make all
	./Lizard