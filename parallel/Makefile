CC = gcc
LD = gcc
CFLAGS = -O3 -Wall -fopenmp
LDFLAGS = -lm -fopenmp -g
RM = /bin/rm -f
OBJS = matrix.o main.o
EXECUTABLE = strassen

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJS)
	$(LD) -o $(EXECUTABLE) $(OBJS) $(LDFLAGS)

main.o: main.c
	$(CC) $(CFLAGS) -c main.c

matrix.o: matrix.c matrix.h
	$(CC) $(CFLAGS) -c matrix.c

clean:
	$(RM) $(EXECUTABLE) $(OBJS)