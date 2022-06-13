CC = gcc
LD = gcc
CFLAGS = -O3 -Wall
LDFLAGS = -lm
RM = /bin/rm -f
OBJS = main.o
EXECUTABLE = strassen

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJS)
	$(LD) -o $(EXECUTABLE) $(OBJS) $(LDFLAGS)

main.o: main.c
	$(CC) $(CFLAGS) -c main.c

clean:
	$(RM) $(EXECUTABLE) $(OBJS)