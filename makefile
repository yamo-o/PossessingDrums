CC = g++
INC += -I./math -I./Audio -I./Algorithm
CFLAGS = -c -g -Wall -Wextra
OBJS=$(patsubst %.cpp,%.o,$(shell find ./ -name "*.cpp" -print))

all: PossessingDrums

%.o: %.cpp
	$(CC) -O3 -c -o $@ $< $(INC) $(CFLAGS)

PossessingDrums.o: PossessingDrums.cpp
	$(CC) -O3 -c PossessingDrums.cpp $(INC) 

clean:
	rm PossessingDrums
	rm *.o