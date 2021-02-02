CC     = g++
CFLAGS = -c -Wall -O3 -std=c++0x -DDEBUG -g
LFLAGS = -g -Wall -O3
EXECUTABLE = SimulateLoci

all: SimulateLoci

OBJECTS = LociSimulator.o SuperPopulation.o Population.o Bacterium.o LocusSequencesSimulator.o LocusTreeStructure.o Auxilaries.o

SimulateLoci : $(OBJECTS)
	$(CC) $(LFLAGS) -o $@ $(OBJECTS)

%.o : %.cpp SuperPopulation.h Population.h Bacterium.h LocusSequencesSimulator.h LocusTreeStructure.h Auxilaries.h
	$(CC) $(CFLAGS) $<

clean :
	rm -f *.o $(EXECUTABLE)

