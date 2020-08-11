CC     = g++
CFLAGS = -c -Wall -O3 -std=c++0x
LFLAGS = -g -Wall -O3
EXECUTABLE = SAPSimulator

all: SAPSimulator

OBJECTS = SAPSimulator.o SuperPopulation.o Population.o Bacterium.o LocusSequencesSimulator.o LocusTreeStructure.o Auxilaries.o WrightFisherProcess.o

SAPSimulator : $(OBJECTS)
	$(CC) $(LFLAGS) -o $@ $(OBJECTS)

%.o : %.cpp SuperPopulation.h Population.h Bacterium.h LocusSequencesSimulator.h LocusTreeStructure.h Auxilaries.h WrightFisherProcess.h
	$(CC) $(CFLAGS) $<

clean :
	rm -f *.o $(EXECUTABLE)

