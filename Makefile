CFLAGS = -std=c99 -Wall -O3
LIBFLAGS = -lm
CC = gcc
OUTFOLDER=Debug
HEADERS = msm4g_bases.h     \
          msm4g_lib.h       \
          msm4g_constants.h \
          msm4g_tests.h     \
          msm4g_types.h

.PHONY: directories

OBJECTS = $(OUTFOLDER)/msm4g_bases.o \
          $(OUTFOLDER)/msm4g_lib.o   \
          
LIBRARY = $(OUTFOLDER)/libmsm4g.a          

EXECUTABLE = $(OUTFOLDER)/msm4g
EXECUTABLECODE = msm4g_main.c
TESTEXECUTABLE = $(OUTFOLDER)/msm4g_tests
TESTEXECUTABLECODE = msm4g_tests.c

$(OUTFOLDER)/%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c -o $@ $<

all: directories $(OBJECTS) $(EXECUTABLE) $(TESTEXECUTABLE) $(LIBRARY)

directories: $(OUTFOLDER)

$(OUTFOLDER):
	mkdir -p $(OUTFOLDER)

$(EXECUTABLE): $(OBJECTS) $(EXECUTABLECODE)
	$(CC) $(CFLAGS) -o $(EXECUTABLE) $(OBJECTS) $(EXECUTABLECODE) $(LIBFLAGS)

$(TESTEXECUTABLE): $(OBJECTS) $(TESTEXECUTABLECODE)
	$(CC) $(CFLAGS) -o $(TESTEXECUTABLE) $(OBJECTS) $(TESTEXECUTABLECODE) $(LIBFLAGS)

$(LIBRARY): $(OBJECTS)
	ar r $@ $^
	
test: $(TESTEXECUTABLE)
	@echo "Running tests"
	$(TESTEXECUTABLE)

demo: $(EXECUTABLE)
	$(EXECUTABLE) data/NaClN8.ini 2 2 2 6 2 4 0 0 0 0
	@echo "--------------------------------------------------"
	@echo "Demo run for NaCl N=8 case ended"
	@echo "--------------------------------------------------"
	@echo "potentialEnergyTotal (see above) should be -6.9903"

help:
	@echo "make all   ---> Create executables and the library"
	@echo "make test  ---> Run unit tests"
	@echo "make demo  ---> Run for NaCl N=8 case"
	@echo "make clean ---> Delete executables and the library"

clean:
	rm -f $(OBJECTS) $(LIBRARY) $(EXECUTABLE) $(TESTEXECUTABLE)

