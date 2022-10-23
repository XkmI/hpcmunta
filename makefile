PROGRAM = newtonUsed newtonNaive newtonQuicc newtonQuiccUC newtonZuperfast
CC = gcc
CFLAGS = -O2 -Wall -fcx-limited-range#-g
LIBS = -lm#-lpthread -lgsl -lgslcblas

.PHONY: all
all: $(PROGRAM)

$(PROGRAM) : %: %.c;
	$(CC) -o $@ $< $(CFLAGS) $(LIBS)

.PHONY: clean
clean:
	rm -f *.o
	rm -f *.ppm
	rm -f $(PROGRAM)
