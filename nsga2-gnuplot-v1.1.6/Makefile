# Makefile for compiling NSGA-II source code
CC=g++
LD=g++
RM=rm -f
LIBS=-lgsl -lgslcblas -lm
LDFLAGS=-Ofast
CFLAGS=-Wall -ansi -pedantic -g
OBJS:=$(patsubst %.cpp,%.o,$(wildcard *.cpp))
OBJS+=$(patsubst %.c,%.o,$(wildcard *.c))
MAIN=nsga2r
all:$(MAIN)
$(MAIN):$(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o $(MAIN) -lm $(LIBS)
%.o: %.c global.h rand.h design_problem.h
	$(CC) $(CFLAGS) -c $<
clean:
	$(RM) $(OBJS)

