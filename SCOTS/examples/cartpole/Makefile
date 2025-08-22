#
# compiler
#
#CC        = g++
CC       	  = clang++
CXXFLAGS 		= -g -Wall -Wextra -std=c++11
CXXFLAGS 		= -Wall -Wextra -std=c++11 -O3 -DNDEBUG

#
# scots 
#
SCOTSROOT		= ../..
SCOTSINC		= -I$(SCOTSROOT)/src  -I$(SCOTSROOT)/utils/

.PHONY: thermal simulate

TARGET = thermal simulate

all: $(TARGET)

thermal: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) thermal.cc -o thermal

simulate: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) thermal.cc -o simulate

clean:
	rm thermal simulate
