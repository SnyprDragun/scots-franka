#
# compiler
#
#CC        = g++
CC       	  = clang++
CXXFLAGS 		= -Wall -Wextra -std=c++11 -O3 -DNDEBUG -DSCOTS_BDD
#
# scots 
#
SCOTSROOT		= ../..
SCOTSINC		= -I$(SCOTSROOT)/src -I$(SCOTSROOT)/utils
#
# cudd 
#
CUDD		    =  -I/opt/local/include -L/opt/local/lib -lcudd


.PHONY: synthesis simulate

TARGET = synthesis simulate

all: $(TARGET)

synthesis: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) $(CUDD) synthesis.cc -o synthesis

simulate: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) $(CUDD) simulate.cc -o simulate

clean:
	rm synthesis simulate
