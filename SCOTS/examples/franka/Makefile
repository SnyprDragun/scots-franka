#
# compiler
#
#CC        = g++
CC       	  = clang++
CXXFLAGS 		= -g -Wall -Wextra -std=c++11
CXXFLAGS 		= -Wall -Wextra -std=c++11 -O3 -DNDEBUG
CXXFLAGS += -I/usr/include/eigen3

 
#
# scots 
#
SCOTSROOT		= ../..
SCOTSINC		= -I$(SCOTSROOT)/src  -I$(SCOTSROOT)/utils/
 
.PHONY: franka
 
TARGET = franka
 
all: $(TARGET)
 
 
franka:
	$(CC) $(CXXFLAGS) $(SCOTSINC) franka.cc -o  franka
 
 
clean:
	rm  franka
