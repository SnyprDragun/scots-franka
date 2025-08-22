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
 
.PHONY: turtlebot
 
TARGET = turtlebot
 
all: $(TARGET)
 
 
turtlebot:
	$(CC) $(CXXFLAGS) $(SCOTSINC) turtlebot.cc -o  turtlebot
 
 
clean:
	rm  turtlebot
