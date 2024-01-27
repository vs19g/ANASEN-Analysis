########################################################################
#
#
#########################################################################

CC	=	g++

#COPTS	=	-fPIC -DLINUX -O2 -std=c++17 -lpthread
COPTS	=	-fPIC -DLINUX -g -O0 -Wall -std=c++17 -lpthread

ROOTLIBS = `root-config --cflags --glibs`

ALL = Mapper 

#########################################################################

all	:	$(ALL)

clean	:
		/bin/rm -f $(OBJS) $(ALL)

Mapper	:	Mapper.cpp  mapping.h
		@echo "--------- making Mapper"
		$(CC) $(COPTS) -o Mapper Mapper.cpp $(ROOTLIBS)
