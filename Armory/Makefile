########################################################################
#
#
#########################################################################

CC	=	g++

#COPTS	=	-fPIC -DLINUX -O2 -std=c++17 -lpthread
COPTS	=	-fPIC -DLINUX -g -O0 -Wall -std=c++17 -lpthread

ROOTLIBS = `root-config --cflags --glibs`

ALL = Mapper EventBuilder#AnasenMS

#########################################################################

all	:	$(ALL)

clean	:
		/bin/rm -f $(OBJS) $(ALL)

Mapper	:	Mapper.cpp  ../mapping.h ClassDet.h
		@echo "--------- making Mapper"
		$(CC) $(COPTS) -o Mapper Mapper.cpp $(ROOTLIBS)

# AnasenMS : constant.h Isotope.h ClassTransfer.h ClassSX3.h ClassPW.h ClassAnasen.h anasenMS.cpp 
# 		@echo "--------- making ANASEN Monte Carlo"
# 		$(CC) $(COPTS) -o AnasenMS anasenMS.cpp $(ROOTLIBS)

EventBuilder : EventBuilder.cpp ClassData.h fsuReader.h Hit.h
		@echo "--------- making EventBuilder"
		$(CC) $(COPTS) -o EventBuilder EventBuilder.cpp $(ROOTLIBS) 
