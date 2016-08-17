CC = g++ -std=c++11
LIBS = -lboost_filesystem -lboost_system
LIBS_PNG = -lpng -lpngwriter -lz -lfreetype


all: 
	$(CC) -I/usr/local/opt/boost/include -L /usr/local/opt/boost/lib $(LIBS) -o divide.o divide.cpp tree.cpp	

init:
	$(CC) -pg -g -ffast-math -o init.o init.cpp `freetype-config --cflags` -I/usr/local/include -L/usr/local/lib -I/usr/local/opt/boost/include -L /usr/local/opt/boost/lib $(LIBS) $(LIBS_PNG) 

init2:
	$(CC) -pg -g -ffast-math -o init2.o init2.cpp `freetype-config --cflags` -I/usr/local/include -L/usr/local/lib -I/usr/local/opt/boost/include -L /usr/local/opt/boost/lib $(LIBS) $(LIBS_PNG) 


