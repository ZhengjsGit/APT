PREFIX = ..
CC = gcc
DIR_INCLUDE=../include
DIR_SRC=./
VPATH = ../include:../src:../lib
LIB_DIR = ../lib
STATIC_LIBS = cgapsio
OBJECTS = lgapsI.o lgapsO.o
HEADERS = lgapsO.h lua.h lualib.h lauxlib.h lgapsI.h 
DBG = 
CFLAGS=-O2
LUA_DIR=..

libgapsio.a: $(OBJECTS)
	$(AR) rcs $@ $(OBJECTS) 

lgapsO.o:	lgapsO.c $(HEADERS)
	$(CC) $(CFLAGS) -c lgapsO.c -I$(DIR_INCLUDE)

lgapsI.o:	lgapsI.c $(HEADERS)
	$(CC) $(CFLAGS) -c lgapsI.c -I$(DIR_INCLUDE)


install :
	cp libgapsio.a $(DIR_INCLUDE)/gapsio.h $(PREFIX)

.PHONY	:clean
clean	:
	rm *.o *.a -rf
