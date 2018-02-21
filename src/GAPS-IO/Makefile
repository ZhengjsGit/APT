PREFIX = ./test/GAPS-IO

INCLUDE_DIR=./include
LUALIB_DIR=./deps/liblua-5.2.3
CC = gcc
CFLAGS = -O2 -I$(INCLUDE_DIR) $(DBG)
AR = ar
ARFLAGS = -rcsv
MAKE=make
VPATH = include:src:lib:deps/liblua-5.2.3

OBJECTS = lgapsI.o lgapsO.o 
HEADERS = lgapsO.h lgapsI.h
DBG =
ALL = lgapsIO liblua

TO_INC= lua.h luaconf.h lualib.h lauxlib.h lua.hpp
TO_LIB= liblua.a

MK_DIR=$(shell pwd)

all : $(ALL)

lgapsIO:
	cd src && $(MAKE) ;\
	$(AR) $(ARFLAGS) libgapsio.a $(OBJECTS) ;\
	cp libgapsio.a $(MK_DIR)

liblua:
	cd $(LUALIB_DIR) && $(MAKE)
	cd $(LUALIB_DIR) && cp $(TO_LIB) $(MK_DIR)

install:
	-mkdir -p $(PREFIX)/include
	-mkdir -p $(PREFIX)/lib
	-mkdir -p $(PREFIX)/tools
	cp libgapsio.a liblua.a $(PREFIX)/lib
	cp  $(INCLUDE_DIR)/gapsio.h $(INCLUDE_DIR)/lua.h  $(INCLUDE_DIR)/luaconf.h $(INCLUDE_DIR)/lualib.h $(INCLUDE_DIR)/lauxlib.h $(INCLUDE_DIR)/lua.hpp $(PREFIX)/include
	cp tools/* $(PREFIX)/tools -r

clean:
	-rm *.a
	cd $(LUALIB_DIR) && make clean
	cd src && make clean

.PHONY : all lgapsIO liblua clean
