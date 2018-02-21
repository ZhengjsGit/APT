VPATH = ../include:../lib
CC = mpicc 
AR = ar

SRC_DIR = $(shell pwd)
LIB_DIR = $(APT_ROOT)/lib
INCLUDE_DIR=$(APT_ROOT)/include
MPI_DIR = /usr/local
DBG = -g
#DBG =
ifeq ($(DATA_FORMAT),HDF5)
DEFINE = -DDATA_FORMAT_HDF5
else
DEFINE = -DDATA_FORMAT_GAPSIO
endif
CFLAGS = $(DBG) $(DEFINE) -std=gnu89 -I$(INCLUDE_DIR) -I$(INCLUDE_DIR)/GeneratedHeaders -I$(HDF5_ROOT)/include -O3 -ffast-math 
HEADERS = $(wildcard ${INCLUDE_DIR}/*.h)
export CC AR SRC_DIR LIB_DIR INCLUDE_DIR DBG DEFINE CFLAGS HEADER


OBJ_SRC = main.o APT_Output.o Initialization.o ParticlePusher.o External_Forces.o RT_ARG.o \
          ParticleStruct.o EM_Field.o APT_PubFunc.o WYLfunc.o output_hdf5.o IO_Tools.o
OBJ_EMField = $(patsubst %.c,%.o,$(wildcard EMField/*.c))
OBJ_ExtForce= $(patsubst %.c,%.o,$(wildcard ExtForce/*.c))
OBJ_Pusher= $(patsubst %.c,%.o,$(wildcard Pusher/*.c))
OBJ_Init_Aclr= $(patsubst %.c,%.o,$(wildcard Initialization/Acceleration/*.c))
OBJ_Init_P= $(patsubst %.c,%.o,$(wildcard Initialization/Momentum/*.c))
OBJ_Init_X= $(patsubst %.c,%.o,$(wildcard Initialization/Position/*.c))

OBJ = $(OBJ_SRC) $(OBJ_EMField) $(OBJ_ExtForce) $(OBJ_Pusher) $(OBJ_Init_X) \
      $(OBJ_Init_P) $(OBJ_Init_Aclr)  
OBJ_IO = $(APT_ROOT)/lib/libgapsio.a $(APT_ROOT)/lib/liblua.a	  


ALL = APT.out Src EMField ExtForce Pusher Initialization 

.PHONY : All
All : $(ALL)

APT.out: $(OBJ) $(HEADERS) libgapsio.a
ifeq ($(DATA_FORMAT),HDF5)
	$(CC) -o APT.out $(OBJ) $(LIB_DIR)/libgapsio.a $(LIB_DIR)/liblua.a -I$(INCLUDE_DIR) -L$(HDF5_ROOT)/lib -lhdf5 -lm -ldl 
else
	$(CC) -o APT.out $(OBJ) $(LIB_DIR)/libgapsio.a $(LIB_DIR)/liblua.a -I$(INCLUDE_DIR) -lm -ldl 
endif
	ctags -R . ../include

Src : $(OBJ_SRC)

EMField :
	$(MAKE) -C $@ 

ExtForce :
	$(MAKE) -C $@

Pusher:
	$(MAKE) -C $@

Initialization :
	$(MAKE) -C $@

libgapsio.a :
	$(MAKE) -C GAPS-IO
	-mkdir $(LIB_DIR)
	-mv GAPS-IO/*.a $(LIB_DIR)

.PHONY	:clean
clean	:
	-rm $(OBJ) APT.out -rf
