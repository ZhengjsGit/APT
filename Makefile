# You can set compiling parameters here

# Set the installation directory
PREFIX = ./DataAnalysis
# Set the directory hdf5 library
HDF5_ROOT = /install/hdf5_1.8.19_openmpi_1.8.1
# Set the Format of output data:(Only HDF5 and GAPSIO are available)
#DATA_FORMAT = HDF5
DATA_FORMAT = GAPSIO

# Please do NOT change anything below!!!
APT_ROOT = $(shell pwd)
export APT_ROOT
export HDF5_ROOT
export DATA_FORMAT	 

.PHONY : All
All : APT

APT:
	make -C src

install :
	-mkdir -p $(PREFIX)
	-cp src/APT.out $(PREFIX)
	-cp Config/* $(PREFIX) -r
	-cp doc $(PREFIX) -r
	-cp tools $(PREFIX) -r

.PHONY	:clean
clean	:
	cd src;make clean
	cd src/GAPS-IO;make clean
	-rm lib/*
