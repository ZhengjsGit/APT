import struct
#import numpy as np

GAPS_IO_GlobalTypeLen=[4,8,2,4,8,2,4,8,16,16,16,1,1]
def get_type_from_int(inp):
	if inp==0:
		return numpy.int32
	elif inp==1:
		return numpy.int64
	elif inp==2:
		return numpy.int16
	elif inp==3:
		return numpy.uint32
	elif inp==4:
		return numpy.uint64
	elif inp==5:
		return numpy.uint16
	elif inp==6:
		return numpy.float32
	elif inp==7:
		return numpy.float64
	elif inp==8:
		return numpy.float128
	elif inp==9:
		return numpy.float128
	elif inp==10:
		return numpy.float128
	elif inp==11:
		return numpy.int8
	elif inp==12:
		return numpy.uint8

def readi64(fp):
	l1=struct.unpack('i', fp.read(4))[0]
	l2=struct.unpack('i', fp.read(4))[0]
	return l1+l2*(1<<32)

def GAPS_IO_Load(filename):
	fp=open(filename,"rb")
	version=readi64(fp)
	tp=readi64(fp)
	numdim=readi64(fp)
	alldimver=range(numdim+1)
	numperstep=1
	for i in xrange(numdim):
		alldimver[i]=readi64(fp)
		numperstep*=alldimver[i]
	datahead=8*(3+numdim)
	#print(version,tp,numdim,alldimver)
	fp.seek(0,2)
	allfilelen=fp.tell()
	numtimestep=(allfilelen-datahead)/numperstep/GAPS_IO_GlobalTypeLen[tp]
	alldimver[numdim]=numtimestep
	alldimver.reverse()
	alldatalen=numperstep*numtimestep
	#print(alldatalen)
	fp.seek(datahead)
	#ret=zeros(alldimver)
	ret=numpy.fromfile(fp,dtype=get_type_from_int(tp),count=alldatalen)
	fp.close()
	return reshape(ret,alldimver)
