#!/usr/bin/env python3
# coding=utf-8
import sys

def delargs(lines):
	pos = 1
	for mem in lines:
		if mem =="--add by python\n":
			print(pos)
			break
		pos =pos +1
	
	em = lines[0:pos]
	return em


arg1 = "Pusher_Type"      + "="   + '"' + str(sys.argv[1]) + '"' + ';' 
arg2 = "EMField_Type"     + "="   + '"' + str(sys.argv[2]) + '"' + ';'
arg3 = "dT"               + "="         + str(sys.argv[3])       + ';'
arg4 = "num_steps"        + "="         + str(sys.argv[4])       + ';'
arg5 = "SavePerNSteps"    + "="         + str(sys.argv[5])       + ';'

conf = open('./template.lua','r')
lines = conf.readlines();
conf.close
empty = delargs(lines);
empty.append(arg1+'\n')
empty.append(arg2+'\n')
empty.append(arg3+'\n')
empty.append(arg4+'\n')
empty.append(arg5+'\n')

conf = open('./template.lua','w+')
conf.writelines(empty)
conf.close()
