#!/bin/bash

APT_ROOT=$1

Dirh1=include
Dirh2=include/GeneratedHeaders

Dirc1=src
Dirc2=src/EMField
Dirc3=src/ExtForce
Dirc4=src/Initialization/Acceleration
Dirc5=src/Initialization/Momentum
Dirc6=src/Initialization/Position
Dirc7=src/Pusher

Dirsh1=script
Dirsh2=script/auxscripts

FILES="$APT_ROOT/$Dirh1/*.h $APT_ROOT/$Dirh2/*.h $APT_ROOT/$Dirc1/*.c $APT_ROOT/$Dirc2/*.c $APT_ROOT/$Dirc3/*.c $APT_ROOT/$Dirc4/*.c $APT_ROOT/$Dirc5/*.c $APT_ROOT/$Dirc6/*.c $APT_ROOT/$Dirc7/*.c $APT_ROOT/$Dirsh1/*.sh $APT_ROOT/$Dirsh2/*.sh $APT_ROOT/$Dirsh1/ADD $APT_ROOT/tools/mathematica/*.m $APT_ROOT/tools/mathematica/*.nb "

cat $FILES | wc -l
