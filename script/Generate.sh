#!/bin/bash
CONFIG_DIR="../Config/pkg/lib"
DOC_DIR="../doc"
HEADER_DIR="../include/GeneratedHeaders"
C_DIR="../src"


#Modify codes below unless you know what to do
AUX_DIR="./auxscripts"

test -e ${AUX_DIR}/GenerateHeaders.sh && rm ${AUX_DIR}/GenerateHeaders.sh
#Murging
cat ${AUX_DIR}/funADD.sh >> ${AUX_DIR}/GenerateHeaders.sh
cat ADD >> ${AUX_DIR}/GenerateHeaders.sh
echo "CONFIG_DIR=\"${CONFIG_DIR}\"" >>${AUX_DIR}/GenerateHeaders.sh
echo "DOC_DIR=\"${DOC_DIR}\"" >>${AUX_DIR}/GenerateHeaders.sh
echo "HEADER_DIR=\"${HEADER_DIR}\"" >>${AUX_DIR}/GenerateHeaders.sh
echo "C_DIR=\"${C_DIR}\"" >>${AUX_DIR}/GenerateHeaders.sh
cat ${AUX_DIR}/lefts.sh >> ${AUX_DIR}/GenerateHeaders.sh

#Running
chmod +x ${AUX_DIR}/GenerateHeaders.sh

${AUX_DIR}/GenerateHeaders.sh


