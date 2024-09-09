#!/bin/bash
#
# This script will download and install gouppy
# params:
# 1 guppy version

if [ x"$1" == x ]; then
        GUPPY_VER='3.4.5'
else
	GUPPY_VER=$1
fi


wget https://cdn.oxfordnanoportal.com/software/analysis/ont-guppy_${GUPPY_VER}_linux64.tar.gz 
if [ $? -eq 0 ]; then
	echo "INSTALLING GUPPY VERSION ${GUPPY_VER}"
else
    echo "GUPPY VERSION ${GUPPY_VER} is not found" 
    exit
fi

tar -zvxf ont-guppy_${GUPPY_VER}_linux64.tar.gz
mv ont-guppy mop_preprocess/bin/
cd mop_preprocess/bin
ln -sf ont-guppy/bin/guppy_* .
ln -sf ont-guppy/lib/* .
cd ../../
if [ ! -e "mop_preprocess/bin/ont-guppy/lib/libz.so" ] ; then
        unlink mop_preprocess/bin/ont-guppy/lib/libz.so
        cd mop_preprocess/bin/ont-guppy/lib/
        ln -sf libz.so.1 libz.so
        cd ../../../../
fi
rm ont-guppy_${GUPPY_VER}_linux64.tar.gz
