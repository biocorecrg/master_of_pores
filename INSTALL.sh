#!/bin/bash
#
# This script will download and install gouppy
# params:
# 1 guppy version

if [ -e "./mop_preprocess/bin/guppy_basecaller" ] ; then
	echo "unlinking previously installed versions"
	cd mop_preprocess/bin; find . -maxdepth 1 -type l | xargs rm; cd ../../
fi

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

wget https://biocore.crg.eu/public/mop3_pub/models.tar
mv models.tar mop_preprocess/guppy_models/
cd mop_preprocess/guppy_models; tar -xvf  models.tar; rm models.tar; cd ../../

mkdir -p mop_preprocess/bin/ont-guppy_${GUPPY_VER}
mv ont-guppy/* mop_preprocess/bin/ont-guppy_${GUPPY_VER}
for i in mop_preprocess/guppy_models/*.gz; do gzip -cd $i > mop_preprocess/bin/ont-guppy_${GUPPY_VER}/data/`basename $i .gz`; done
rmdir ont-guppy

cd mop_preprocess/bin
ln -s ont-guppy_${GUPPY_VER}/bin/guppy_* .
ln -s ont-guppy_${GUPPY_VER}/lib/* .
cd ../../

if [ ! -e "./mop_preprocess/bin/ont-guppy_${GUPPY_VER}/lib/libz.so" ] ; then
	if [ -e "./mop_preprocess/bin/ont-guppy_${GUPPY_VER}/lib/libz.so.1" ] ; then
        	unlink mop_preprocess/bin/ont-guppy_${GUPPY_VER}/lib/libz.so
        	cd mop_preprocess/bin/ont-guppy_${GUPPY_VER}/lib/
        	ln -s libz.so.1 libz.so
        	cd ../../../../
	fi
fi

rm ./ont-guppy_${GUPPY_VER}_linux64.tar.gz
