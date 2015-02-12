#!/bin/bash
os=$1
if [ "$os" = linux ]; then
    odir=../../ribomap-linux
elif [ "$os" = osx ]; then
    odir=../../ribomap-osx
else
    echo "Usage: ./make_executables.sh os_type"
    echo "os_type = [ linux | osx ]"
    exit
fi
cd ../src
make riboprof
make install
mkdir -p ${odir}
scripts="run_ribomap.sh include_prerequisites.sh hela_ribo_analysis.sh offset.txt build_contaminant.py filter_gencode_transcript.py translation.py"
cp -r ../bin ${odir}
cp -r ../lib ${odir}
mkdir -p ${odir}scripts
for s in ${scripts}; do
    cp ../scripts/$s ${odir}/scripts
done
cp ../README.md ${odir}
cd ../../
odir=${odir##*/}
tar -cvzf ${odir}.tar.gz ${odir}
rm -r ${odir}
