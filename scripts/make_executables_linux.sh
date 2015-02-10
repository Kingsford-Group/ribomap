#!/bin/bash
cd ../src
make riboprof
make install
mkdir -p ../../ribomap-linux/
ribomap_dir=../../ribomap-linux/
scripts="run_ribomap.sh include_prerequisites.sh hela_ribo_analysis.sh offset.txt build_contaminant.py filter_gencode_transcript.py translation.py"
cp -r ../bin ${ribomap_dir}
cp -r ../lib ${ribomap_dir}
mkdir -p ${ribomap_dir}scripts
for s in ${scripts}; do
    echo "cp ../scripts/$s ${ribomap_dir}scripts"
done
cp ../README.md ${ribomap_dir}
cd ../../
tar -cvzf ribomap-linux.tar.gz ribomap-linux
rm -r ribomap-linux
