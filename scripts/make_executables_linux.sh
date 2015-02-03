#!/bin/bash
mkdir -p ../../ribomap-linux/
mkdir -p ../../ribomap-linux/data
mkdir -p ../../ribomap-linux/src
cp ribomap merge_fq_to_fa run_ribomap.sh get_data.sh transcript_codon_check.py transcript_filter.py ../../ribomap-linux/src
cp ../README.md ../../ribomap-linux/
cp ../data/rrna_human.fasta ../../ribomap-linux/data
cd ../../
tar -cvzf ribomap-linux.tar.gz ribomap-linux
rm -r ribomap-linux
