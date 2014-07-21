preprocess_fq:
	g++ -std=c++11 -O3 -o merge_fq_to_fa preprocess_fq.cpp
compile_mapper:
	/opt/local/stow/gcc-4.8.1/bin/g++ -std=c++11 -I/home/hw1/.local/include/ -lz -lbz2 -DSEQAN_HAS_ZLIB=1 -DSEQAN_HAS_BZIP2=1 -Wno-deprecated-declarations -Wno-virtual-move-assign -I/opt/local/include -L/opt/local/lib -O3 -o ribomap read_mapper.cpp ribomap_profiler.cpp bam_parser.cpp gencode_parser.cpp utils.cpp
