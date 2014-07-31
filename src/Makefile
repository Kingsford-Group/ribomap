COMPILE.c = $(CXX) $(CXXFLAGS)
CXX = g++
CXXFLAGS = -std=c++11 -O3 -Wno-deprecated-declarations -DSEQAN_HAS_ZLIB=1 -DSEQAN_HAS_BZIP2=1 
INC = -I/home/hw1/ribomap/include/ -I/opt/local/include
CXXFLAGS += $(INC)
LDFLAGS = -lz -lbz2 
OUTPUT_OPTION = -o $@

all: ribomap merge_fq_to_fa
clean:
	rm -f *.o ribomap merge_fq_to_fa

merge_fq_to_fa: preprocess_fq.cpp 
	$(COMPILE.c) $(OUTPUT_OPTION) $^
	chmod u+x $@

ribomap: read_mapper.o ribomap_profiler.o abundance_rank.o bam_parser.o gencode_parser.o utils.o
	$(CXX) $(OUTPUT_OPTION) $^ $(LDFLAGS)
	chmod u+x $@

utils.cpp: utils.hpp
gencode_parser.cpp: gencode_parser.hpp utils.hpp
bam_parser.cpp: gencode_parser.hpp ribomap_profiler.hpp bam_parser.hpp utils.hpp
ribomap_profiler.cpp: gencode_parser.hpp ribomap_profiler.hpp bam_parser.hpp
abundance_rank.cpp: gencode_parser.hpp ribomap_profiler.hpp abundance_rank.hpp
preprocess_fq.cpp: ezOptionParser.hpp
read_mapper.cpp: ezOptionParser.hpp gencode_parser.hpp ribomap_profiler.hpp bam_parser.hpp abundance_rank.hpp

