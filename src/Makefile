CXX = g++
CXXFLAGS = -std=c++11 -O3 -Wno-deprecated-declarations -DSEQAN_HAS_ZLIB=1 -DSEQAN_HAS_BZIP2=1 # -fopenmp
INC = -I/home/hw1/.local/include/ -I/opt/local/include
#LIB = -L/opt/local/lib -lboost_system -lboost_filesystem
CXXFLAGS += $(INC)
COMPILE.c = $(CXX) $(CXXFLAGS)
LDFLAGS = -lz -lbz2 # -pthread -fopenmp -lpthread
OUTPUT_OPTION = -o $@

all: riboprof
install:
	mkdir -p ../bin/
	cp riboprof ../bin/
clean:
	rm -f *.o riboprof ref_checker footprint_generator

riboprof: read_mapper.o ribomap_profiler.o bam_parser.o reference_info_builder.o utils.o abundance_rank.o profile_writer.o
	$(COMPILE.c) $(OUTPUT_OPTION) $^ $(LDFLAGS)
	chmod u+x $@

footprint_generator: footprint_generator.o reference_info_builder.o ribomap_profiler.o transcript_model.o utils.o tasep_ode.o
	$(COMPILE.c) $(OUTPUT_OPTION) $^ $(LDFLAGS)
	chmod u+x $@

ref_checker: reference_info_builder.o utils.o
	$(COMPILE.c) $(OUTPUT_OPTION) $^ $(LDFLAGS)
	chmod u+x $@

read_mapper.o: read_mapper.cpp ezOptionParser.hpp reference_info_builder.hpp ribomap_profiler.hpp bam_parser.hpp utils.hpp abundance_rank.hpp profile_writer.hpp
utils.o: utils.cpp utils.hpp 
reference_info_builder.o: reference_info_builder.cpp reference_info_builder.hpp utils.hpp
bam_parser.o: bam_parser.cpp reference_info_builder.hpp ribomap_profiler.hpp bam_parser.hpp utils.hpp
ribomap_profiler.o: ribomap_profiler.cpp reference_info_builder.hpp ribomap_profiler.hpp bam_parser.hpp
abundance_rank.o: abundance_rank.cpp reference_info_builder.hpp ribomap_profiler.hpp abundance_rank.hpp
profile_writer.o: profile_writer.cpp profile_writer.hpp ribomap_profiler.hpp reference_info_builder.hpp utils.hpp
footprint_generator.o: footprint_generator.cpp reference_info_builder.hpp ribomap_profiler.hpp transcript_model.hpp utils.hpp
