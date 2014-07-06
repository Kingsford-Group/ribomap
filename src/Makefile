compile_mapper: 
	/opt/local/stow/gcc-4.8.1/bin/g++ -std=c++11 -I/home/hw1/.local/include -lz -lbz2 -DSEQAN_HAS_ZLIB=1 -DSEQAN_HAS_BZIP2=1 -Wno-deprecated-declarations -Wno-virtual-move-assign -I/opt/local/include -L/opt/local/lib -lboost_system -lboost_filesystem -O3 -o ribomap read_mapper.cpp footprint_assigner.cpp tasep_ode.cpp bam_parser.cpp gencode_parser.cpp utils.cpp

compile_footprint_generator:
	/opt/local/stow/gcc-4.8.1/bin/g++ -std=c++11 -I/home/hw1/.local/include -lz -lbz2 -DSEQAN_HAS_ZLIB=1 -DSEQAN_HAS_BZIP2=1 -Wno-deprecated-declarations -Wno-virtual-move-assign -I/opt/local/include -L/opt/local/lib -lboost_system -lboost_filesystem -O3 -o footprint_generator footprint_generator.cpp footprint_assigner.cpp tasep_ode.cpp bam_parser.cpp gencode_parser.cpp utils.cpp

run_mapper:
	./ribomap /home/hw1/ribomap/tmp all_expressed_10m /data/iGenomes/Homo_sapiens/GenCode/gencode.v18.pc_transcripts_filter.fa /data/iGenomes/Homo_sapiens/GenCode/gencode.v18.annotation.gtf /home/hw1/riboseq/data/human_mouse_rp_Guo/sailfish_quant/quant.sf /home/hw1/ribomap/outputs 15

generate_all_expressed:
	./footprint_generator ../synthetic_data all_expressed /data/iGenomes/Homo_sapiens/GenCode/gencode.v18.pc_transcripts_filter.fa /data/iGenomes/Homo_sapiens/GenCode/gencode.v18.annotation.gtf ../sf_quant/quant.sf

generator_tester:
	./footprint_generator ../synthetic_data test_gene /data/iGenomes/Homo_sapiens/GenCode/gencode.v18.pc_transcripts_filter.fa /data/iGenomes/Homo_sapiens/GenCode/gencode.v18.annotation.gtf gene_ENSG00000145592.9.txt

