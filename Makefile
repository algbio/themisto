CXX = g++
STD = -std=c++14

.PHONY: tests build_index pseudoalign count_pseudoalignments statistics test_delta_coding count_hits_to_clusters KMC_test KMC_tester
all: tests build_index pseudoalign count_pseudoalignments statistics test_delta_coding count_hits_to_clusters KMC_test
essentials: tests build_index pseudoalign

libraries = BD_BWT_index/lib/*.a sdsl-lite/build/lib/libsdsl.a sdsl-lite/build/external/libdivsufsort/lib/libdivsufsort64.a
includes = -I BD_BWT_index/include -I sdsl-lite/include -I KMC -I KMC/kmer_counter/ 
KMC_objects = KMC/kmer_counter/mmer.o KMC/kmer_counter/mem_disk_file.o KMC/kmer_counter/rev_byte.o KMC/kmer_counter/bkb_writer.o KMC/kmer_counter/cpu_info.o KMC/kmer_counter/bkb_reader.o KMC/kmer_counter/fastq_reader.o KMC/kmer_counter/timer.o KMC/kmer_counter/develop.o KMC/kmer_counter/kb_completer.o KMC/kmer_counter/kb_storer.o KMC/kmer_counter/kmer.o KMC/kmer_counter/splitter.o KMC/kmer_counter/kb_collector.o KMC/kmer_counter/raduls_sse2.o KMC/kmer_counter/raduls_sse41.o KMC/kmer_counter/raduls_avx2.o KMC/kmer_counter/raduls_avx.o KMC/kmc_api/kmer_api.o KMC/kmc_api/kmc_file.o #KMC/kmer_counter/libs/libz.a KMC/kmer_counter/libs/libbz2.a KMC/kmer_counter/stdafx.o

tests : build/KMC_wrapper.o build/KMC.a 
	$(CXX) $(STD) tests.cpp build/KMC_wrapper.o build/KMC.a $(libraries) -lz -lbz2 -o tests -Wall -Wno-sign-compare -Wextra $(includes) -g -fopenmp -pthread

build_index: build/KMC_wrapper.o build/KMC.a
	$(CXX) $(STD) build_index.cpp build/KMC_wrapper.o build/KMC.a $(libraries) -lz -lbz2 -o build_index -Wall -Wno-sign-compare -Wextra $(includes) -g -O3 -fopenmp -pthread

pseudoalign:
	$(CXX) $(STD) pseudoalign.cpp build/KMC_wrapper.o build/KMC.a $(libraries) -lz -lbz2 -o pseudoalign -Wall -Wno-sign-compare -Wextra $(includes) -g -O3 -fopenmp -pthread
	
count_hits_to_clusters:
	$(CXX) $(STD) count_hits_to_clusters.cpp $(libraries) -o count_hits_to_clusters -Wall -Wno-sign-compare -Wextra $(includes) -g -O3 -pthread

count_pseudoalignments:
	$(CXX) $(STD) count_pseudoalignments.cpp $(libraries) -o count_pseudoalignments -Wall -Wno-sign-compare -Wextra $(includes) -g -O3 -pthread

statistics:
	$(CXX) $(STD) statistics.cpp $(libraries) -o statistics -Wall -Wno-sign-compare -Wextra $(includes) -g -O3 -pthread	

test_delta_coding:
	$(CXX) $(STD) test_delta_coding.cpp $(libraries) -o test_delta_coding -Wall -Wno-sign-compare -Wextra $(includes) -g -O3 -pthread

build/KMC.a: # NOT a .PHONY target
	ar r build/KMC.a $(KMC_objects)

build/KMC_counter.o: # NOT a .PHONY target
	$(CXX) $(STD) KMC/kmer_counter/kmer_counter.cpp -lz -lbz2 -o build/KMC_counter.o -c -Wall -Wno-sign-compare -Wno-implicit-fallthrough -Wextra $(includes) -g -fopenmp -pthread -no-pie

build/KMC_wrapper.o:  # NOT a .PHONY target
	$(CXX) $(STD) KMC_wrapper.cpp -lz -lbz2 -o build/KMC_wrapper.o -c -Wall -Wno-sign-compare -Wno-implicit-fallthrough -Wextra $(includes) -g -fopenmp -pthread -no-pie
	
KMC_test: build/KMC_wrapper.o build/KMC.a 
	$(CXX) $(STD) KMC_test.cpp build/KMC_wrapper.o build/KMC.a  $(libraries) -lz -lbz2 -o KMC_test -Wall -Wno-sign-compare -Wno-implicit-fallthrough -Wextra $(includes) -g -fopenmp -pthread -no-pie
	
KMC_tester: build/KMC_wrapper.o build/KMC.a 
	$(CXX) $(STD) KMC_tester.cpp build/KMC_wrapper.o build/KMC.a  $(libraries) -lz -lbz2 -o KMC_tester -Wall -Wno-sign-compare -Wno-implicit-fallthrough -Wextra $(includes) -g -fopenmp -pthread -no-pie

