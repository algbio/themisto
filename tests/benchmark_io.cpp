#include "../include/test_tools.hh"
#include "../include/ReadBatch.hh"
#include "../include/input_reading.hh"
#include <vector>
#include <string>
#include <fstream>

using namespace std;

int main(){

    get_temp_file_manager().set_dir("./temp");

    // Generate 10 megabases of reads of length 100
    vector<string> reads;
    int64_t total_len = 0;
    for(int64_t i = 0; i < 100000; i++){
        reads.push_back(get_random_dna_string(100, 4));
        total_len += 100;
    }

    // Create a new fasta file for each so avoid caching effects
    string fastq1 = get_temp_file_manager().create_filename();
    string fastq2 = get_temp_file_manager().create_filename();
    string fastq3 = get_temp_file_manager().create_filename();

    write_as_fastq(reads, fastq1);
    write_as_fastq(reads, fastq2);
    write_as_fastq(reads, fastq3);

    // Unbuffered

    Sequence_Reader sr(fastq1, FASTQ_MODE);
    int64_t sr_t0 = cur_time_millis();
    while(!sr.done()){
        sr.get_next_query_stream().get_all();
    }
    int64_t sr_t1 = cur_time_millis();
    double sr_seconds = (sr_t1 - sr_t0) / 1000.0;
    cout << total_len / 1e6 / sr_seconds << " Mbp/s" << endl;


    // Buffered

    Sequence_Reader_Buffered srb(fastq2, FASTQ_MODE);
    int64_t srb_t0 = cur_time_millis();
    while(true){
        if(srb.get_next_read_to_buffer() == 0) break;
    }
    int64_t srb_t1 = cur_time_millis();
    double srb_seconds = (srb_t1 - srb_t0) / 1000.0;
    cout << total_len / 1e6 / srb_seconds << " Mbp/s" << endl;

}