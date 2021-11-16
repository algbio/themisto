#pragma once

#include <gtest/gtest.h>
#include "../globals.hh"
#include "setup_tests.hh"

void check_sequence_reader_output(const vector<string>& seqs, LL mode, string fastafile){
    Sequence_Reader sr(fastafile, mode);
    for(string seq : seqs){
        ASSERT_FALSE(sr.done());        
        ASSERT_EQ(sr.get_next_query_stream().get_all(), seq);
    }
    ASSERT_TRUE(sr.done());
}

TEST(INPUT_PARSING, fasta_basic){
    vector<string> seqs = {"AAGTGCTGTANAYA","ACGTURYKMSWBDHVN-"};
    string fasta = ">\n" + seqs[0] + "\n>\n" + seqs[1] + "\n";
    logger << fasta << endl << seqs << endl;
    string filename = string_to_temp_file(fasta);
    check_sequence_reader_output(seqs, FASTA_MODE, filename);

}

TEST(INPUT_PARSING, fasta_multiple_lines){
    vector<string> seqs = {"AAGTGCTGTANAYA","ACGTURYKMSWBDHVN-"};
    string fasta;

    // Write 3 chars per line
    for(string seq : seqs){
        fasta += ">\n";
        for(LL i = 0; i < (LL)seq.size(); i += 3){
            fasta += seq.substr(i,3) + "\n";
        }
    }
    logger << fasta << endl << seqs << endl;
    string filename = string_to_temp_file(fasta);
    check_sequence_reader_output(seqs, FASTA_MODE, filename);
}

TEST(INPUT_PARSING, fasta_upper_case){
    vector<string> seqs = {"AagTGCtGTaNAYA","AcGTURYKmSWbDHVn-"};
    string fasta;

    for(string seq : seqs) fasta += ">\n" + seq + "\n";
    
    logger << fasta << endl << seqs << endl;
    string filename = string_to_temp_file(fasta);

    for(string& seq : seqs) for(char& c : seq) c = toupper(c); // Upper case for validation
    
    check_sequence_reader_output(seqs, FASTA_MODE, filename);
}

TEST(INPUT_PARSING, fasta_super_long_line){
    vector<string> seqs;
    seqs.push_back(string(1e6, 'A'));
    seqs.push_back(string(1e5, 'G'));

    string fasta;
    for(string seq : seqs) fasta += ">\n" + seq + "\n";
    string filename = string_to_temp_file(fasta);
    check_sequence_reader_output(seqs, FASTA_MODE, filename);
}

TEST(INPUT_PARSING, fasta_headers){
    vector<string> seqs;
    seqs.push_back(string(1e6, 'A'));
    seqs.push_back(string(1e5, 'G'));

    vector<string> headers;
    headers.push_back(string(1e5, 'h'));
    headers.push_back(string(1e6, 'H'));

    string fasta;
    for(LL i = 0; i < seqs.size(); i++) fasta += ">" + headers[i] + "\n" + seqs[i] + "\n";
    string filename = string_to_temp_file(fasta);
    Sequence_Reader sr(filename, FASTA_MODE);
    Read_stream rs = sr.get_next_query_stream();
    ASSERT_EQ(rs.header, headers[0]);
    rs.get_all();
    rs = sr.get_next_query_stream();
    ASSERT_EQ(rs.header, headers[1]);
    rs.get_all();
}



TEST(INPUT_PARSING, fastq_basic){
    vector<string> seqs =  {"AAGTGCTGTANAYAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA","ACGTURYKMSWBDHVN-"};
    vector<string> quals = {"!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~", "IIIIIIIIIIIIIIIII"};
    string fastq = "@\n" + seqs[0] + "\n+\n" + quals[0] + "\n" +
                   "@\n" + seqs[1] + "\n+\n" + quals[1] + "\n";
    logger << fastq << endl << seqs << " " << quals << endl;
    string filename = string_to_temp_file(fastq);
    check_sequence_reader_output(seqs, FASTQ_MODE, filename);

}


TEST(INPUT_PARSING, fastq_upper_case){
    vector<string> seqs =  {"AAGtGcygTANAynAAaAAAAAAAAAAAAAAAAAAAAAAAAAAaaaAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA","ACGTURYKMSWBDHVN-"};
    vector<string> quals = {"!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~", "IIIIIIIIIIIIIIIII"};
    string fastq = "@\n" + seqs[0] + "\n+\n" + quals[0] + "\n" +
                   "@\n" + seqs[1] + "\n+\n" + quals[1] + "\n";
    logger << fastq << endl << seqs << " " << quals << endl;
    string filename = string_to_temp_file(fastq);
    for(string& seq : seqs) for(char& c : seq) c = toupper(c); // Upper case for validation
    check_sequence_reader_output(seqs, FASTQ_MODE, filename);
}


TEST(INPUT_PARSING, fastq_super_long_line){
    vector<string> seqs;
    seqs.push_back(string(1e6, 'A'));
    seqs.push_back(string(1e5, 'G'));
    vector<string> quals;
    quals.push_back(string(1e6, 'I'));
    quals.push_back(string(1e5, 'I'));

    string fastq = "@\n" + seqs[0] + "\n+\n" + quals[0] + "\n" +
                   "@\n" + seqs[1] + "\n+\n" + quals[1] + "\n";
    string filename = string_to_temp_file(fastq);
    check_sequence_reader_output(seqs, FASTQ_MODE, filename);
}


TEST(INPUT_PARSING, fastq_headers){
    vector<string> seqs;
    seqs.push_back(string(1e6, 'A'));
    seqs.push_back(string(1e5, 'G'));

    vector<string> quals;
    quals.push_back(string(1e6, 'I'));
    quals.push_back(string(1e5, 'I'));

    vector<string> headers;
    headers.push_back(string(1e5, 'h'));
    headers.push_back(string(1e6, 'H'));

    string fastq = "@" + headers[0] + "\n" + seqs[0] + "\n+\n" + quals[0] + "\n" +
                   "@" + headers[1] + "\n" + seqs[1] + "\n+\n" + quals[1] + "\n";
    string filename = string_to_temp_file(fastq);
    Sequence_Reader sr(filename, FASTQ_MODE);
    Read_stream rs = sr.get_next_query_stream();
    ASSERT_EQ(rs.header, headers[0]);
    rs.get_all();
    rs = sr.get_next_query_stream();
    ASSERT_EQ(rs.header, headers[1]);
    rs.get_all();
}

TEST(INPUT_PARSING, fastq_things_after_plus){
    vector<string> seqs =  {"AAGTGCTGTANAYAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA","ACGTURYKMSWBDHVN-"};
    vector<string> quals = {"!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~", "IIIIIIIIIIIIIIIII"};
    string fastq = "@\n" + seqs[0] + "\n+SOMETHING\n" + quals[0] + "\n" +
                   "@\n" + seqs[1] + "\n+SOMETHING2\n" + quals[1] + "\n";
    logger << fastq << endl << seqs << " " << quals << endl;
    string filename = string_to_temp_file(fastq);
    check_sequence_reader_output(seqs, FASTQ_MODE, filename);
}


/*
TEST(INPUT_PARSING, fastq_multiple_lines){
    // We don't support multi-line sequences in FASTQ
}
*/
