#include "stdafx.h"
/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 3.1.1
  Date   : 2019-05-19
*/

#include "kmc_api/kmc_file.h"
#include <fstream>
#include <string>
#include <vector>
#include <time.h>
#include <functional>
#include "timer.h"
#include "kmc.h"
//#include "meta_oper.h"
//#include "globals.hh"
#include "Argv.hh"
#include "TempFileManager.hh"

using namespace std;

extern Temp_File_Manager temp_file_manager;

uint64 total_reads, total_fastq_size;

void usage();
bool parse_parameters(int argc, char *argv[]);

CKMCParams Params;

//----------------------------------------------------------------------------------
// Application class
// Template parameters:
//    * SIZE     - maximal size of the k-mer (divided by 32)
template<unsigned SIZE> class CApplication
{
	CApplication<SIZE - 1> *app_1;
	CKMC<SIZE> *kmc;
	int p_k;
	bool is_selected;

public:
	CApplication(CKMCParams &Params) {
		p_k = Params.p_k;
		is_selected = p_k <= (int32) SIZE * 32 && p_k > ((int32) SIZE-1)*32;

		app_1 = new CApplication<SIZE - 1>(Params);
		if(is_selected)
		{			
			kmc = new CKMC<SIZE>;
			kmc->SetParams(Params);
		}
		else
		{
			kmc = nullptr;
		}
	};
	~CApplication() {
		delete app_1;
		if (kmc)
			delete kmc;
	}

	void SaveStatsInJSON(bool was_small_k_opt)
	{
		if (is_selected)
			kmc->SaveStatsInJSON(was_small_k_opt);
		else
			app_1->SaveStatsInJSON(was_small_k_opt);
	}

	void GetStats(double &time1, double &time2, double &time3, uint64 &_n_unique, uint64 &_n_cutoff_min, uint64 &_n_cutoff_max, uint64 &_n_total, uint64 &_n_reads, uint64 &_tmp_size, uint64 &_tmp_size_strict_mem, uint64 &_max_disk_usage, uint64& _n_total_super_kmers, bool& _was_small_k_opt) {
		if (is_selected)
		{
			kmc->GetStats(time1, time2, time3, _n_unique, _n_cutoff_min, _n_cutoff_max, _n_total, _n_reads, _tmp_size, _tmp_size_strict_mem, _max_disk_usage, _n_total_super_kmers, _was_small_k_opt);
		}
		else
			app_1->GetStats(time1, time2, time3, _n_unique, _n_cutoff_min, _n_cutoff_max, _n_total, _n_reads, _tmp_size, _tmp_size_strict_mem, _max_disk_usage, _n_total_super_kmers, _was_small_k_opt);
	}

	bool Process() {
		if (is_selected)
		{
			return kmc->Process();
		}
		else
			return app_1->Process();
	}
};

//----------------------------------------------------------------------------------
// Specialization of the application class for the SIZE=1
template<> class CApplication<1>
{
	CKMC<1> *kmc;
	int p_k;
	bool is_selected;

public:
	CApplication(CKMCParams &Params) {
		is_selected = Params.p_k <= 32;
		if(is_selected)
		{
			kmc = new CKMC<1>;
			kmc->SetParams(Params);
		}
		else
		{
			kmc = nullptr;
		}
	};
	~CApplication() {
		if(kmc)
			delete kmc;
	};

	void GetStats(double &time1, double &time2, double &time3, uint64 &_n_unique, uint64 &_n_cutoff_min, uint64 &_n_cutoff_max, uint64 &_n_total, uint64 &_n_reads, uint64 &_tmp_size, uint64 &_tmp_size_strict_mem, uint64 &_max_disk_usage, uint64& _n_total_super_kmers, bool& _was_small_k_opt) {
		if (is_selected)
		{
			if(kmc)
				kmc->GetStats(time1, time2, time3, _n_unique, _n_cutoff_min, _n_cutoff_max, _n_total, _n_reads, _tmp_size, _tmp_size_strict_mem, _max_disk_usage, _n_total_super_kmers, _was_small_k_opt);
		}
	}

	void SaveStatsInJSON(bool was_small_k_opt)
	{
		if (is_selected)
			kmc->SaveStatsInJSON(was_small_k_opt);
	}

	bool Process() {
		if (is_selected)
		{
			return kmc->Process();
		}
		return false;
	}
};


//----------------------------------------------------------------------------------
// Show execution options of the software
void usage()
{
	cout << "K-Mer Counter (KMC) ver. " << KMC_VER << " (" << KMC_DATE << ")\n"
		 << "Usage:\n kmc [options] <input_file_name> <output_file_name> <working_directory>\n"
		 << " kmc [options] <@input_file_names> <output_file_name> <working_directory>\n"
		 << "Parameters:\n"
		 << "  input_file_name - single file in specified (-f switch) format (gziped or not)\n"
		 << "  @input_file_names - file name with list of input files in specified (-f switch) format (gziped or not)\n"
		 << "Options:\n"
		 << "  -v - verbose mode (shows all parameter settings); default: false\n"
		 << "  -k<len> - k-mer length (k from " << MIN_K << " to " << MAX_K << "; default: 25)\n"
		 << "  -m<size> - max amount of RAM in GB (from 1 to 1024); default: 12\n"
		 << "  -sm - use strict memory mode (memory limit from -m<n> switch will not be exceeded)\n"
		 << "  -p<par> - signature length (5, 6, 7, 8, 9, 10, 11); default: 9\n"
		 << "  -f<a/q/m/bam> - input in FASTA format (-fa), FASTQ format (-fq), multi FASTA (-fm) or BAM (-fbam); default: FASTQ\n"	
		 << "  -ci<value> - exclude k-mers occurring less than <value> times (default: 2)\n"
		 << "  -cs<value> - maximal value of a counter (default: 255)\n"
		 << "  -cx<value> - exclude k-mers occurring more of than <value> times (default: 1e9)\n"
		 << "  -b - turn off transformation of k-mers into canonical form\n"
		 << "  -r - turn on RAM-only mode \n"
		 << "  -n<value> - number of bins \n"
		 << "  -t<value> - total number of threads (default: no. of CPU cores)\n"
		 << "  -sf<value> - number of FASTQ reading threads\n"
		 << "  -sp<value> - number of splitting threads\n"
		 << "  -sr<value> - number of threads for 2nd stage\n"
		 << "  -j<file_name> - file name with execution summary in JSON format\n"
		 << "  -w - without output\n"
		 << "Example:\n"
		 << "kmc -k27 -m24 NA19238.fastq NA.res /data/kmc_tmp_dir/\n"
	     << "kmc -k27 -m24 @files.lst NA.res /data/kmc_tmp_dir/\n";
}

bool CanCreateFile(const string& path)
{
	FILE* f = fopen(path.c_str(), "wb");
	if (!f)
		return false;
	fclose(f);
	remove(path.c_str());
	return true;
}
bool CanCreateFileInPath(const string& path)
{
	static const string name = "kmc_test.bin"; //Some random name
	if (path.back() == '\\' || path.back() == '/')
		return CanCreateFile(path + name);
	else
		return CanCreateFile(path + '/' + name);
}
//----------------------------------------------------------------------------------
// Parse the parameters
bool parse_parameters(int argc, char *argv[])
{
	int i;
	int tmp;

	if(argc < 4)
		return false;

	for(i = 1 ; i < argc; ++i)
	{
		if(argv[i][0] != '-')
			break;
		// Number of threads
		if(strncmp(argv[i], "-t", 2) == 0)
			Params.p_t = atoi(&argv[i][2]);
//		else 
		// k-mer length
		if (strncmp(argv[i], "-k", 2) == 0)
		{
			tmp = atoi(&argv[i][2]);
			if (tmp < MIN_K || tmp > MAX_K)
			{
				cerr << "Wrong parameter: k must be from range <" << MIN_K << "," << MAX_K << ">\n";
				return false;
			}
			else
				Params.p_k = tmp;
		}
		// Memory limit
		else if (strncmp(argv[i], "-m", 2) == 0)
		{
			tmp = atoi(&argv[i][2]);
			if (tmp < MIN_MEM)
			{
				cerr << "Wrong parameret: min memory must be at least " << MIN_MEM << "GB\n";
				return false;
			}
			else
				Params.p_m = tmp;
		}
		// Minimum counter threshold
		else if (strncmp(argv[i], "-ci", 3) == 0)
			Params.p_ci = atoi(&argv[i][3]);
		// Maximum counter threshold
		else if (strncmp(argv[i], "-cx", 3) == 0)
			Params.p_cx = atoll(&argv[i][3]);
		// Maximal counter value
		else if (strncmp(argv[i], "-cs", 3) == 0)
			Params.p_cs = atoll(&argv[i][3]);		
		// Set p1
		else if (strncmp(argv[i], "-p", 2) == 0)
		{
			tmp = atoi(&argv[i][2]);
			if (tmp < MIN_SL || tmp > MAX_SL)
			{
				cerr << "Wrong parameter: p must be from range <" << MIN_SL << "," << MAX_SL << ">\n";
				return false;
			}
			else
				Params.p_p1 = tmp;
		}
		// FASTA input files
		else if (strncmp(argv[i], "-fa", 3) == 0)
			Params.p_file_type = fasta;
		// FASTQ input files
		else if (strncmp(argv[i], "-fq", 3) == 0)
			Params.p_file_type = fastq;
		else if (strncmp(argv[i], "-fm", 3) == 0)
			Params.p_file_type = multiline_fasta;
		else if (strncmp(argv[i], "-fbam", 5) == 0)
			Params.p_file_type = bam;
#ifdef DEVELOP_MODE
		else if (strncmp(argv[i], "-vl", 3) == 0)
			Params.p_verbose_log = true;
#endif
		else if (strncmp(argv[i], "-v", 2) == 0)
			Params.p_verbose = true;
		else if (strncmp(argv[i], "-sm", 3) == 0 && strlen(argv[i]) == 3)
			Params.p_strict_mem = true;
		else if (strncmp(argv[i], "-r", 2) == 0)
			Params.p_mem_mode = true;
		else if(strncmp(argv[i], "-b", 2) == 0)
			Params.p_both_strands = false;
		// Number of reading threads
		else if(strncmp(argv[i], "-sf", 3) == 0)
		{
			tmp = atoi(&argv[i][3]);
			if(tmp < MIN_SF || tmp > MAX_SF)
			{
				cerr << "Wrong parameter: number of reading thread must be from range <" << MIN_SF << "," << MAX_SF << ">\n";
				return false;
			}
			else
				Params.p_sf = tmp;
		}
		// Number of splitting threads
		else if(strncmp(argv[i], "-sp", 3) == 0)
		{
			tmp = atoi(&argv[i][3]);
			if(tmp < MIN_SP || tmp > MAX_SP)
			{
				cerr << "Wrong parameter: number of splitting threads must be in range <" << MIN_SP << "," << MAX_SP << "<\n";
				return false;
			}
			else
				Params.p_sp = tmp;
		}		
		// Number of internal threads per 2nd stage
		else if(strncmp(argv[i], "-sr", 3) == 0)
		{
			tmp = atoi(&argv[i][3]);
			if(tmp < MIN_SR || tmp > MAX_SR)
			{
				cerr << "Wrong parameter: number of threads for 2nd stage must be in range <" << MIN_SR << "," << MAX_SR << "\n";
				return false;
			}
			else
				Params.p_sr = tmp;
		}
		else if (strncmp(argv[i], "-n", 2) == 0)
		{
			tmp = atoi(&argv[i][2]);
			if (tmp < MIN_N_BINS || tmp > MAX_N_BINS)
			{
				cerr << "Wrong parameter: number of bins must be in range <" << MIN_N_BINS << "," << MAX_N_BINS << "\n";
				return false;
			}
			else
				Params.p_n_bins = tmp;
		}
		else if (strncmp(argv[i], "-j", 2) == 0)
		{
			Params.json_summary_file_name = &argv[i][2];
			if (Params.json_summary_file_name == "")
				cerr << "Warning: file name for json summary file missed (-j switch)\n";			
		}
		else if (strncmp(argv[i], "-w", 2) == 0)	
			Params.p_without_output = true;		

		if (strncmp(argv[i], "-smso", 5) == 0)
		{
			tmp = atoi(&argv[i][5]);
			if (tmp < MIN_SMSO || tmp > MAX_SMSO)
			{
				cerr << "Wrong parameter: number of sorting threads per sorter in strict memory mode must be in range <" << MIN_SMSO << "," << MAX_SMSO << "\n";
				return false;
			}
			else
				Params.p_smso = tmp;
		}
		if (strncmp(argv[i], "-smun", 5) == 0)
		{
			tmp = atoi(&argv[i][5]);
			if (tmp < MIN_SMUN || tmp > MAX_SMUN)
			{
				cerr << "Wrong parameter: number of uncompactor threads in strict memory mode must be in range <" << MIN_SMUN << "," << MAX_SMUN << "\n";
				return false;
			}
			else
				Params.p_smun = tmp;
		}		
		if (strncmp(argv[i], "-smme", 5) == 0)
		{
			tmp = atoi(&argv[i][5]);
			if (tmp < MIN_SMME || tmp > MAX_SMME)
			{
				cerr << "Wrong parameter: number of merger threads in strict memory mode must be in range <" << MIN_SMME << "," << MAX_SMME << "\n";
				return false;
			}
			else
				Params.p_smme = tmp;
		}
	}

	if(argc - i < 3)
		return false;

	string input_file_name = string(argv[i++]);
	Params.output_file_name = string(argv[i++]);
	Params.working_directory = string(argv[i++]);

	Params.input_file_names.clear();
	if(input_file_name[0] != '@')
		Params.input_file_names.push_back(input_file_name);
	else
	{
		ifstream in(input_file_name.c_str()+1);
		if(!in.good())
		{
			cerr << "Error: No " << input_file_name.c_str()+1 << " file\n";
			return false;
		}

		string s;
		while(getline(in, s))
			if(s != "")
				Params.input_file_names.push_back(s);

		in.close();
		random_shuffle(Params.input_file_names.begin(), Params.input_file_names.end());
	}

	if (Params.p_t > Params.p_m * 64)
	{
		Params.p_t = Params.p_m * 64;
		cerr << "Warning: number of threads is reduced to " << Params.p_t << " (maximun numer of threads equals 64 * value of the -m parameter)\n";		
	}
	//Validate and resolve conflicts in parameters
	if (Params.p_strict_mem && Params.p_mem_mode)
	{
		cerr << "Error: -sm can not be used with -r\n";
		return false;
	}
	
	if (Params.p_k > 9)
	{
		if ((uint64)Params.p_cx > ((1ull << 32) - 1))
		{
			cerr << "Warning: for k > 9 maximum value of -cx is 4294967295\n";
			Params.p_cx = 4294967295;
		}
		if ((uint64)Params.p_cs > ((1ull << 32) - 1))
		{
			cerr << "Warning: for k > 9 maximum value of -cs is 4294967295\n";
			Params.p_cs = 4294967295;
		}
	}

	//Check if output files may be created and if it is possible to create file in specified tmp location
	if(!Params.p_without_output)
	{
		string pre_file_name = Params.output_file_name + ".kmc_pre";
		string suff_file_name = Params.output_file_name + ".kmc_suf";
		if (!CanCreateFile(pre_file_name))
		{
			cerr << "Error: Cannot create file: " << pre_file_name << "\n";
			return false;
		}
		if (!CanCreateFile(suff_file_name))
		{
			cerr << "Error: Cannot create file: " << suff_file_name << "\n";
			return false;
		}
	}
	if (!CanCreateFileInPath(Params.working_directory))
	{
		cerr << "Error: Cannot create file in specified working directory: " << Params.working_directory << "\n";
		return false;
	}
	return true;
}

//----------------------------------------------------------------------------------
// Check if --help or --version was used
bool help_or_version(int argc, char** argv)
{
	const string version = "--version";
	const string help = "--help";
	for (int i = 1; i < argc; ++i)
	{
		if (argv[i] == version || argv[i] == help)
			return true;
	}
	return false;
}

//----------------------------------------------------------------------------------
// Main function
void call_kmc(int argc, _TCHAR* argv[])
{

#ifdef WIN32
	_setmaxstdio(2040);
#endif

	if(!parse_parameters(argc, argv))
	{
		cerr << "Error calling KMC" << endl;
        exit(1);
	}
	
	CApplication<KMER_WORDS> *app = new CApplication<KMER_WORDS>(Params);

	if(!app->Process())
	{
		cerr << "KMC: Not enough memory or some other error" << endl;
		delete app;
		exit(1);
	}
}

// Only works for alphabet {a,c,g,t,A,C,G,T}
void KMC_wrapper(int64_t k, int64_t ram_gigas, int64_t n_threads, string fastafile, string outfile, string tempdir){

	// Check that the alphabet is {a,c,g,t,A,C,G,T} (otherwise k-mers would be dropped silently)
	ifstream fasta_input(fastafile);
	string line;
	while(getline(fasta_input,line)){
		if(line.size() > 0 && line[0] != '>'){
			for(char c : line){
				if(c != 'a' && c != 'c' && c != 'g' && c != 't' &&
				   c != 'A' && c != 'C' && c != 'G' && c != 'T'){
					cerr << "Error calling KMC: all characters in the sequences in the fasta-file must be from the alphabet {a,c,g,t,A,C,G,T}. Found invalid character: " << c << endl;
					exit(1);
				}
			}
		}
	}

	string KMC_database_file = temp_file_manager.get_temp_file_name("KMC");
	
	vector<string> argv = {"kmc", 
	                      "-fm", 
						  "-k" + to_string(k),
						  "-b",
						  "-m" + to_string(ram_gigas),
						  "-ci1", // Exclude k-mers occurring less than 1 times = do not exclude anything
						  "-cs1", // Maximum value of counter: 1
						  "-cx4294967295", // Excluding k-mers that occur more than times 4294967295 (the largest value it supports)
						  "-t" + to_string(n_threads),
						  fastafile, 
						  KMC_database_file, 
						  tempdir};

	cerr << "Calling KMC with: ";
	for(string S : argv) cerr << S << " ";
	cerr << endl;
	
	Argv A(argv);
    call_kmc(argv.size(), A.array);

	cerr << "Dumping k-mers to disk" << endl;

	// Dump the database to text
	ofstream out(outfile);

	CKMCFile kmer_database;
	kmer_database.OpenForListing(KMC_database_file);

	CKmerAPI kmer_object(k);
	string str;
	float counter;
	
	while(kmer_database.ReadNextKmer(kmer_object, counter)){
		kmer_object.to_string(str);
		out << str << "\n";
	}

}

