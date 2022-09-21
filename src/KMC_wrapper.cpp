
#include "KMC_code.hh"
#include "globals.hh"
#include "SeqIO.hh"
#include "kmc_api/kmc_file.h"
#include "include/kmc_runner.h"

pair<string, int64_t> run_kmc(const vector<string>& input_files, LL k, LL n_threads, LL ram_gigas, int64_t min_abundance, int64_t max_abundance){

    write_log("Running KMC k-mer counter", LogLevel::MAJOR);

    string KMC_db_file_prefix = get_temp_file_manager().create_filename("kmers");

    KMC::Stage1Params stage1Params;

    string f = input_files[0]; // First input file
    SeqIO::FileFormat format = SeqIO::figure_out_file_format(f);

    for(string f2 : input_files){
        SeqIO::FileFormat format2 = SeqIO::figure_out_file_format(f2);
        if(format.format != format2.format || format.gzipped != format2.gzipped){
            throw std::runtime_error("Error: all input files must have the same format");
        }
    }

    stage1Params.SetInputFiles(input_files)
        .SetKmerLen(k)
        .SetNThreads(n_threads)
        .SetMaxRamGB(ram_gigas)
        .SetInputFileType(format.format == SeqIO::FASTA ? KMC::InputFileType::MULTILINE_FASTA : KMC::InputFileType::FASTQ)
        .SetCanonicalKmers(false)
        .SetTmpPath(get_temp_file_manager().get_dir());

    KMC::Runner kmc;

    auto stage1Results = kmc.RunStage1(stage1Params);

    uint32_t ramForStage2 = ram_gigas;
    KMC::Stage2Params stage2Params;
    stage2Params.SetNThreads(n_threads)
        .SetMaxRamGB(ramForStage2)
        .SetCutoffMin(min_abundance)
        .SetCutoffMax(max_abundance)
        .SetOutputFileName(KMC_db_file_prefix)
        .SetStrictMemoryMode(true);

    auto stage2Results = kmc.RunStage2(stage2Params);

    int64_t n_kmers = stage2Results.nUniqueKmers - stage2Results.nBelowCutoffMin - stage2Results.nAboveCutoffMax;

    // Clean up the KMC global singleton config state because it seems that it's left
    // in a partial state sometimes, which messes up our code if we call KMC again later.
    CConfig::GetInstance().input_desc.clear();
    CConfig::GetInstance().headers.clear();
    CConfig::GetInstance().simple_output_desc.clear();
    CConfig::GetInstance().transform_output_desc.clear();

    return {KMC_db_file_prefix, n_kmers};
}
