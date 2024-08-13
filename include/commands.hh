#pragma once

int build_index_main(int argc, char** argv);
int pseudoalign_main(int argc, char** argv);
int extract_unitigs_main(int argc, char** argv);
int dump_index_main(int argc, char** argv);
int stats_main(int argc, char** argv);
int dump_color_matrix_main(int argc, char** argv);

int color_set_diagnostics_main(int argc, char** argv); // Undocumented developer feature
int make_d_equal_1_main(int argc, char** argv); // Undocumented developer feature
int dump_distinct_color_sets_to_binary_main(int argc, char** argv); // Undocumented developer feature


/*
int lookup_kmer_main(int argc, char** argv);
int lookup_color_main(int argc, char** argv);
*/