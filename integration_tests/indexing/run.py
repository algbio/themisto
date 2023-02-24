import subprocess
import time
import sys
import os
import uuid
import signal
import random
import csv

if sys.version_info < (3, 0):
    sys.stdout.write("Error: Python3 required\n")
    sys.exit(1)

def run_get_output(command):
    # Command can have pipes
    sys.stderr.write(command + "\n")
    return subprocess.run(command, shell=True, stdout=subprocess.PIPE).stdout.decode("utf-8").strip()

def run(command):
    sys.stderr.write(command + "\n")
    return subprocess.run(command, shell=True).returncode

def check_outputs(themisto_outfile, ref_outfile):
    themisto_lines = open(themisto_outfile).read().splitlines()
    ref_lines = open(ref_outfile).read().splitlines()

    assert(len(themisto_lines) == len(ref_lines))
    for i in range(len(themisto_lines)):
        T = themisto_lines[i]
        R = ref_lines[i]
        assert(len(T) == len(R))
        assert(T[0] == R[0]) # K-mer string
        T_colors = map(int, T.split()[1:])
        R_colors = map(int, R.split()[1:])
        assert(sorted(T_colors) == sorted(R_colors)) # The colors are in arbitrary order so we sort
    print("\n===\nOK: {}\n===\n".format(themisto_outfile))

themisto_binary = "../../build/bin/themisto"
ref_binary = "../reference_implementation/themisto_reference_implementation"
temp_dir = "./temp"
out_dir = "./out"
infile_list = "file_list.txt"
manual_colorfile = "../ref_sequences/colorfile.txt"
concat_all_file = temp_dir + "/all.fasta.gz"

run("mkdir -p {}".format(out_dir))
run("mkdir -p {}".format(temp_dir))
assert(run("find ../ref_sequences -type f | grep fasta.gz > " + infile_list) == 0)
assert(run("cat " + infile_list +" | xargs cat > " + concat_all_file) == 0)

print("Input files:")
print(open(infile_list).read())

def build_index(k, d, input, rc, color_input_mode, outfile, color_set_type):
    print("Color set type", color_set_type)
    assert(run("{} build --n-threads 4 -k {} -i {} -o {} --temp-dir {} {} -d {} {} --coloring-structure-type {}".format(
        themisto_binary, k, input, outfile, temp_dir, "--reverse-complements" if rc else "", d, color_input_mode, color_set_type))
    == 0)

def dump_color_matrix(indexfile, outfile):
    assert(run("{} dump-color-matrix -i {} -o {} --sparse".format(
        themisto_binary, indexfile, outfile))
    == 0)

def dump_reference_color_matrix(k, inputfile, rc, color_input_mode, outfile):
    assert(run("{} dump-color-matrix -k {} -i {} -o {} {} {}".format(
        ref_binary, k, inputfile, outfile, "--rc" if rc else "", color_input_mode))
    == 0)

def test_file_colors_with_load_dbg(input):
    out_prefix = out_dir + "/load_dbg_test"

    # Build DBG without colors
    assert(run("{} build --n-threads 4 -k 31 -i {} -o {} --temp-dir {} --reverse-complements --no-colors".format(
        themisto_binary, input, out_prefix, temp_dir
    )) == 0)

    # Build colors with --file-colors
    assert(run("{} build --n-threads 4 -k 31 -i {} -o {} --temp-dir {} --reverse-complements --file-colors --load-dbg".format(
        themisto_binary, input, out_prefix, temp_dir
    )) == 0)

    # Compare to reference implementation
    dump_reference_color_matrix(31, input, True, "--file-colors", out_dir + "/load_dbg_test.colordump")
    dump_color_matrix(out_prefix, out_dir + "/load_dbg_test.colordump.ref")
    check_outputs(out_dir + "/load_dbg_test.colordump", out_dir + "/load_dbg_test.colordump.ref")
    

def test_transform_index(input):
    out1_prefix = out_dir + "/test_transform_index"
    out2_prefix = out_dir + "/test_transform_index"
    out3_prefix = out_dir + "/test_transform_index"

    # Build sdsl-hybrid index
    assert(run("{} build --n-threads 4 -k 31 -i {} -o {} --temp-dir {} --reverse-complements --coloring-structure-type sdsl-hybrid".format(
        themisto_binary, input, out1_prefix, temp_dir
    )) == 0)

    # Transform to Roaring
    assert(run("{} build --n-threads 4 --from-index {} -o {} --temp-dir {} --reverse-complements --coloring-structure-type roaring".format(
        themisto_binary, out1_prefix, out2_prefix, temp_dir
    )) == 0)

    # Transform back to sdsl-hybrid
    assert(run("{} build --n-threads 4 --from-index {} -o {} --temp-dir {} --reverse-complements --coloring-structure-type sdsl-hybrid".format(
        themisto_binary, out2_prefix, out3_prefix, temp_dir
    )) == 0)

    # Check that the color matrix dump matches with the dump from the original
    dump_color_matrix(out1_prefix, out1_prefix + ".colordump") 
    dump_color_matrix(out3_prefix, out3_prefix + ".colordump")
    check_outputs(out1_prefix + ".colordump", out3_prefix + ".colordump")

test_transform_index(infile_list)
test_file_colors_with_load_dbg(infile_list) # Tests --no-colors and --load-dbg


# Testing different coloring options, different k, diffenent d, rc or no-rc.
runs = [
    [31, 1, concat_all_file, False, "--manual-colors " + manual_colorfile, out_dir + "/manual-colors", "sdsl-hybrid"],
    [100, 8, concat_all_file, False,  "--manual-colors " + manual_colorfile, out_dir + "/manual-colors", "sdsl-hybrid"], # Large k
    [31, 2, concat_all_file, True,  "--manual-colors " + manual_colorfile, out_dir + "/manual-colors-rc", "sdsl-hybrid"],
    [31, 3, infile_list, False, "--sequence-colors", out_dir + "/seq-colors", "sdsl-hybrid"],
    [31, 4, infile_list, True,  "--sequence-colors", out_dir + "/seq-colors-rc", "sdsl-hybrid"],
    #[31, 5, infile_list, False, "--file-colors",     out_dir + "/file-colors", "sdsl-hybrid"], # Can't do file colors without rc because of how ggcat is called
    [31, 6, infile_list, True,  "--file-colors", out_dir + "/file-colors-rc", "sdsl-hybrid"],
    [31, 7, infile_list, True,  "--file-colors", out_dir + "/file-colors-rc", "roaring"], # Roaring
]

for k, d, input, rc, color_input_mode, outfile, color_set_type in runs:
    build_index(k, d, input, rc, color_input_mode, outfile, color_set_type)
    dump_color_matrix(outfile, outfile + ".colordump")
    dump_reference_color_matrix(k, infile_list, rc, color_input_mode, outfile + ".colordump.ref")
    check_outputs(outfile + ".colordump", outfile + ".colordump.ref")