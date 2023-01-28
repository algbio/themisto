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
    return subprocess.run(command, shell=True)

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
k = 31
temp_dir = "./temp"
out_dir = "./out"
infile_list = "file_list.txt"

run("mkdir -p {}".format(out_dir))
run("mkdir -p {}".format(temp_dir))
run("find ../ref_sequences -type f | grep fasta.gz > " + infile_list)

def build_index(k, input, rc, color_input_mode, outfile, color_set_type):
    print("Color set type", color_set_type)
    run("{} build --n-threads 4 -k {} -i {} -o {} --temp-dir {} {} -d 5 {} --coloring-structure-type {}".format(
        themisto_binary, k, input, outfile, temp_dir, "--reverse-complements" if rc else "", color_input_mode, color_set_type)
    )

def dump_color_matrix(indexfile, outfile):
    run("{} dump-color-matrix -i {} -o {} --sparse".format(
        themisto_binary, indexfile, outfile)
    )

def dump_reference_color_matrix(k, inputfile, rc, outfile):
    run("{} dump-color-matrix -k {} -i {} -o {} {}".format(
        ref_binary, k, inputfile, outfile, "--rc" if rc else "")
    )

runs = [
    [31, infile_list, False, "--sequence-colors", out_dir + "/seq-colors"],
    [31, infile_list, True,  "--sequence-colors", out_dir + "/seq-colors-rc"],
    #[31, infile_list, False, "--file-colors",     out_dir + "/file-colors"], # Can't do file colors without rc because of how ggcat is called
    [31, infile_list, True,  "--file-colors",     out_dir + "/file-colors-rc"]
]

for k, input, rc, color_input_mode, outfile in runs:
    for color_set_type in ["sdsl-hybrid", "roaring"]:
        build_index(k, input, rc, color_input_mode, outfile, color_set_type)
        dump_color_matrix(outfile, outfile + ".colordump")
        dump_reference_color_matrix(k, infile_list, rc, outfile + ".colordump.ref")
        check_outputs(outfile + ".colordump", outfile + ".colordump.ref")
