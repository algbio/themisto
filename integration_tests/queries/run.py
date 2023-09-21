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
        print(T)
        print(R)
        print("--")
        assert(len(T) == len(R))
        assert(T.split()[0] == R.split()[0]) # Sequence id
        T_hits = list(map(int, T.split()[1:]))
        R_hits = list(map(int, R.split()[1:]))
        assert(T_hits == sorted(R_hits)) # The reference hits are in arbitrary order so we sort
    print("OK:", themisto_outfile)

themisto_binary = "../../build/bin/themisto"
ref_binary = "../reference_implementation/themisto_reference_implementation"
k = 9
temp_dir = "./temp"
index_prefix = temp_dir + "/index"
query_file = "generated_queries.fasta.gz"
out_dir = "./out"

run("mkdir -p {}".format(out_dir))
run("mkdir -p {}".format(temp_dir))

assert(run("python3 gen_queries.py " + query_file) == 0)
assert(run("echo -n > file_list.txt") == 0) # Clear file
assert(run("find ../ref_sequences -type f | grep fasta.gz >> file_list.txt") == 0)

# Build index
assert(run("{} build -k {} -i {} -o {} --temp-dir {} --sequence-colors -d 5".format(
    themisto_binary, k, "file_list.txt", index_prefix, temp_dir)) 
== 0)

with open('parameters.csv') as csvfile:
    reader = csv.reader(csvfile, delimiter=',', quotechar='"')
    data = list(csv.reader(csvfile, delimiter=',', quotechar='"'))
    for line in data[1:]: # Skip the header line
        threshold, ignore, revcomp = line
        run_id = "{}-{}-{}".format(threshold, ignore, revcomp)
        themisto_outfile = out_dir + "/" + run_id + ".txt"
        ref_outfile = out_dir + "/ref-" + run_id + ".txt"

        # Query Themisto
        assert(run("{} pseudoalign -q {} -i {} -o {} --temp-dir {} --threshold {} {} {} --sort-hits --sort-output-lines -t {} --buffer-size-megas {}".format(
            themisto_binary, query_file, index_prefix, themisto_outfile, temp_dir, threshold, # Very small buffer to expose parallelism bugs
            "--include-unknown-kmers" if ignore == "no" else "", 
            "--rc" if revcomp == "yes" else "",
             4, 0.00001))
        == 0)

        # Query reference implementation
        assert(run("{} query -k {} -i {} -q {} -o {} --threshold {} {} {}".format(
            ref_binary, k, "file_list.txt", query_file, ref_outfile,
            threshold,
            "--include-unknown-kmers" if ignore == "no" else "", 
            "--rc" if revcomp == "yes" else ""))
        == 0)

        # Check that outputs match
        check_outputs(themisto_outfile, ref_outfile)

