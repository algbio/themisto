import gzip
import random
import sys

random.seed(42)

rhinovirus_C = open("rhinovirus_C.txt").read().strip()
# ^ This is the first virus genome from ../ref_sequences/rhinovirus_c.fasta.gz

filename = sys.argv[1]
sample = list(rhinovirus_C[1000:3000]) # Slice of 2000bp from the middle

# Let's add mutations one by one
ACGT = "ACGT"

out = gzip.open(filename, 'wb')
for i in range(0,100):

	# Add a mutation
	p = random.randint(0,len(sample)-1)
	sample[p] = ACGT[random.randint(0,3)]

	# Write out
	out.write(">{}\n".format(i).encode("utf-8")) # Header
	out.write("{}\n".format("".join(sample)).encode("utf-8")) # Sequence
	
	
     
