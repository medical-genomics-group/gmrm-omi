import numpy as np
import argparse
import struct
from array import array
import os

# This script calculate posterior inclusion probabilities from binary .bet file

print("...Calculating Posterior Inclusion Probabilities", flush=True)

# Initialize parser
parser = argparse.ArgumentParser()
parser.add_argument("-bet", "--bet", help = "Path to beta binary file")
parser.add_argument("-out-name", "--out-name", help = "Output file name")
parser.add_argument("-iterations", "--iterations", help = "MCMC iterations range start:end")
args = parser.parse_args()

# Input parameters
betfile = args.bet
out_name = args.out_name
iterations = args.iterations

print("Input arguments: \n")
print("--bet", betfile)
print("--out-name", out_name)
print("--iterations", iterations)
print("\n", flush=True)

# Parse iterations range
it_start = int(iterations.split(":")[0])
it_end = int(iterations.split(":")[1])
iterations = (it_end - it_start) + 1 # total number of iterations

# Get basename of file from input bet file
basename = os.path.basename(betfile)
basename = basename.split('.')[0]
dirpath = os.path.dirname(betfile)

# open beta file
f = open(betfile, "rb")

# read first 4 bytes. Number of markers in bet file
buffer = f.read(4)

# unpack to int variable.
[m] = struct.unpack('I', buffer)

# Posterior Inclusion Probability 
pip = np.zeros(m)

# for each iterations to it_end
for i in range(1, it_end + 1):
    # read and unpack 4 bytes. Iteration number
    buffer = f.read(4)
    [it] = struct.unpack('I', buffer)

    # read all markers in current iteration
    buffer = f.read(m * 8)

    # if current iteration is in range
    if it >= it_start:
        # unpack vector of betas as doubles
        beta = struct.unpack(str(m)+'d', buffer)
        beta = np.array(beta)

        # sum all non zero betas
        beta[np.abs(beta) > 0] = 1
        pip += beta
# number of non zero betas devided by number of iterations. This results in probability of including marker to the model
pip /= iterations

# save pip to file
print("...Saving Posterior Inclusion Probabilities to binary file")
out_fpath = os.path.join(dirpath, "%s.bin" % out_name)
print(out_fpath, flush=True)

output_file = open(out_fpath, 'wb')
float_array = array('d', pip)
float_array.tofile(output_file)
output_file.close()