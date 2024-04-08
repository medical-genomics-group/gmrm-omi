import numpy as np
from sklearn.metrics import r2_score
import argparse
import os
import csv

# This script calculates R2 metric
print("...Calculating R2 metric", flush=True)

# Initialize parser
parser = argparse.ArgumentParser()
parser.add_argument("-est", "--est", help = "Path to y estimates file")
parser.add_argument("-true", "--true", help = "Path to true phen file")
parser.add_argument("-out_name", "--out-name", help = "Output file name")
args = parser.parse_args()
est_fpath = args.est
phen_fpath = args.true
out_name = args.out_name

print("Input arguments:")
print("--est", est_fpath)
print("--true", phen_fpath)
print("--out-name", out_name)
print("\n", flush=True)

est_dirpath = os.path.dirname(est_fpath)

# load data from file
def load_file(path, col):
    y = []
    file = open(path, "rb")
    for row in file:
        row = row.split()
        y.append(float(row[col]))
    return np.array(y)

# There is only one column in y_est file
y_est = load_file(est_fpath, 0)

# true phen file is in plink format
y_true = load_file(phen_fpath, 2)

# calculate R2
r2 = r2_score(y_true, y_est)

print("R2 = %0.4f\n" % r2, flush=True)

fout = os.path.join(est_dirpath, out_name+'.csv')
print("...Saving R2 to CSV file")
print(fout)
print("\n", flush=True)

csv_file = open(fout, 'w', newline="")
csv_writer = csv.writer(csv_file, delimiter='\t')
row = [r2]
csv_writer.writerow(row)
csv_file.close()