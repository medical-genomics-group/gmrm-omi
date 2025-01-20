import numpy as np
import argparse
import os
import struct
import random

print("---- Simulating example i.i.d. data ----\n", flush=True)

# Initialize parser
parser = argparse.ArgumentParser()
parser.add_argument("-out_dir", "--out-dir", help = "Output directory")
parser.add_argument("-out_name", "--out-name", help="Output file name", default="example")
parser.add_argument("-N", "--N", help="Number of individuals", default=1000)
parser.add_argument("-M", "--M", help="Number of markers", default=2000)
parser.add_argument("-C", "--C", help="Number of covariates", default=2)
parser.add_argument("-lam", "--lam", help="Sparsity", default=0.1)
parser.add_argument("-epi_var", "--epi-var", help="Epigenetic component variance", default=0.5)
parser.add_argument("-cov_var", "--cov-var", help="Covariate component variance", default=0.2)
args = parser.parse_args()

out_name = args.out_name
out_dir = args.out_dir
N = int(args.N)
M = int(args.M)
C = int(args.C)
lam = float(args.lam)
epi_var = float(args.epi_var)
cov_var = float(args.cov_var)

print("Input arguments:")
print("--out_dir", out_dir)
print("--out_name", out_name)
print("--N", N)
print("--M", M)
print("--C", C)
print("--lam", lam)
print("--cov-var", cov_var)
print("--epi-var", epi_var, flush=True)

print("\n...Simulating design matrix", flush=True)
X = np.random.normal(0,1,N*M).reshape((N,M))
Z = np.random.normal(0,1,N*C).reshape((N,C)) #covariate matrix

print("\n...Simulating marker effects", flush=True)
CM = int(M * lam) # number of cuasal markers
sigma2 = epi_var / CM
idx = random.sample(range(M), CM)
beta = np.zeros(M)
beta[idx] = np.random.normal(0.0, np.sqrt(sigma2), CM)
print("sigma2 = ", sigma2, flush=True)

print("\n...Simulating covariate effects", flush=True)
sigma2_cov = cov_var / C
delta = np.random.normal(0.0, np.sqrt(sigma2_cov), C)

print("\n...Computing outcome variable", flush=True)
g = np.matmul(X,beta)
c = np.matmul(Z,delta)
y = g + c + np.random.normal(0, np.sqrt( 1 - np.var(g + c)),  N) # adding Gaussian noise

print("Var(y) = ", np.var(y))
print("Var(g) = ", np.var(g))
print("Var(c) = ", np.var(c))

# Saving data in binary format
print("\n...Saving design matrix to bin file")
bin_fpath = os.path.join(out_dir, "%s.bin" % out_name)
print(bin_fpath, flush=True)
binf = open(bin_fpath, "wb")
b = struct.pack(str(N*M)+'d', *X.transpose().ravel().squeeze())
binf.write(b)
binf.close()

# Saving phenotype data
print("\n...Saving outcome to .phen file")
phen_fpath = os.path.join(out_dir, "%s.phen" % out_name)
print(phen_fpath, flush=True)
phenf = open(phen_fpath, "w")
for i, pheno in enumerate(y):
    line = "%d %d %0.10f\n" % (i, i, pheno)
    phenf.write(line)
phenf.close()

# Saving covariates
print("\n...Saving covariate matrix to .cov file")
cov_fpath = os.path.join(out_dir, "%s.cov" % out_name)
print(cov_fpath, flush=True)
covf = open(cov_fpath, "w")

header = "IID FID"
for i in range(C):
    header += " cov"+str(i+1)
covf.write(header+"\n")

for i in range(N):
    line = "%d %d" % (i, i)
    for j in range(C):
        line += " %0.10f" % (Z[i,j])
    line+="\n"
    covf.write(line)
covf.close()

# Save true signals to file
print("\n...Saving true signals to bin file")
ts_bin_fpath = os.path.join(out_dir, "%s_ts.bin" % out_name)
print(ts_bin_fpath, flush=True)
ts_binf = open(ts_bin_fpath, "wb")
ts_binf.write(struct.pack(str(M)+'d', *beta.squeeze()))
ts_binf.close()

# Save true covariate effects to file
print("\n...Saving true covariate effects to file")
tc_fpath = os.path.join(out_dir, "%s_true_cov.csv" % out_name)
print(tc_fpath, flush=True)
tcf = open(tc_fpath, "w")
for d in delta:
    line = "%0.10f\n" % (d)
    tcf.write(line)
tcf.close()