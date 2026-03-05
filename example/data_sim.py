import numpy as np
import argparse
import os
import struct
import random
from scipy.stats import norm
import pandas as pd

print("---- Simulating example i.i.d. data ----\n", flush=True)

# Initialize parser
parser = argparse.ArgumentParser()
parser.add_argument("-out_dir", "--out-dir", help = "Output directory")
parser.add_argument("-out_name", "--out-name", help="Output file name", default="example")
parser.add_argument("-N", "--N", help="Number of individuals", default=2000)
parser.add_argument("-N_test", "--N-test", help="Number of test individuals", default=500)
parser.add_argument("-M", "--M", help="Number of markers", default=2000)
parser.add_argument("-C", "--C", help="Number of covariates", default=2)
parser.add_argument("-lam", "--lam", help="Sparsity", default=0.1)
parser.add_argument("-prevalence", "--prevalence", help="Prevalence", default=0.5)
parser.add_argument("-epi_var", "--epi-var", help="Epigenetic component variance", default=0.8)
parser.add_argument("-cov_var", "--cov-var", help="Covariate component variance", default=0.1)
args = parser.parse_args()

out_name = args.out_name
out_dir = args.out_dir
N = int(args.N)
N_test = int(args.N_test)
M = int(args.M)
C = int(args.C)
lam = float(args.lam)
epi_var = float(args.epi_var)
cov_var = float(args.cov_var)
prevalence = float(args.prevalence)

print("Input arguments:")
print("--out_dir", out_dir)
print("--out_name", out_name)
print("--N", N)
print("--N-test", N_test)
print("--M", M)
print("--C", C)
print("--lam", lam)
print("--prevalence", prevalence)
print("--cov-var", cov_var)
print("--epi-var", epi_var, flush=True)

print("\n...Simulating design matrix", flush=True)
X = np.random.normal(0,1,(N + N_test)*M).reshape(((N + N_test),M))

print("\n...Simulating marker effects", flush=True)
CM = int(M * lam) # number of cuasal markers
sigma2 = epi_var / CM
idx = random.sample(range(M), CM)
beta = np.zeros(M)
beta[idx] = np.random.normal(0.0, np.sqrt(sigma2), CM)
print("sigma2 = ", sigma2, flush=True)

print("\n...Computing outcome variable", flush=True)
g = np.matmul(X,beta)

y = g + np.random.normal(0, np.sqrt( 1 - np.var(g)),  (N+ N_test)) # adding Gaussian noise

print("Var(y) = ", np.var(y))
print("Mean(y) = ", np.mean(y))
print("Var(g) = ", np.var(g))

intercept = norm.ppf(prevalence)
print("Intercept = ", intercept)
ybin = ((y + intercept) >= 0).astype(int)
print("Observed prevalence = ", np.mean(ybin), flush=True)

# Saving data in binary format
print("\n...Saving train design matrix to bin file")
X_train = X[:N,:]
bin_fpath = os.path.join(out_dir, "%s.bin" % out_name)
print(bin_fpath, flush=True)
binf = open(bin_fpath, "wb")
b = struct.pack(str(N*M)+'d', *X_train.transpose().ravel().squeeze())
binf.write(b)
binf.close()

print("\n...Saving test design matrix to bin file")
X_test = X[N:,:]
bin_fpath = os.path.join(out_dir, "%s_test.bin" % out_name)
print(bin_fpath, flush=True)
binf = open(bin_fpath, "wb")
b = struct.pack(str(N_test*M)+'d', *X_test.transpose().ravel().squeeze())
binf.write(b)
binf.close()

# Saving phenotype data
print("\n...Saving outcome to .phen file")
phen_fpath = os.path.join(out_dir, "%s.phen" % out_name)
phen_fpath_bin = os.path.join(out_dir, "%s_binary.phen" % out_name)
print(phen_fpath)
print(phen_fpath_bin, flush=True)
df_phen = pd.DataFrame({"IID": range(N), "FID": range(N), "Y": y[:N], "Y_BIN": ybin[:N]})
df_phen[["IID", "FID", "Y"]].to_csv(phen_fpath, index=False, header=None, sep="\t")
df_phen[["IID", "FID", "Y_BIN"]].to_csv(phen_fpath_bin, index=False, header=None, sep="\t")

print("\n...Saving test outcome to .phen file")
phen_fpath_test = os.path.join(out_dir, "%s_test.phen" % out_name)
phen_fpath_bin_test = os.path.join(out_dir, "%s_binary_test.phen" % out_name)
print(phen_fpath_test)
print(phen_fpath_bin_test, flush=True)
df_phen_test = pd.DataFrame({"IID": range(N_test), "FID": range(N_test), "Y": y[N:], "Y_BIN": ybin[N:]})
df_phen_test[["IID", "FID", "Y"]].to_csv(phen_fpath_test, index=False, header=None, sep="\t")
df_phen_test[["IID", "FID", "Y_BIN"]].to_csv(phen_fpath_bin_test, index=False, header=None, sep="\t")

# Save true signals to file
print("\n...Saving true signals to bin file")
ts_bin_fpath = os.path.join(out_dir, "%s_ts.bin" % out_name)
print(ts_bin_fpath, flush=True)
ts_binf = open(ts_bin_fpath, "wb")
ts_binf.write(struct.pack(str(M)+'d', *beta.squeeze()))
ts_binf.close()

gri = np.zeros(M).astype(int)
df = pd.DataFrame({"GRI": gri})
df.to_csv(os.path.join(out_dir, "%s.gri" % out_name), header=None, sep=" ")

dim = "%d %d" % (N, M)
f = open(os.path.join(out_dir, "%s.dim" % out_name), "w")
f.write(dim)
f.close()

grm = "0.0000 0.0001 0.0010 0.0100 0.1000 1.0000"
f = open(os.path.join(out_dir, "%s.grm" % out_name), "w")
f.write(grm)
f.close()