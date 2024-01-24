# gmrm-omi
Bayesian grouped mixture of regressions model for analysing omics data. This software is based on recently developed GMRM for GWAS. More information:
https://github.com/medical-genomics-group/gmrm/wiki

# Compilation

Cloning from github
```
git clone https://github.com/medical-genomics-group/gmrm-omi.git
cd gmrm-omi
```

For gcc and openmpi:

```
cd setup/
cp Make.gcc_mvapich2 Make.gcc_openmpi
cd ..
```

Load modules:

```
module load gcc openmpi
sh configure gcc_openmpi
```

This will create a Makefile that should compile the software

```
make
```

This will create a new folder, which in this instance is called build_gcc_openmpi and within that folder should be a compiled binary.

# Input options

| Option | Description |
| --- | --- |
| `--model` | Regression model `linear` / `probit` . Defualt is `linear`|
| `--bin-file` | Path to .bin file including the data in binary format in double precision |
| `--phen-files` | Comma separated (no space!) list of phenotype files to be processed. At least one phenotype file is expected. |
| `--cov-file` | File including covariates to be processed (for probit model). |
| `--cov-num` | Number of covariates |
| `--dim-file` | Path to the file containing the dimensions of the genotype: expected to be a single line file containing 2 integers: N and M, N being the number of individuals and M the number of markers. |
| `--group-index-file` | Path to group index file. One line for each marker. On each line a label and the group to which the markers belongs to. Each marker must belong to a group. Groups start at 0 and are successive. Like 0,1,2,3 if you have 4 groups. At least a single group with all markers. |
| `--group-mixture-file` | Path to group mixture file. For each group, a line describing the component mixture. Must start with 0.0. Next entry must be strictly greater. For all groups, the same number of components is required. |
| `--verbosity` | Default is 0, printing the bare minimum information on the processing. Level 3 is the debugging mode, printing extensively. |
| `--shuffle-markers` | Shuffling (1) or not (0) the markers. Default is shuffling. |
| `--seed` | Seed for pseudo-random number generator (Mersenne twister). |
| `--iterations` | Number of iterations to run. |
| `--trunc-markers` | Truncate the number of markers to process. |
| `--predict` | When specified, running prediction and produce .yest file. |
| `--out-dir` | Output directory, where to store output files (.bet, .csv, _cov.csv) |
| `--in-name-base` | Used when predicting. Estimates of marker and covariate effects will be loaded from *.bet and *_cov.csv files. |

# Example
Assuming standard HPC environment, we can run inferece for linear model:
```
module load gcc openmpi boost

srun path/to/compiled/binary/gmrm \
    --bin-file path/to/file/example.bin \
    --dim-file path/to/file/example.dim \
    --phen-files path/to/file/example.phen \
    --group-index-file path/to/file/example.gri \
    --group-mixture-file path/to/file/example.grm \
    --shuffle-markers 1 \
    --iterations 1000 \
    --out-dir output/directory/
```