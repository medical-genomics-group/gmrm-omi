# Example
Example files are provided here to show how to run the software. 

## Data simulation
Python script ```data_sim.py``` can be used to generate example i.i.d. design matrix (example.bin file), marker effects (example_ts.bin file) and coresponding continues outcome (example.phen file). Other files (example.dim, example.grm, example.gri) are already included here, suited for default ```data_sim.py``` setup. 

```
python3 data_sim.py --out-dir ...
```

This creates following files in specified output directory:
 - ``.bin`` file with data for each sample. Data are stored as 8 byte DOUBLE precision numbers. Whole file is sequence of M marker blocks, where each block contains N 8 byte DOUBLE precision values corespoding to each sample.
 - ``_ts.bin`` file with true signals in binary format. File contains M 8 byte DOUBLE precision values, coresponding to true signals used for simulation.
 - ``.phen`` file with simulated phenotype in PLINK format.

Other files required for running software:
 - ``.grm`` group mixture file. For each group, a line describing the prior
mixture variances. Must start with 0.0. Next entry must be strictly greater. For all groups, the same number of components is required.
 - ``.gri`` group index file. One line for each marker. On each line a label
and the group to which the markers belongs to. Each marker must
belong to a group. Groups start at 0 and are successive. Like 0,1,2,3 if
you have 4 groups. At least a single group with all markers.
 - ``.dim`` file is expected to be a single line file containing 2 integers: N and M, N being the number of individuals and M the number of markers.

## Running GMRM-omi
Assuming standard HPC environment, we can run inference for linear model as follows. Paths were example data is stored has to be specified.
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

In output directory we can than find following output files:

- ``.csv`` file with the model's hyper-parameters (variance explained by group, residual variance, SNP-heritability and markers in the model per sample). Each column of the file represents a hyperparameter and each row a posterior sample.

- ``.bet`` binary file, which stores the effect sizes per SNP per posterior sample. First 4 bytes INT specifes number of markers M. Than follows iteration number (4 byte INT) and posterior samples of marker effects in current iteration M * (8 byte DOUBLE).
