# Example
Example files are provided here to show how to run the software. 

## Data simulation
Python script ```data_sim.py``` can be used to generate example i.i.d. design matrix (example.bin file), marker effects (example_ts.bin file) and coresponding continues outcome (example.phen file). Other files (example.dim, example.grm, example.gri) are already included here, suited for default ```data_sim.py``` setup. 

```
python3 data_sim.py --out-dir ...
```

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

- .csv file with the model's hyper-parameters (variance explained by group, residual variance, SNP-heritability and markers in the model per sample). Each column of the file represents a hyperparameter and each row a posterior sample.

- .bet binary file, which stores the effect sizes per SNP per posterior sample. First 4 bytes INT specifes number of markers M. Than follows iteration number (4 byte INT) and posterior samples of marker effects in current iteration M * (8 byte DOUBLE).
