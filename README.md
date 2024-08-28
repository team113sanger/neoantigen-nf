# neoantigen-nf

## Step 1 - Create a conda environment
```
conda env create -n neoantigen-pipeline-env --file conda.yml
```

## Step 2 - Add the following line to the nextflow.config file
```
process.conda = 'path/to/neoantigen-pipeline-env'
```
If the line above already exists, update the path attributed to `process.conda` to reflect the location of the `neoantigen-pipeline-env` environment on your system.

## Step 3 - Run the neoantigen pipeline
```
module load nextflow
nextflow run main.nf -params-file params.json -c nextflow.config -profile farm22
```