# neoantigen-nf

## Step 1 - Create a conda environment
```
conda env create -n neoantigen-nf --file conda.yml
```

## Step 2 - Run the neoantigen pipeline
```
module load nextflow
nextflow run main.nf -params-file params.json -c nextflow.config -profile farm22
```