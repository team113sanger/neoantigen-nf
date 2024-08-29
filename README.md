# neoantigen-nf

## Step 0 - Install Miniconda
Go to the [Miniconda website](https://docs.anaconda.com/miniconda/) and download Miniconda3. Choose the installer link that best suits your system. You can download the installer via the command line by running the command below (please note that the link used in the example below corresponds to Miniconda3 for Linux).
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

## Step 1 - Create a conda environment
The configuration file `conda.yml` used in the command below is provided in this repository.
```
conda env create -n neoantigen-pipeline-env --file conda.yml
```

## Step 2 - Add the following line to the nextflow.config file
If the line below already exists, update the path attributed to `process.conda` to reflect the location of the `neoantigen-pipeline-env` environment on your system.
```
process.conda = 'path/to/neoantigen-pipeline-env'
```

## Step 3 - Enable conda to be used
Depending on how you installed conda on your system, it may be that conda is not automatically initialised when you get into a terminal and therefore running `conda activate` or any other conda-related commands will not work. If that is the case, run the following command on your command line before triggering the `neoantigen-nf` pipeline.
```
CONDA_PATH=path/to/miniconda3/bin/conda
eval "$(${CONDA_PATH} shell.bash hook)"
```

## Step 4 - Run the neoantigen pipeline
```
module load nextflow
nextflow run main.nf -params-file params.json -c nextflow.config -profile farm22
```