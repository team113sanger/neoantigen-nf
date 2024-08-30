# Identifying neoantigens with the `neoantigen-nf` pipeline!

The `neoantigen-nf` pipeline leverages the R-based package [Neoantimon](https://academic.oup.com/bioinformatics/article/36/18/4813/5906504) to identify neoantigens. Please follow the steps below to get started!

## Step 0 - Install Miniconda
Go to the [Miniconda website](https://docs.anaconda.com/miniconda/) and download Miniconda3. Choose the installer link that best suits your system. You can download the installer via the command line by running the command below (please note that the link used in the example below corresponds to Miniconda3 for Linux).
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

## Step 1 - Create a conda environment
The configuration file `conda.yml` used in the command below is provided in this repository.
```
conda env create -n neoantigen-pipeline-env --file conda.yml
```

## Step 2 - Add the following line to the `nextflow.config` file
If the line below already exists, update the path attributed to `process.conda` to reflect the location of the `neoantigen-pipeline-env` environment on your system.
```
process.conda = 'path/to/neoantigen-pipeline-env'
```
This will allow the `neoantigen-nf` pipeline to use the conda environment created in **Step 1**.

## Step 3 - Enable conda usage
Depending on how conda is set up on your system, it may be that conda is not automatically initialised when you get into a terminal and therefore running `conda activate` or any other conda-related commands will not work. If that is the case, run the following command on your command line before triggering the `neoantigen-nf` pipeline.
```
CONDA_PATH=path/to/miniconda3/bin/conda
eval "$(${CONDA_PATH} shell.bash hook)"
```
Please note that you will have to run this command every time you start a new terminal/shell.

## Step 4 - Download the NetMHCpan software
Follow the steps described in the [Neoantimon repository](https://github.com/hase62/Neoantimon/tree/master) to download and set up the NetMHCpan software. You can place NetMHCpan in the location that best suits you. After following that tutorial, you should have something like this:
```
netMHCpan-4.1
├── data
├── data.tar.gz
├── Linux_x86_64
├── netMHCpan
├── netMHCpan.1
├── netMHCpan-4.1.readme
├── test
└── tmp
```
Take note of the path `path/to/netMHCpan-4.1/netMHCpan` (or whatever it looks like for you on your own machine!), as you will need this path in **Step 6**.

## Step 5 - Create a file with the paths to your data
You need to create a file - say, `data_paths.csv` - containing the IDs of your samples and the paths to your VCF, BAM and BAI files. Follow the example below:
```
sample,vcf,bam,bai
sample1,path/to/vcf1,path/to/bam1,path/to/bai1
sample2,path/to/vcf2,path/to/bam2,path/to/bai2
sample3,path/to/vcf3,path/to/bam3,path/to/bai3
```
The path `path/to/data_paths.csv` will be used later on, in **Step 6**.

## Step 6 - Modify the file `params.json`
Another thing you have to do before triggering the `neoantigen-nf` pipeline is tweak the file `params.json`, so that it contains the paths to the inputs the pipeline will need. That file should be available as part of this repository, but you can also create it from scratch following the example below:
```
{
    "data_files": "path/to/data_paths.csv",
    "bed_file": "path/to/HLA.bed",
    "net_mhc_pan": "path/to/netMHCpan-4.1/netMHCpan",
    "outdir": "path/to/outdir",
    "copy_mode": "symlink"
}
```

## Step 7 - Run the `neoantigen-nf` pipeline
```
module load nextflow
nextflow run main.nf -params-file params.json -c nextflow.config -profile farm22
```