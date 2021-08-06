# Somatic Variant & CNA Calling Pipeline for human Whole Exom/targeted NGS data

* Fully customizable pipeline for calling of somatic variants and copy number alterations 
* Runs on SGE / SLURM cluster environments

## Setup

### Prerequisites
* linux environment
* package manager conda is installed (see [official conda documentation](https://conda.io/projects/conda/en/latest/user-guide/install/linux.html) for instructions)

### Install pipeline
You need access to a cluster environment (documentation and updates only for SLURM workload management system!)

* enter cluster environment and move to folder you want to install the pipeline code into
* clone the somVarWES code from github and run the init.sh to install packages and conda environment

```
$ git clone https://github.com/Mar111tiN/somVarWES.git && cd somVarWES && source setup/init.sh
```

## Test the Pipeline
* for testing, I provide testdata as dropbox links for you to download by running simple scripts
* the downloaded bam folder contains bam files of AML tumor normal pairs in three different sizes (500MB to 25GB)
* the bam files are derived from exom-sequenced samples that have been prepped using the Agilent SureSelect XT-HS kit with HumanAllExome_v7 baits

### download bam files of test samples
* prepare a data folder with appr. 70GB of space and run the download_testbams.sh with the folder path as first argument
```
$ . setup/download_testbams.sh <path-to-data-folder>
```
* adjust \<TESTDATA\> and provide a working directory (\<WKDIR\>) in the yaml config file configs/config_test.yaml

### download static files
* for the pipeline to work for bam files created with this specific exon prep, several database files and other accessory data have to be prepared
* for testing (and for general use with SureSelect HAEv7 data), static data is provided as dropbox links
* run the script download_static.sh with the path to the desired static folder (50GB of space is required) as argument:
```
$ . setup/download_static.sh <path-to-static-folder>
```
* adjust \<STATIC\> to the desired path and provide a working directory (\<WKDIR\>) in the yaml config file configs/general/config_testgeneral.yaml

### get the annovar executables
* most tools used by the pipeline are provided as packages and are installed via the init.sh script
* However, the tool annovar used for populating mutation lists with data from multiple databases has to be downloaded [here](https://www.openbioinformatics.org/annovar/annovar_download_form.php) upon signing a user agreement
* after downloading the annovar executables, provide the path to the folder containing the annovar perl scripts (\<ANNOPATH\>) in the TOOLS section of configs/general/config_testgeneral.yaml 

### run the test
* for testing all bam files (small to large) make the config/config_test.yaml the active config which is used by the pipeline masterfile Snakefile:
```
$ cp configs/config_test.yaml configs/active_config.yaml
```
* now, you can test the pipeline either interactively or via job submission:
  + from an interactive slurm-session with >20 cores go to the somVarWES root folder and run: 
    ```
    $ conda activate snake-env
    $ snakemake -prj 20 --use-conda
    ```
  + as batch job from the somVarWES root folder:
    ```
    sbatch -J testing SLURM.sh
    ```
* either way, the pipeline will first install all conda environments into your workdir in .snakemake folder and then start the pipeline
* depending on cluster traffic, this can initially take 30 min or more!
