# GATK Best Practices - Somatic Variant Calling for human Whole Exom Data

trim > map(BWA) > dedup(picard) > indelRealign(gatk) > baseRecalibrate(gatk) > variantCalling(varscan/muTect2) > annotation(annovar) > filter(custom)
## Setup

### Setup local
* copy files from setup/local into your local setup 
* copy stuff into your own .bash*
* adjust for OS dialect
* produce your .ssh folder following the BIH guideline
enter the BIH server directly to computation node

```
$ ssh bihcluster qrsh
```
* provide password
* ready
* copy .bash_profile .bashrc .condarc to your $HOME

### SETUP Conda
* install miniconda:
```
$ wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
```
* set to your path
```
$ bash Miniconda3-latest-Linux-x86_64.sh -b -f -p $HOME/work/miniconda
```
* add $HOME/work/miniconda/bin to PATH:
```
$ 'export PATH=$PATH:$HOME/work/miniconda/bin' >> ~/.bashrc
```
* refresh bash:
```
$ source ~/.bashrc
```
* activate conda base environment:
```
$ . activate base # or conda activate base
```

### Create Runtime Environment
* create or move to project folder <folder/> for cloning pipeline repository into
* git clone this repo into that folder:
```
$ git clone https://github.com/Mar111tiN/somVarGATK.git
```
* use env-file (hwes.yml) to create new environment (rename in first line to <your_env>)
```
$ conda env create -f setup/hwes.yml
$ . activate <your_env>
$ . deactivate
```
* env is stored centrally in ../miniconda/envs/<your_env>
* executables are stored/symlinked in ../miniconda/envs/<your_env>/bin

### INSTALL execs to your ENV
!! while in that environment!!!:
* add missing conda channels to ~.condarc file
* install what you want using conda install ...
* [ place symlinks to your files into .../miniconda/envs/<your_env> ]
    ( = $CONDA_PREFIX/envs/<your_env>)
* 
* add shell-scripts to PATH:
```
export $PATH:<folder>/scripts
```
* ready to go

## Configure Pipeline Configuration
* a sample config_blueprint.yaml is provided to guide you through the configuration of the tools
* make sure all reference and annotation files are existing (snakemake will tell you if not)




