# necessary for grabbing all packages and files for annotation pipeline

# creates basic python 3.6 environment
conda create -n tempus python=3.6 anaconda
conda config --add channels bioconda
conda install -n tempus ensembl-vep

# this yml file is uploaded to github so that someone doesn't need to run build.sh, they can just create the environment from environment.yml
conda env export > environment.yml
