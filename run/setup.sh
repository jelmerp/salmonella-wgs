## Plasmidfinder
mamba create -p /fs/project/PAS0471/jelmer/conda/plasmidfinder-2.1.6 -c bioconda plasmidfinder=2.1.6
#Please run download-db.sh to download the PlasmidFinder database to /fs/project/PAS0471/jelmer/conda/plasmidfinder-2.1.6/share/plasmidfinder-2.1.6/database.
#If you have a database in custom path, please use plasmidfinder.py with the option -p.
download-db.sh

## Get the Busco DB for QUAST
wget http://busco.ezlab.org/v2/datasets/bacteria_odb9.tar.gz
mv bacteria_odb9.tar.gz /fs/project/PAS0471/jelmer/conda/quast-5.0.2/lib/python3.7/site-packages/quast_libs/busco/bacteria.tar.gz

## Subset data for testing
bash mcic-scripts/misc/fqsub_dir.sh -i data/fastq -o data/subset_fastq

## Install software
mamba create -y -p /fs/project/PAS0471/jelmer/conda/abricate-1.0.1 -c bioconda abricate=1.0.1
mamba create -y -p /fs/project/PAS0471/jelmer/conda/knsp-3.1 -c hcc ksnp=3.1
mamba create -p /fs/project/PAS0471/jelmer/conda/prokka-1.14.6 -c conda-forge -c bioconda -c defaults prokka=1.14.6
mamba create -p /fs/project/PAS0471/jelmer/conda/amrfinderplus-3.10.30 -c bioconda ncbi-amrfinderplus=3.10.30
mamba create -y -n treetime-0.9.4 -c bioconda treetime=0.9.4 fasttree iqtree
mamba create -y -n snippy-4.6.0 -c bioconda snippy=4.6.0
mamba create -y -n gubbins-3.2.1 -c bioconda gubbins=3.2.1
