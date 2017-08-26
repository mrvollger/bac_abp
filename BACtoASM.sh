#!/bin/bash

# these commands load the proper env for snakemake and the pipeline
module purge
. /etc/profile.d/modules.sh
module load modules modules-init modules-gs/prod modules-eichler
module load anaconda/20161130
base2="/net/eichler/vol21/projects/bac_assembly/nobackups/scripts"


NPROC=$(nproc)
echo "number of cores: $NPROC"

if [ "$1" == "unlock" ]; then
    snakemake --unlock -s $base2/BACtoASM.py
fi

#snakemake --help
snakemake --config MINCOV=1.5 --cores $NPROC -s $base2/BACtoASM.py



