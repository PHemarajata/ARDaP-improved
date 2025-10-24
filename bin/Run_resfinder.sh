#!/bin/bash

#variables?

baseDir=$1
forward=$2
reverse=$3
id=$4

# run resfinder inside the process workdir and write a per-sample output filename
# Set PYTHONPATH to include the resfinder directory for proper module imports
export PYTHONPATH="${baseDir}/bin/resfinder:$PYTHONPATH"
python3 ${baseDir}/bin/resfinder/run_resfinder.py -ifq ${forward} ${reverse} -acq -db_res ${baseDir}/Databases/Resfinder_general

# resfinder writes pheno_table.txt — move it to a per-sample file
if [ -s pheno_table.txt ]; then
	mv pheno_table.txt ${id}_resfinder.txt
else
	# create an empty per-sample file so downstream processes see a file
	touch ${id}_resfinder.txt
fi