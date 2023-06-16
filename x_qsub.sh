#! /bin/bash

PBS_path=${1}
echo ${PBS_path}
PBS_files=`ls ${PBS_path}/*.pbs 2>/dev/null`
for Sub_file in ${PBS_files}; do
    qsub $Sub_file
done
