#! /usr/bin/env bash

# Need to remove JH_17[78]_ prefix from fastq names when creating symlinks

fq_dir=$HOME/../rbilab/data/tamburini/200109_A00405_0194_AHNM2VDSXX

for file in $fq_dir/*.fastq.gz
do
    name=$(basename $file | sed 's/^JH_17[78]_//g')

    ln -s $file $name
done

