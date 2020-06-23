#! /bin/usr/env bash

VV_REF=/beevol/home/walshs/analysis/references/vaccinia_mm10_ref
GTF=$VV_REF/NC_006998.gtf
REF=$HOME/Projects/antigen-tracking/ref
GENES=$REF/VACV_genes.txt

cat $GTF \
    | grep -Eo "gene_name \"[a-zA-Z0-9]+\"" \
    | grep -Eo "\"[a-zA-Z0-9]+\"$" \
    | awk '{gsub("\"", "", $1); print}' \
    > $GENES


