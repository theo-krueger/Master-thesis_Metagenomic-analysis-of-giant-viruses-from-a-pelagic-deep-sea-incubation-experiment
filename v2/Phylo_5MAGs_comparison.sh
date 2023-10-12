#!/bin/sh

cp $WORK/bird_bags/YuChens_pipeline/Phylogeny/polB_contigs_YC_ESOM_shortheaders.faa \
    $WORK/bird_bags/YuChens_pipeline/Phylogeny/comparison/v2_5MAGs

names=$(grep ">" $WORK/bird_bags/YuChens_pipeline/Phylogeny/comparison/v2_5MAGs/polB_contigs_YC_ESOM_shortheaders.faa | \
    sed 's/>//g' | \
    cut -f2 -d "-")
for name in $names; do
    echo $name
    grep "$name" $WORK/bird_bags/YuChens_pipeline/Comparison_Schulz/v1_23/Comparison_YC_v1_23_polB_vs_Schulz2020_wtitle.diamond \
        >> $WORK/bird_bags/YuChens_pipeline/Phylogeny/comparison/v2_5MAGs/comparison_Schulz.diamond
    grep -h "$name" $WORK/bird_bags/YuChens_pipeline/03_annotated/v2_polB_bins/NR/* \
        >> $WORK/bird_bags/YuChens_pipeline/Phylogeny/comparison/v2_5MAGs/comparison_NCBInr.diamond
done