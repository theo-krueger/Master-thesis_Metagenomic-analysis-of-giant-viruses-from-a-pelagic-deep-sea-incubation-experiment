#!/bin/sh

base_dir=$WORK/bird_bags/YuChens_pipeline

# get number of proteins
for file in $(ls $base_dir/02_predicted/v2/); do
  bin=$(echo $file | cut -f1 -d '-')
  proteins=$(grep ">" $base_dir/02_predicted/v2/$file | wc -l)
  echo $bin $proteins >> $base_dir/00_stats/v2/genomes_n_proteins.txt
done

# n annotated eggnog
for file in $(ls $base_dir/03_annotated/v2/eggnog/annotations/clean); do
    bin=$(echo $file | cut -f1 -d '-')
    proteins=$(cat $base_dir/03_annotated/v2/eggnog/annotations/clean/$file | wc -l )
    echo $bin "$(($proteins-1))" >> $base_dir/00_stats/v2/genomes_n_annotated_eggnog.txt
done

# n annotated NCBI
# local R: annotation_stats.Rmd

# n pdb
for file in $(ls $WORK/bird_bags/YuChens_pipeline/03_annotated/v2_polB_bins/pdb/); do
  awk '{ a[$1]++ } END { for (b in a) { print b } }' $WORK/bird_bags/YuChens_pipeline/03_annotated/v2_polB_bins/pdb/$file | wc -l
done

# n swissprot
for file in $(ls $WORK/bird_bags/YuChens_pipeline/03_annotated/v2_polB_bins/sprot/); do
  awk '{ a[$1]++ } END { for (b in a) { print b } }' $WORK/bird_bags/YuChens_pipeline/03_annotated/v2_polB_bins/sprot/$file | wc -l
done
