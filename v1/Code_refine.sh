# Importing bins to anvi’o
# https://merenlab.org/2016/06/22/anvio-tutorial-v2/#anvi-import-collection
# Build a binning results file
mkdir anvio
# Gets the contig names
grep '>' 01_DNA_sequences/*fasta | sed 's/>//g' | cut -f2 -d ':' > anvio/name_contig.list
#this gets the bin name
grep '>' 01_DNA_sequences/*fasta | cut -f1 -d ':' | sed 's/.fasta//' | cut -f2 -d '/' > anvio/name_bin.list
#this pastes them together
paste anvio/name_contig.list anvio/name_bin.list | sed 's/\./_/g' > anvio/bin_contigs_collection_select_NCLDV
#make different collections
while read line; do
    bin=$(echo $line | cut -f2 -d $' ')
  if [[ "$bin" == *"coAssem_MaxBin"* ]]; then
    echo "$line" >> anvio/collections/collection_select_NCLDV_coAssem_MaxBin
  elif [[ "$bin" == *"coAssem_metabat2"* ]]; then
    echo "$line" >> anvio/collections/collection_select_NCLDV_coAssem_metabat2
  elif [[ "$bin" == *"coAssem_output_concoct"* ]]; then
    echo "$line" >> anvio/collections/collection_select_NCLDV_coAssem_concoct
  elif [[ "$bin" == *"ESOM"* ]]; then
    echo "$line" >> anvio/collections/collection_select_NCLDV_ESOM
  elif [[ "$bin" == *"concoct"* ]]; then
    echo "$line" >> anvio/collections/collection_select_NCLDV_concoct
  elif [[ "$bin" == *"metabat2"* ]]; then
    echo "$line" >> anvio/collections/collection_select_NCLDV_metabat2
  elif [[ "$bin" == *"MaxBin"* ]]; then
    echo "$line" >> anvio/collections/collection_select_NCLDV_MaxBin
  else
    echo 'Failure'
  fi
done < anvio/bin_contigs_collection_select_NCLDV
# metabat2 issue
while read line; do
  bin=$(echo $line | cut -f2 -d $' ')
  echo "$line" >> anvio/collections/metabat2/collection_select_NCLDV_metabat2_$bin
done < anvio/collections/collection_select_NCLDV_metabat2

# import “bin_contigs_collection_original” to anvi’o
for collection in $(ls anvio/collections/); do
  anvi-import-collection \
    -c anvio/rename.5k.all.spades.contigs.fasta.db \
    -p anvio/merged_samples/PROFILE.db \
    -C $collection \
    --contigs-mode anvio/collections/$collection
done
# importing metabat2 extra
for collection in $(ls anvio/collections/metabat2/); do
  anvi-import-collection \
    -c anvio/rename.5k.all.spades.contigs.fasta.db \
    -p anvio/merged_samples/PROFILE.db \
    -C $collection \
    --contigs-mode anvio/collections/metabat2/$collection
done

############ export collection ######################
# start of refinement loop
# summarize abundances and coverage across different samples
mkdir anvio/summaries

for collection in $(ls anvio/collections/c* | cut -f3 -d '/'); do
  anvi-summarize \
      -c $WORK/bird_bags/Comparison_YuChens_pipeline/anvio/rename.5k.all.spades.contigs.fasta.db \
      -p $WORK/bird_bags/Comparison_YuChens_pipeline/anvio/merged_samples/PROFILE.db \
      -C $collection \
      -o $WORK/bird_bags/Comparison_YuChens_pipeline/anvio/summaries/v1_refined/$collection-output_v1_refined
done

# for metabat2
for collection in $(ls anvio/collections/metabat2/c* | cut -f4 -d '/'); do
  anvi-summarize \
      -c $WORK/bird_bags/Comparison_YuChens_pipeline/anvio/rename.5k.all.spades.contigs.fasta.db \
      -p $WORK/bird_bags/Comparison_YuChens_pipeline/anvio/merged_samples/PROFILE.db \
      -C $collection \
      -o $WORK/bird_bags/Comparison_YuChens_pipeline/anvio/summaries/v1_refined/metabat2/$collection-output_v1_refined
done

# add abundances files together
for collection in $(ls anvio/collections/collection* | cut -f3 -d '/'); do
  cat anvio/summaries/v1_refined/$collection-output_v1_refined/bins_across_samples/abundance.txt \
    >> anvio/summaries/abundances/summary_abundances_v2.txt
done
for collection in $(ls anvio/collections/metabat2/collection* | cut -f4 -d '/'); do
  cat anvio/summaries/v1_refined/metabat2/$collection-output_v1_refined/bins_across_samples/abundance.txt \
    >> anvio/summaries/abundances/summary_abundances_v2.txt
done
awk 'NR==1 || $0 !~ /^bins/ {print}' anvio/summaries/abundances/summary_abundances_v2.txt \
  >> anvio/summaries/abundances/summary_abundances_v2_clean.txt


############ anvi-refine ######################
# ssh -L 8092:localhost:8092 -X smomw539@nesh-fe.rz.uni-kiel.de
# ssh -L 8093:localhost:8093 -X smomw539@nesh-fe.rz.uni-kiel.de

# working with bins from Bins_abundances_v1
anvio_dir=$WORK/bird_bags/YuChens_pipeline/anvio

anvi-refine \
    -c $anvio_dir/rename.5k.all.spades.contigs.fasta.db \
    -p $anvio_dir/merged_samples/PROFILE.db \
    -C collection_select_NCLDV_MaxBin \
    --server-only -P 8092 \
    -b BB_2016_1040m_12C_4m_20_L002_MaxBin_bins_013_1


########### export again ######################
for collection in $(ls anvio/collections/c* | cut -f3 -d '/'); do
anvi-summarize \
      -c $WORK/bird_bags/Comparison_YuChens_pipeline/anvio/rename.5k.all.spades.contigs.fasta.db \
      -p $WORK/bird_bags/Comparison_YuChens_pipeline/anvio/merged_samples/PROFILE.db \
      -C $collection \
      -o $WORK/bird_bags/Comparison_YuChens_pipeline/anvio/summaries/v1_refined/$collection-output_v1_refined
done

# for metabat2
for collection in $(ls anvio/collections/metabat2/c* | cut -f4 -d '/'); do
  anvi-summarize \
      -c $WORK/bird_bags/Comparison_YuChens_pipeline/anvio/rename.5k.all.spades.contigs.fasta.db \
      -p $WORK/bird_bags/Comparison_YuChens_pipeline/anvio/merged_samples/PROFILE.db \
      -C $collection \
      -o $WORK/bird_bags/Comparison_YuChens_pipeline/anvio/summaries/v1_refined/metabat2/$collection-output_v1_refined
done

# add abundances files together
for collection in $(ls anvio/collections/collection* | cut -f3 -d '/'); do
  cat anvio/summaries/v1_refined/$collection-output_v1_refined/bins_across_samples/abundance.txt \
    >> anvio/summaries/abundances/summary_abundances_v2.txt
done
for collection in $(ls anvio/collections/metabat2/collection* | cut -f4 -d '/'); do
  cat anvio/summaries/v1_refined/metabat2/$collection-output_v1_refined/bins_across_samples/abundance.txt \
    >> anvio/summaries/abundances/summary_abundances_v2.txt
done
awk 'NR==1 || $0 !~ /^bins/ {print}' anvio/summaries/abundances/summary_abundances_v2.txt \
  >> anvio/summaries/abundances/summary_abundances_v2_clean.txt

############ virsort to v2 ##################

# copying input files
base_dir=$WORK/bird_bags/Comparison_YuChens_pipeline/anvio/summaries/v1_refined
out_dir=$WORK/bird_bags/Comparison_YuChens_pipeline/virsorter2/v2/01_DNA_sequences

while read bin; do
  if [[ "$bin" == *"coAssem_MaxBin"* ]]; then
    collection=collection_select_NCLDV_coAssem_MaxBin-output_v1_refined
  elif [[ "$bin" == *"coAssem_metabat2"* ]]; then
    collection=collection_select_NCLDV_coAssem_metabat2-output_v1_refined
  elif [[ "$bin" == *"coAssem_output_concoct"* ]]; then
    collection=collection_select_NCLDV_coAssem_concoct-output_v1_refined
  elif [[ "$bin" == *"ESOM"* ]]; then
    collection=collection_select_NCLDV_ESOM-output_v1_refined
  elif [[ "$bin" == *"concoct"* ]]; then
    collection=collection_select_NCLDV_concoct-output_v1_refined
  elif [[ "$bin" == *"metabat2"* ]]; then
    intermediate=$(echo $bin | cut -f1,2,3,4,5,6,7,8,9,10 -d "_")
    collection=metabat2/collection_select_NCLDV_metabat2_"$intermediate"-output_v1_refined
  elif [[ "$bin" == *"MaxBin"* ]]; then
    collection=collection_select_NCLDV_MaxBin-output_v1_refined
  else
    echo 'Failure'
  fi

  full_in_dir=$base_dir/$collection/bin_by_bin/$bin

  cp $full_in_dir/$bin-contigs.fa $out_dir/$bin-contigs.fa

done < $WORK/bird_bags/Comparison_YuChens_pipeline/anvio/summaries/abundances/names_forv2.txt

# running virsorter2

#!/bin/sh

#SBATCH --job-name=vs2
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --qos=long
#SBATCH --time=5-00:00:00
#SBATCH --output=reports/vs2_out
#SBATCH --error=reports/vs2_err
#SBATCH --partition=cluster
#SBATCH --mail-type=END,FAIL,TIME-LIMIT
#SBATCH --mail-user=theokrueger.infomail@gmail.com

source $HOME/miniconda3/bin/activate virsorter2

in_dir=$WORK/bird_bags/Comparison_YuChens_pipeline/virsorter2/v2/01_DNA_sequences
out_dir=$WORK/bird_bags/Comparison_YuChens_pipeline/virsorter2/v2/02_output

for bin in $(cat $WORK/bird_bags/Comparison_YuChens_pipeline/anvio/summaries/abundances/names_forv2.txt); do
  virsorter run \
  -w $out_dir/"$bin".vs2.out \
  -i $in_dir/$bin-contigs.fa \
  --include-groups 'dsDNAphage, ssDNA, NCLDV, RNA, lavidaviridae' \
  --min-score 0.5 all
done

# join outputs together
echo -e 'seqname\tdsDNAphage\tssDNA\tNCLDV\tRNA\tlavidaviridae\tmax_score\tmax_score_group\tlength\thallmark\tviral\tcellular\tbinname' \
  > $base_dir/06_virsorter2/v2/summary.table

for bin in $(ls $base_dir/06_virsorter2/v2/02_output/); do
  awk -v binname=$(echo $bin | cut -f1 -d '-') '{print $0, binname}' $base_dir/06_virsorter2/v2/02_output/$bin/final-viral-score.tsv \
  >> $base_dir/06_virsorter2/v2/summary.table
done

awk '{
  if (NR==1) 
    {
      print $0
    } 
  else if ($0 !~ /^seqname/ && $0 !~ /^dsDNAphage/) 
    {
      print $0
    }
  }' $base_dir/06_virsorter2/v2/summary.table \
> $base_dir/06_virsorter2/v2/virsort2_v1tov2_summary.table

while read line; do
  bin=$(echo $line | cut -f13 -d $'\t')
  if [[ "$bin" == *"coAssem_MaxBin"* ]]; then
    method="coAssem_MaxBin"
  elif [[ "$bin" == *"coAssem_metabat2"* ]]; then
    method="coAssem_metabat2"
  elif [[ "$bin" == *"coAssem_output_concoct"* ]]; then
    method="coAssem_concoct"
  elif [[ "$bin" == *"ESOM"* ]]; then
    method="ESOM"
  elif [[ "$bin" == *"concoct"* ]]; then
    method="concoct"
  elif [[ "$bin" == *"metabat2"* ]]; then
    method="metabat2"
  elif [[ "$bin" == *"MaxBin"* ]]; then
    method="MaxBin"
  else
    method="Method"
  fi

  echo -e "$line\t$method" >> $base_dir/06_virsorter2/v2/virsort2_v1tov2_summary_clean.table
done < $base_dir/06_virsorter2/v2/virsort2_v1tov2_summary.table

# continue analysis in R


###### going back
# refinement only done for increasing, this is finding bins above 25 uniq NCVOG that are increasing or decreasing
# refinement on these (find list of all above 25 before refining in R: Methods_find_possible_NCLDV_v1.Rmd)
