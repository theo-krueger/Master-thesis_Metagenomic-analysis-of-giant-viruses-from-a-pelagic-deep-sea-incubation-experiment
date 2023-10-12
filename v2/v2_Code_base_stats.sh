
base_dir=$WORK/bird_bags/YuChens_pipeline

### bbmap
in_dir=$WORK/bird_bags/YuChens_pipeline/01_DNA_sequences/v2
out_dir=$WORK/bird_bags/YuChens_pipeline/00_stats/v2/bbmap
for bin in $(ls $in_dir); do
  stats.sh \
	in=$in_dir/$bin \
	out=$out_dir/"$bin"_report.txt
done

for file in $(ls $base_dir/00_stats/v2/bbmap/); do
  bin=$(echo $file | cut -f1 -d '-')
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

    awk -v group=$method -v file=$file 'NR > 22 && NF {printf ("%s\t%s\t%s\n", $0, file, group)}' $base_dir/00_stats/v2/bbmap/$file >> $base_dir/00_stats/v2/bbmap_combined.table
done

### overview

conda activate bioawk

# avg gc

# for file in $(ls $base_dir/01_DNA_sequences/v2/); do
#   bin=$(echo $file | sed 's/-contigs.fa//g')
#   gc=$(bioawk -c fastx '{ print gc($seq), $name}' $base_dir/01_DNA_sequences/v2/$file | \
#   awk -F ',' '{s+=$1;} END{print s/NR}')
#   echo $bin $gc >> $base_dir/00_stats/v2/genomes_gc.txt
# done

for file in $(ls $base_dir/00_stats/v2/bbmap/); do
  bin=$(echo $file | sed 's/-contigs.fa_report.txt//g')
  gc_avg=$(head -n2 $base_dir/00_stats/v2/bbmap/$file | tail -n1 | awk '{print $8}')
  gc_stddev=$(head -n2 $base_dir/00_stats/v2/bbmap/$file | tail -n1 | awk '{print $9}')
  # careful with N50, bbmap drops trailing 0s
  N50=$(head -n9 $base_dir/00_stats/v2/bbmap/$file | tail -n1 | cut -f3 -d '/' | cut -f1 -d ' ' | sed 's/\./,/g')
  echo $bin $gc_avg $gc_stddev $N50 >> $base_dir/00_stats/v2/genomes_gc_N50.txt
done

# n contigs
for file in $(ls $base_dir/01_DNA_sequences/v2/) ; do
  bin=$(echo $file | sed 's/-contigs.fa//g')
  contig_n=$(grep ">" $base_dir/01_DNA_sequences/v2/$file | wc -l)
  echo -e "$bin \t $contig_n" >> $base_dir/00_stats/v2/genomes_n_contigs.txt
done

# calculate genome size
for file in $(ls $base_dir/01_DNA_sequences/v2/) ; do
  bin=$(echo $file | sed 's/-contigs.fa//g')
  genome_length=$(bioawk -c fastx '{print length ($seq), $name}' $base_dir/01_DNA_sequences/v2/$file | \
  awk -F',' '{sum+=$1;} END{print sum;}')
  echo -e "$bin \t $genome_length" >> $base_dir/00_stats/v2/genomes_lengths.txt
done

# n proteins
for file in $(ls $base_dir/02_predicted/v2/); do
  bin=$(echo $file | cut -f1 -d '-')
  proteins=$(grep ">" $base_dir/02_predicted/v2/$file | wc -l)
  echo $bin $proteins >> $base_dir/00_stats/v2/genomes_n_proteins.txt
done

# NCVOGs
for file in $(ls $base_dir/03_annotated/v2/NCVOGs/); do
  bin=$(echo $file | cut -f1 -d "-")
  count=$(grep '^NCVOG' $base_dir/03_annotated/v2/NCVOGs/$file | wc -l)
  count_uniq=$(grep '^NCVOG' $base_dir/03_annotated/v2/NCVOGs/$file | cut -f1 -d ' ' | sort | uniq | wc -l)

  echo $bin $count $count_uniq >> $base_dir/00_stats/v2/genomes_NCVOG_stats.txt 
done

#join
join -j 1 -o 1.1,1.2,1.3,1.4,2.2 \
  <(sort -k1 $base_dir/00_stats/v2/genomes_gc_N50.txt) \
  <(sort -k1 $base_dir/00_stats/v2/genomes_n_contigs.txt) \
  > $base_dir/00_stats/v2/gc_contigs.temp
join -j 1 -o 1.1,1.2,1.3,1.4,1.5,2.2 \
  <(sort -k1 $base_dir/00_stats/v2/gc_contigs.temp) \
  <(sort -k1 $base_dir/00_stats/v2/genomes_lengths.txt) \
  > $base_dir/00_stats/v2/gc_contigs_length.temp
join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,2.2 \
  <(sort -k1 $base_dir/00_stats/v2/gc_contigs_length.temp) \
  <(sort -k1 $base_dir/00_stats/v2/genomes_n_proteins.txt) \
  > $base_dir/00_stats/v2/gc_contigs_length_proteins.temp
join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,2.2,2.3 \
  <(sort -k1 $base_dir/00_stats/v2/gc_contigs_length_proteins.temp) \
  <(sort -k1 $base_dir/00_stats/v2/genomes_NCVOG_stats.txt) \
  > $base_dir/00_stats/v2/possible_NCLDV_summary.txt
rm $base_dir/00_stats/v2/*.temp
# add header
echo -e 'bin gc_avg gc_stddev N50 n_contigs genome_length n_proteins n_NCVOG n_NCVOG_uniq' | \
  cat - $base_dir/00_stats/v2/possible_NCLDV_summary.txt \
  > $base_dir/00_stats/v2/possible_NCLDV_summary_wtitle.txt
# visual analysis in R


# only for 32 bins after filter and NCVOG cutoff
for bin in $(cat $base_dir/00_stats/v2/v2_bins_afterfilter_above_25_uniqNCVOG.txt); do
  grep "$bin" $base_dir/00_stats/v2/possible_NCLDV_summary.txt \
    >> $base_dir/00_stats/v2/possible_NCLDV_afterfilter_summary.txt
done
sort \
  -k1 \
  -o $base_dir/00_stats/v2/possible_NCLDV_afterfilter_summary.txt \
  $base_dir/00_stats/v2/possible_NCLDV_afterfilter_summary.txt
# add header
echo -e 'bin gc_avg gc_stddev N50 n_contigs genome_length n_proteins n_NCVOG n_NCVOG_uniq' | \
  cat - $base_dir/00_stats/v2/possible_NCLDV_afterfilter_summary.txt \
  > $base_dir/00_stats/v2/possible_NCLDV_afterfilter_summary_wtitle.txt





# NCVOG table
# separate folder
mkdir $base_dir/03_annotated/v2_32/
mkdir $base_dir/03_annotated/v2_32/NCVOGs

for bin in $(cat $base_dir/00_stats/v2/v2_bins_afterfilter_above_25_uniqNCVOG.txt); do
    cp $base_dir/03_annotated/v2/NCVOGs/"$bin"-contigs.fa.faa.NCVOG_tblout $base_dir/03_annotated/v2_32/NCVOGs/
done
# download this
# run R scipt on local: Create_NCVOG_abundance_table.Rmd

## for single contigs
# contig length
for file in $(ls $base_dir/01_DNA_sequences/v2/) ; do
  bin=$(echo $file | sed 's/-contigs.fa//g')
  bioawk -v bin=$bin -t -c fastx '{print bin, $name, length($seq)}' $base_dir/01_DNA_sequences/v2/$file \
    >> $base_dir/00_stats/v2/genomes_contig_lengths.txt
done
# contig gc
for file in $(ls $base_dir/01_DNA_sequences/v2/); do
  bin=$(echo $file | sed 's/-contigs.fa//g')
  bioawk -v bin=$bin -t -c fastx '{ print bin, $name, gc($seq)}' $base_dir/01_DNA_sequences/v2/$file \
  >> $base_dir/00_stats/v2/genomes_contig_gc.txt
done

# for 32
for bin in $(cat $base_dir/00_stats/v2/v2_bins_afterfilter_above_25_uniqNCVOG.txt); do
  grep "$bin" $base_dir/00_stats/v2/genomes_contig_lengths.txt \
    >> $base_dir/00_stats/v2/genomes_contig_lengths_afterfilter_over25uniqNCVOG.txt
  grep "$bin" $base_dir/00_stats/v2/genomes_contig_gc.txt \
    >> $base_dir/00_stats/v2/genomes_contig_gc_afterfilter_over25uniqNCVOG.txt
done

### ANI
for bin in $(ls 01_DNA_sequences/v2/); do
    echo $WORK/bird_bags/YuChens_pipeline/01_DNA_sequences/v2/$bin >> 00_stats/v2/genomes_pathlist.txt
done
fastANI --ql stats/v2/genomes_pathlist.txt --rl stats/v2/genomes_pathlist.txt -o stats/v2/fastANI_allVSall_out.txt

### 10 core NCVOG

for file in $(ls 03_annotated/v2/NCVOGs/); do
    bin=$(echo $file | sed 's/-contigs.fa.faa.NCVOG_tblout//g')

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

    for NCVOG in NCVOG0023 NCVOG0037 NCVOG0038 NCVOG0076 NCVOG0249 NCVOG0262 NCVOG0271 NCVOG0274 NCVOG1117 NCVOG1192; do
        count=$(grep "$NCVOG" 03_annotated/v2/NCVOGs/$file | wc -l)
        echo $method $bin $NCVOG $count >> 00_stats/v2/core_NCVOGs/$bin-stats.txt
        echo $method $bin $NCVOG $count >> 00_stats/v2/core_NCVOGs/NCVOG_10core_summary-stats.txt
    done
done

### abundance

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
    echo 'Failed to read collection'
  fi
  
  cat 05_anvio/summaries/v1_refined/$collection/bin_by_bin/$bin/$bin-abundance.txt |\
    sed "s/abundance/$bin/g" \
    >> 00_stats/v2/abundances_bins_after_filter.txt

done < $WORK/bird_bags/YuChens_pipeline/00_stats/v2/v2_bins_remaining_afterfilter.txt
awk 'NR==1 || $0 !~ /^bin/ {print}' 00_stats/v2/abundances_bins_after_filter.txt \
  > temp && rm 00_stats/v2/abundances_bins_after_filter.txt && mv temp 00_stats/v2/abundances_bins_after_filter.txt

 