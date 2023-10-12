base_dir=$WORK/bird_bags/YuChens_pipeline

# change names
cd 01_DNA_sequences
for bin in $(ls *fasta); do
  mv $bin BB_$bin
done
cd ..


### bbmap
mkdir $base_dir/00_stats/v1/bbmap/
conda activate bbmap
for bin in $(ls $base_dir/01_DNA_sequences/v1_raw/$bin); do
  stats.sh \
	in=$base_dir/01_DNA_sequences/v1_raw/$bin \
	out=$base_dir/00_stats/v1/bbmap/"$bin"_report.txt
done

for file in $(ls $base_dir/00_stats/v1/bbmap/); do
  bin=$(echo $file | sed 's/.fasta_report.txt//g')
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

    awk -v group=$method -v bin=$bin 'NR > 22 && NF {printf ("%s\t%s\t%s\n", $0, bin, group)}' $base_dir/00_stats/v1/bbmap/$file >> $base_dir/00_stats/v1/bbmap_combined.table
done


### overview
conda activate bioawk

# avg gc
for file in $(ls $base_dir/00_stats/v1/bbmap/); do
  bin=$(echo $file | sed 's/.fasta_report.txt//g')
  gc_avg=$(head -n2 $base_dir/00_stats/v1/bbmap/$file | tail -n1 | awk '{print $8}')
  gc_stddev=$(head -n2 $base_dir/00_stats/v1/bbmap/$file | tail -n1 | awk '{print $9}')
  # careful with N50, bbmap drops trailing 0s
  N50=$(head -n9 $base_dir/00_stats/v1/bbmap/$file | tail -n1 | cut -f3 -d '/' | cut -f1 -d ' ' | sed 's/\./,/g')
  echo $bin $gc_avg $gc_stddev $N50 >> $base_dir/00_stats/v1/genomes_gc_N50.txt
done

# n contigs
for file in $(ls $base_dir/01_DNA_sequences/v1_raw/) ; do
  bin=$(echo $file | sed 's/.fasta//g')
  contig_n=$(grep ">" $base_dir/01_DNA_sequences/v1_raw/$file | wc -l)
  echo -e "$bin \t $contig_n" >> $base_dir/00_stats/v1/genomes_n_contigs.txt
done

# calculate genome size
for file in $(ls $base_dir/01_DNA_sequences/v1_raw/) ; do
  bin=$(echo $file | sed 's/.fasta//g')
  genome_length=$(bioawk -c fastx '{print length ($seq), $name}' $base_dir/01_DNA_sequences/v1_raw/$file | \
  awk -F',' '{sum+=$1;} END{print sum;}')
  echo -e "$bin \t $genome_length" >> $base_dir/00_stats/v1/genomes_lengths.txt
done

for file in $(ls $base_dir/02_predicted/v1/); do
  bin=$(echo $file | sed 's/.fasta.faa//g')
  proteins=$(grep ">" $base_dir/02_predicted/v1/$file | wc -l)
  echo $bin $proteins >> $base_dir/00_stats/v1/genomes_n_proteins.txt
done

# NCVOGs
for file in $(ls $base_dir/03_annotated/v1/NCVOGs/); do
  bin=$(echo $file | sed 's/.fasta.faa.NCVOG_tblout//g')
  count=$(grep '^NCVOG' $base_dir/03_annotated/v1/NCVOGs/$file | wc -l)
  count_uniq=$(grep '^NCVOG' $base_dir/03_annotated/v1/NCVOGs/$file | cut -f1 -d ' ' | sort | uniq | wc -l)

  echo $bin $count $count_uniq >> $base_dir/00_stats/v1/genomes_NCVOG_stats.txt 
done

#join
join -j 1 -o 1.1,1.2,1.3,1.4,2.2 \
  <(sort -k1 $base_dir/00_stats/v1/genomes_gc_N50.txt) \
  <(sort -k1 $base_dir/00_stats/v1/genomes_n_contigs.txt) \
  > $base_dir/00_stats/v1/gc_contigs.temp
join -j 1 -o 1.1,1.2,1.3,1.4,1.5,2.2 \
  <(sort -k1 $base_dir/00_stats/v1/gc_contigs.temp) \
  <(sort -k1 $base_dir/00_stats/v1/genomes_lengths.txt) \
  > $base_dir/00_stats/v1/gc_contigs_length.temp
join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,2.2 \
  <(sort -k1 $base_dir/00_stats/v1/gc_contigs_length.temp) \
  <(sort -k1 $base_dir/00_stats/v1/genomes_n_proteins.txt) \
  > $base_dir/00_stats/v1/gc_contigs_length_proteins.temp
join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,2.2,2.3 \
  <(sort -k1 $base_dir/00_stats/v1/gc_contigs_length_proteins.temp) \
  <(sort -k1 $base_dir/00_stats/v1/genomes_NCVOG_stats.txt) \
  > $base_dir/00_stats/v1/possible_NCLDV_summary.txt
rm $base_dir/00_stats/v1/*.temp
# add header
echo -e 'bin gc_avg gc_stddev N50 n_contigs genome_length n_proteins n_NCVOG n_NCVOG_uniq' | \
  cat - $base_dir/00_stats/v1/possible_NCLDV_summary.txt \
  > $base_dir/00_stats/v1/possible_NCLDV_summary_wtitle.txt
# visual analysis in R



## for single contigs
# contig length
for file in $(ls $base_dir/01_DNA_sequences/v1_raw/*.fasta) ; do
  bin=$(echo $file | cut -f9)
  echo $bin
  bioawk -v bin=$bin -t -c fastx '{print bin, $name, length($seq)}' $file \
    >> $base_dir/00_stats/v1/genomes_contig_lengths.txt
done

### ANI

find 01_DNA_sequences/ -maxdepth 1 -type f -not -path '*/\.*' > stats/genomes_pathlist.txt
fastANI --ql stats/genomes_pathlist.txt --rl stats/genomes_pathlist.txt -o stats/fastANI_allVSall_out.txt

