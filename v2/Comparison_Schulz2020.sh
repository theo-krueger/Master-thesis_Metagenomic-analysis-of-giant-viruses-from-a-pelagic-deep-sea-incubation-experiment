# use db from before
base_dir=$WORK/bird_bags/YuChens_pipeline/

mkdir $base_dir/Comparison_Schulz/
mkdir $base_dir/Comparison_Schulz/v2_32

# get own polB contigs
for bin in $(cat $base_dir/00_stats/v2/v2_bins_afterfilter_above_25_uniqNCVOG.txt); do
    file="$bin-contigs.fa.faa.NCVOG_tblout"
    awk -v bin=$bin '$0 ~ /^NCVOG0038/ {print bin, $3}' $base_dir/03_annotated/v2/NCVOGs/$file >> $base_dir/Comparison_Schulz/v2_32/polB_contigs_YC_32.txt
done

# get sequences
conda activate seqtk
while IFS= read line; do
    bin=$(echo $line | cut -f1 -d " ")
    contig=$(echo $line | cut -f2 -d " ")
    grep -w "$contig" $base_dir/02_predicted/v2/"$bin"-contigs.fa.faa | sed 's/>//g' | seqtk subseq $base_dir/02_predicted/v2/"$bin"-contigs.fa.faa - | awk -v pat="$contig" -v rep="${bin}-${contig}" '{gsub(pat, rep)} 1' >> $base_dir/Comparison_Schulz/v2_32/polB_contigs_YC_32.faa
done < $base_dir/Comparison_Schulz/v2_32/polB_contigs_YC_32.txt
awk -F " " '{if ($0 ~ /^>/) { print $1;} else { print $0}}' $base_dir/Comparison_Schulz/v2_32/polB_contigs_YC_32.faa \
    > $base_dir/Comparison_Schulz/v2_32/polB_contigs_YC_32_shortheaders.faa

# blast against existing db
conda activate diamond-0.9.14
sbatch \
    --job-name=diamond_Schulz_polB \
    --nodes=1 \
    --tasks-per-node=32 \
    --cpus-per-task=1 \
    --mem=120G \
    --time=01:00:00 \
    --output=$base_dir/reports/diamondSchulzpolB_out \
    --error=$base_dir/reports/diamondSchulzpolB_err \
    --mail-type=END,FAIL,TIME_LIMIT \
    --mail-user=theokrueger.infomail@gmail.com \
    --wrap="\
        diamond blastp \
        --query $base_dir/Comparison_Schulz/v2_32/polB_contigs_YC_32_shortheaders.faa \
        --db $WORK/bird_bags/db/Schulz_2020/04_db/Schulz_polB.dmnd \
        --out $base_dir/Comparison_Schulz/v2_32/Comparison_YC_32_polB_vs_Schulz2020.diamond \
        --threads 32 \
        --evalue 1e-04 \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
        --max-target-seqs 10 \
        "
# add header
echo -e 'qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstitle' | \
    cat - $base_dir/Comparison_Schulz/v2_32/Comparison_YC_32_polB_vs_Schulz2020.diamond \
    > $base_dir/Comparison_Schulz/v2_32/Comparison_YC_32_polB_vs_Schulz2020_wtitle.diamond

# graphic analysis in R: v2_Comparison_plot_Schulz2020.Rmd

# blast ALL polB from tree
conda activate diamond-0.9.14
diamond blastp \
    --query $base_dir/Phylogeny/tree07_allpolB_v1_unique/polB_contigs_YC_V1_unique_shortheaders.faa \
    --db $WORK/bird_bags/db/Schulz_2020/04_db/Schulz_polB.dmnd \
    --out $base_dir/Comparison_Schulz/v1_23/Comparison_YC_v1_23_polB_vs_Schulz2020.diamond \
    --threads 32 \
    --evalue 1e-04 \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
    --max-target-seqs 5
# add header
echo -e 'qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstitle' | \
    cat - $base_dir/Comparison_Schulz/v1_23/Comparison_YC_v1_23_polB_vs_Schulz2020.diamond \
    > $base_dir/Comparison_Schulz/v1_23/Comparison_YC_v1_23_polB_vs_Schulz2020_wtitle.diamond
# graphic analysis in R: v2_Comparison_plot_Schulz2020.Rmd
