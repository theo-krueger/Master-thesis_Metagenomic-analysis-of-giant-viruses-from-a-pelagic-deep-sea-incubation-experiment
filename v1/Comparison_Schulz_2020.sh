# use db from before
base_dir=$WORK/bird_bags/YuChens_pipeline/

mkdir $base_dir/Comparison_Schulz/v1

# get own polB contigs
for file in $(ls $base_dir/03_annotated/v1/NCVOGs/); do
    bin=$(echo $file | sed 's/.fasta.faa.NCVOG_tblout//g')
    awk -v bin=$bin '$0 ~ /^NCVOG0038/ {print bin, $3}' $base_dir/03_annotated/v1/NCVOGs/$file \
        >> $base_dir/Comparison_Schulz/v1/polB_contigs_YC_V1.txt
done

# get sequences
conda activate seqtk
while IFS= read line; do
    bin=$(echo $line | cut -f1 -d " ")
    contig=$(echo $line | cut -f2 -d " ")
    grep -w "$contig" $base_dir/02_predicted/v1/"$bin".fasta.faa \
        | sed 's/>//g' \
        | seqtk subseq $base_dir/02_predicted/v1/"$bin".fasta.faa - \
        | awk -v pat="$contig" -v rep="${bin}-${contig}" '{gsub(pat, rep)} 1' \
        >> $base_dir/Comparison_Schulz/v1/polB_contigs_YC_V1.faa
done < $base_dir/Comparison_Schulz/v1/polB_contigs_YC_V1.txt
awk -F " " '{if ($0 ~ /^>/) { print $1;} else { print $0}}' $base_dir/Comparison_Schulz/v1/polB_contigs_YC_V1.faa \
    > $base_dir/Comparison_Schulz/v1/polB_contigs_YC_V1_shortheaders.faa

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
        --query $base_dir/Comparison_Schulz/v1/polB_contigs_YC_V1_shortheaders.faa \
        --db $WORK/bird_bags/db/Schulz_2020/04_db/Schulz_polB.dmnd \
        --out $base_dir/Comparison_Schulz/v1/Comparison_YC_V1_polB_vs_Schulz2020.diamond \
        --threads 32 \
        --evalue 1e-04 \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
        --max-target-seqs 10 \
        "

# graphic analysis in R: Comparison_plot_Schulz2020.Rmd

