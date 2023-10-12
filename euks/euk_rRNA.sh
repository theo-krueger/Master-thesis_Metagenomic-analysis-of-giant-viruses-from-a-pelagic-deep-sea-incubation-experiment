base_dir=$WORK/bird_bags/YuChens_pipeline

# get rRNA seq files from PhyloFlash output
mkdir $base_dir/Euk_abundance/rRNA_seqs
for folder in $(ls -d $base_dir/Euk_abundance/10_eukaryote_pf/*); do
    sample=$(echo $folder | rev | cut -f1 -d "/" | rev)
    echo $sample
    cp "$folder"/LIB.all.final.fasta $base_dir/Euk_abundance/rRNA_seqs/"$sample"_rRNA.fasta
done

# get pr2 database -> euks_pr2_database.Rmd

# blast
conda activate blast
for sample in $(ls $base_dir/Euk_abundance/rRNA_seqs/); do
    sbatch \
        --job-name=euk_rRNA_blast_$sample \
        --nodes=1 \
        --tasks-per-node=32 \
        --cpus-per-task=1 \
        --mem=120G \
        --time=02:00:00 \
        --output=$WORK/bird_bags/YuChens_pipeline/reports/euk_rRNA_blast_"$sample"_out \
        --error=$WORK/bird_bags/YuChens_pipeline/reports/euk_rRNA_blast_"$sample"_err \
        --wrap="\
            blastn \
                -query $base_dir/Euk_abundance/rRNA_seqs/$sample \
                -subject $WORK/bird_bags/db/pr2_5.0_20230619/pr2_database.fasta \
                -out $base_dir/Euk_abundance/rRNA_seqs/$sample.vs.pr2.tbl \
                -outfmt 6 \
                -max_target_seqs 1 \
                -num_threads 32 \
                "
done

for file in $(ls $base_dir/Euk_abundance/rRNA_seqs/*.tbl); do
    sample=$(echo $file | cut -f1 -d "." | rev | cut -f1 -d "/" | rev)
    awk -v sample=$sample '{print sample, $0}' $file \
        >> $base_dir/Euk_abundance/summary_rRNA.fasta.vs.pr2.tbl
done

# get coverage from phyloflash output
mkdir $base_dir/Euk_abundance/coverage_files
for folder in $(ls -d $base_dir/Euk_abundance/10_eukaryote_pf/*); do
    sample=$(echo $folder | rev | cut -f1 -d "/" | rev | cut -f1 -d ".")
    echo $sample
    awk -F ',' -v sample=$sample 'NR!=1 {print sample,$1,$3}' "$folder"/LIB.phyloFlash.extractedSSUclassifications.csv \
        >> $base_dir/Euk_abundance/summary_rRNA_coverage.tbl
done

# coverage for ssufinder from anvio on local
echo -e "sample\tcontig\tcoverage" > anvio_coverage/coverage_summary.txt
for file in $(ls anvio_coverage/anvio_*); do
    sample=$(echo $file | cut -f2 -d "/" | sed 's/anvio_//g')
    awk -v sample="$sample" 'BEGIN {OFS="\t"} NR>1 {print sample, $0}' $file \
        >> anvio_coverage/coverage_summary.txt
done