YC_base_dir=$WORK/bird_bags/YuChens_pipeline/
for DNA_file in $(ls $YC_base_dir/01_DNA_sequences/v2_5MAGs/); do
    bin=$(echo $DNA_file | cut -f1 -d ".")
    name=$(echo $bin | cut -f1 -d "-" | rev | cut -f1,2 -d "_" | rev)
    name="BB_$name"

    python3 $YC_base_dir/circos/circos/circos_main.py \
        -ig $YC_base_dir/01_DNA_sequences/v2_5MAGs/$DNA_file \
        -ip $YC_base_dir/02_predicted/v2_polB_bins/"$DNA_file".faa \
        -ia $YC_base_dir/03_annotated/v2_polB_bins/NR/"$DNA_file".faa.vs.NRwTax.diamond \
        -itdb $WORK/bird_bags/db/NCBI_vir_20220927_WithTaxonomy/ncbi_lineages_2023-04-04.csv \
        -ic $YC_base_dir/03_annotated/v2_polB_bins/crispr/"$bin"_CRISPRS.fasta \
        -in $YC_base_dir/03_annotated/v2_polB_bins/NCVOGs/$DNA_file.faa.NCVOG_tblout\
        -n $name \
        -o $YC_base_dir/circos
done

