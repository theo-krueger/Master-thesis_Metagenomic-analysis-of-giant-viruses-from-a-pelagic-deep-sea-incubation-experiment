######### anvio ###########
## contigs db
data_YC=$WORK/bird_bags/YuChens_pipeline
for file in $(ls $data_YC/01_DNA_sequences/v2); do
    bin=$(echo $file | cut -f1 -d '-')
    sbatch --job-name=contigsdb \
        --nodes=1 \
        --tasks-per-node=32 \
        --cpus-per-task=1 \
        --mem=120G \
        --time=02:00:00 \
        --output=$WORK/bird_bags/YuChens_pipeline/reports/contigsdb_out \
        --error=$WORK/bird_bags/YuChens_pipeline/reports/contigsdb_err \
        --wrap="\
            anvi-gen-contigs-database \
            -f $data_YC/01_DNA_sequences/v2/$file \
            -o $data_YC/05_anvio/contigs_db/"$bin"_contigs.db \
            -T 32"
done

# setup annot db
anvi-setup-kegg-kofams

# annotate
data_YC=$WORK/bird_bags/YuChens_pipeline
for file in $(ls $data_YC/05_anvio/contigs_db/); do
    sbatch --job-name=anvi_annot \
        --nodes=1 \
        --tasks-per-node=32 \
        --cpus-per-task=1 \
        --mem=120G \
        --time=02:00:00 \
        --output=$WORK/bird_bags/YuChens_pipeline/reports/anviannot_out \
        --error=$WORK/bird_bags/YuChens_pipeline/reports/anviannot_err \
        --wrap="\
            anvi-run-kegg-kofams \
            -c $data_YC/05_anvio/contigs_db/$file \
            -T 32"
done

# estimate
data_YC=$WORK/bird_bags/YuChens_pipeline
for file in $(ls $data_YC/05_anvio/contigs_db/); do
    bin=$(echo $file | sed 's/_contigs.db//g')
    sbatch \
        --job-name=anvi_estimate \
        --nodes=1 \
        --tasks-per-node=32 \
        --cpus-per-task=1 \
        --mem=120G \
        --time=02:00:00 \
        --output=$WORK/bird_bags/YuChens_pipeline/reports/anviestimate_"$file"_out \
        --error=$WORK/bird_bags/YuChens_pipeline/reports/anviestimate_"$file"_err \
        --wrap="\
            anvi-estimate-metabolism \
            -c $data_YC/05_anvio/contigs_db/$file \
            -O $data_YC/03_annotated/v2/anvio/output/$bin \
            "
done

# remove those that were kicked out by v2 filter
while read bin; do
    echo $bin"_modules.txt" >> $data_YC/03_annotated/v2/anvio/v2_bins_remaining_afterfilter_modules.txt
done < $data_YC/00_stats/v2/v2_bins_remaining_afterfilter.txt
for file in $(ls $data_YC/03_annotated/v2/anvio/output/); do
    if ! grep -qxFe "$file" $data_YC/03_annotated/v2/anvio/v2_bins_remaining_afterfilter_modules.txt; then
        echo "Deleting: $file"
        rm $data_YC/03_annotated/v2/anvio/output/$file
    fi
done

# add together divided by method
for file in $(ls $data_YC/03_annotated/v2/anvio/output/); do
    if [[ "$file" == *"coAssem_MaxBin"* ]]; then
        method=coAssem_MaxBin
    elif [[ "$file" == *"coAssem_metabat2"* ]]; then
        method=coAssem_metabat2
    elif [[ "$file" == *"coAssem_output_concoct"* ]]; then
        method=coAssem_concoct
    elif [[ "$file" == *"ESOM"* ]]; then
        method=ESOM
    elif [[ "$file" == *"concoct"* ]]; then
        method=concoct
    elif [[ "$file" == *"metabat2"* ]]; then
        method=metabat2
    elif [[ "$file" == *"MaxBin"* ]]; then
        method=MaxBin
    else
        echo 'Failed to read file: ' $file
        echo '  '
    fi 

    cat $data_YC/03_annotated/v2/anvio/output/$file \
        >> $data_YC/03_annotated/v2/anvio/output_combined/pathways_all_"$method".txt

done
for file in $(find $data_YC/03_annotated/v2/anvio/output_combined/ -maxdepth 1 -type f); do
    awk 'NR==1 || $0 !~ /^unique_id/ {print}' $file > temp && mv temp $file
done

# using own script to make bubble plots offline
# online doesnt work again

Rscript /Users/tkrueger/Google\ Drive/My\ Drive/AG_Worden/BirdBags/Code/Pathway_plots_new.R \
    -i /Users/tkrueger/Google\ Drive/My\ Drive/AG_Worden/BirdBags/YuChens_pipeline/v2_ESOM/anvio_pathways/pathways_ESOM_selected_increasing.txt \
    -o /Users/tkrueger/Google\ Drive/My\ Drive/AG_Worden/BirdBags/YuChens_pipeline/v2_ESOM/anvio_pathways/
