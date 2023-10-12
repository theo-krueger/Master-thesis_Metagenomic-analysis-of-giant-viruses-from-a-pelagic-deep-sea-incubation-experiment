############# CRISPR #############

crt_dir=$HOME/miniconda3/envs/crt-mod-2.0/bin
cd $crt_dir/crt-mod-master/
source $HOME/miniconda3/bin/activate crt-mod-2.0

in_dir=$WORK/bird_bags/Comparison_YuChens_pipeline/01_DNA_sequences/v2
out_dir=$WORK/bird_bags/Comparison_YuChens_pipeline/03_annotated/v2/crispr
for file in $(ls $in_dir); do
    bin=$(echo $file | cut -f1 -d '.')
    python crt-mod.py \
        -i $in_dir/$file fasta \
        -o $out_dir/txt/"$bin"_CRISPRS.txt
done

# summarizing
echo -e "bin\tposition\trepeat\tspacer" >$WORK/bird_bags/YuChens_pipeline/03_annotated/v2/crispr/crispr_summary.table
for file in $(ls $WORK/bird_bags/YuChens_pipeline/03_annotated/v2/crispr/txt/); do
    bin=$(echo $file | sed 's/-contigs_CRISPRS.txt//g')
    while read line; do
        echo $line
        awk -v bin=$bin '{
        if ($0 ~ /^[0-9]/) 
            {
            printf("%s\t%s\n", bin, $0)
            }
        }' \
            >>$WORK/bird_bags/YuChens_pipeline/03_annotated/v2/crispr/crispr_summary.table
    done <$WORK/bird_bags/YuChens_pipeline/03_annotated/v2/crispr/txt/$file
done

# get only repeats for further analysis
for file in $(ls $WORK/bird_bags/YuChens_pipeline/03_annotated/v2/crispr/txt/); do
    bin=$(echo $file | sed 's/_CRISPRS.txt//g')
    while read line; do
        awk -v bin=$bin '{
        if ($0 ~ /^[0-9]/) 
            {
            print $2
            }
        }' \
            >>$WORK/bird_bags/YuChens_pipeline/03_annotated/v2/crispr/repeats/"$bin"_repeats.txt
    done <$WORK/bird_bags/YuChens_pipeline/03_annotated/v2/crispr/txt/$file
done

# continue in R with CRISPRclassify
# command line doesn't work, libraries not being updated
# reupload output to server

# summarizing
echo "Bin,Repeat,Subtype,Probability,Closest_Strain,Closest_Strain_Accession,Closest_Strain_Location,Closest_Strain_Subtype,Closest_Strain_Repeat,Edit_Dist" \
    >$WORK/bird_bags/YuChens_pipeline/03_annotated/v2/crispr/classification_summary.table
for file in $(ls $WORK/bird_bags/YuChens_pipeline/03_annotated/v2/crispr/classification/); do
    bin=$(echo $file | sed 's/_repeats.txt.crclass//g')
    while read line; do
        awk -v bin=$bin '{printf("%s,%s\n", bin, $0)}' | sed 's/"//g' \
            >>$WORK/bird_bags/YuChens_pipeline/03_annotated/v2/crispr/classification_summary.csv
    done <$WORK/bird_bags/YuChens_pipeline/03_annotated/v2/crispr/classification/$file
done
