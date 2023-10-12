# list of contigs
grep ">" $WORK/bird_bags/YuChens_pipeline/01_DNA_sequences/v2/BB_coAssem_metabat2_bins_94_1-contigs.fa > random_smaller_files/GroupE/contigs_coAssem_metabat2_94_1
# same with other

# find non-overlapping elements in two files
awk 'NR==FNR{a[$0];next}!($0 in a)' \
    random_smaller_files/GroupE/contigs_coAssem_metabat2_94_1 \
    random_smaller_files/GroupE/contigs_ESOM_4_1

# get diamond of both contigs
grep "BB_2016_1040m_12C_4m_19_L002_000000001099" 03_annotated/v2/NR/BB_2016_1040m_12C_4m_19_L002_ESOM_4_1-contigs.fa.faa.vs.NR.diamond

for contig in $(cat random_smaller_files/GroupB/contigs_in_coA_metabat_not_on_ESOM); do
    new=$(echo $contig | sed 's/>//g')
    grep $new 03_annotated/v2/NR/BB_coAssem_metabat2_bins_281_1-contigs.fa.faa.vs.NR.diamond \
        >> random_smaller_files/GroupB/contigs_in_coA_metabat_not_in_ESOM_annotation.diamond
done

# making a folder to have a better overview
base_dir=$WORK/bird_bags/YuChens_pipeline/

mkdir $WORK/bird_bags/YuChens_pipeline/03_annotated/v2_ESOM
mkdir $WORK/bird_bags/YuChens_pipeline/03_annotated/v2_ESOM/anvio
mkdir $WORK/bird_bags/YuChens_pipeline/03_annotated/v2_ESOM/crispr
mkdir $WORK/bird_bags/YuChens_pipeline/03_annotated/v2_ESOM/NCVOGs
mkdir $WORK/bird_bags/YuChens_pipeline/03_annotated/v2_ESOM/NR
mkdir $WORK/bird_bags/YuChens_pipeline/03_annotated/v2_ESOM/eggnog

mkdir $WORK/bird_bags/YuChens_pipeline/02_predicted/v2_ESOM
mkdir $WORK/bird_bags/YuChens_pipeline/01_DNA_sequences/v2_ESOM

for file in $(ls $base_dir/01_DNA_sequences/v2/); do
    if [[ "$file" == *"ESOM"* ]]; then
        cp $base_dir/01_DNA_sequences/v2/$file $base_dir/01_DNA_sequences/v2_ESOM/
        cp $base_dir/02_predicted/v2/$file.faa $base_dir/02_predicted/v2_ESOM/
    fi
done

# getting 32 output for sharing
mkdir $WORK/NR_32
for bin in $(cat $WORK/bird_bags/YuChens_pipeline/00_stats/v2/v2_bins_afterfilter_above_25_uniqNCVOG.txt); do
    cp $WORK/bird_bags/YuChens_pipeline/03_annotated/v2/NR/$bin-contigs.fa.faa.vs.NR.diamond \
        $WORK/NR_32/
done

