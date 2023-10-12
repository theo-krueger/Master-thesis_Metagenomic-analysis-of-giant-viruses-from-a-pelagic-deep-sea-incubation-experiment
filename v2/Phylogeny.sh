#!/bin/bash

base_dir=$WORK/bird_bags/YuChens_pipeline
phylo_dir=$WORK/bird_bags/YuChens_pipeline/Phylogeny


# copy own polB extracted before
cp $WORK/bird_bags/YuChens_pipeline/Comparison_Schulz/v2_32/polB_contigs_YC_32_shortheaders.faa \
    $phylo_dir/

# copy Davids polB
cp $WORK/bird_bags/Fabians_pipeline/Phylogenomics/Needham2019_TargetedMetagenomics_sequences/OG0000022_DNA_polymerase_elongation_subunit_family_B_NCVOG0038_linsi_aligned.fasta \
    $phylo_dir/polB_Needham_2019_aligned.fa

mkdir $phylo_dir/tree01

# add own to Needham alignment
conda activate mafft
sbatch \
    --job-name=align_mafft \
    --nodes=1 \
    --tasks-per-node=32 \
    --cpus-per-task=1 \
    --mem=10G \
    --time=0:30:00 \
    --output=$base_dir/reports/align_out \
    --error=$base_dir/reports/align_err \
    --mail-type=END,FAIL,TIME_LIMIT \
    --mail-user=theokrueger.infomail@gmail.com \
    --wrap="\
        mafft \
        --add $phylo_dir/polB_contigs_YC_32_shortheaders.faa \
        --reorder $phylo_dir/polB_Needham_2019_aligned.fa \
        > $phylo_dir/tree01/polB_own_Needham2019_aligned.fa \
        "


# trimming
conda activate trimal
trimal \
    -in $phylo_dir/tree01/polB_own_Needham2019_aligned.fa \
    -out $phylo_dir/tree01/polB_own_Needham2019_aligned_trimmed.fa \
    -gappyout \
    -htmlout $phylo_dir/tree01/report_polB_own_Needham2019_aligned_trimmed.html

# treebuilding
conda activate iqtree
sbatch \
    --job-name=iqtree \
    --nodes=1 \
    --tasks-per-node=4 \
    --cpus-per-task=1 \
    --mem=150G \
    --time=24:00:00 \
    --output=$base_dir/reports/iqtree_out \
    --error=$base_dir/reports/iqtree_err \
    --mail-type=END,FAIL,TIME_LIMIT \
    --mail-user=theokrueger.infomail@gmail.com\
    --wrap="\
        iqtree \
        -s $phylo_dir/tree01/polB_own_Needham2019_aligned_trimmed.fa \
        -T AUTO \
        -B 1000 \
        "


### tree2
# same alignment as tree1
conda activate trimal
trimal \
    -in $phylo_dir/tree02/polB_own_Needham2019_aligned.fa \
    -out $phylo_dir/tree02/polB_own_Needham2019_aligned_trimmed.fa \
    -gt 0.8 \
    -htmlout $phylo_dir/tree02/report_polB_own_Needham2019_aligned_trimmed.html

# treebuilding
conda activate iqtree
sbatch \
    --job-name=iqtree2 \
    --nodes=1 \
    --tasks-per-node=4 \
    --cpus-per-task=1 \
    --mem=150G \
    --time=12:00:00 \
    --output=$base_dir/reports/iqtree_2_out \
    --error=$base_dir/reports/iqtree_2_err \
    --mail-type=END,FAIL,TIME_LIMIT \
    --mail-user=theokrueger.infomail@gmail.com\
    --wrap="\
        iqtree \
        -s $phylo_dir/tree02/polB_own_Needham2019_aligned_trimmed.fa \
        -T AUTO \
        -B 1000 \
        "

### tree3
# only ESOM

# add to Needham
conda activate mafft
sbatch \
    --job-name=align_mafft \
    --nodes=1 \
    --tasks-per-node=32 \
    --cpus-per-task=1 \
    --mem=10G \
    --time=0:30:00 \
    --output=$base_dir/reports/align_out \
    --error=$base_dir/reports/align_err \
    --mail-type=END,FAIL,TIME_LIMIT \
    --mail-user=theokrueger.infomail@gmail.com \
    --wrap="\
        mafft \
        --add $phylo_dir/polB_contigs_YC_ESOM_shortheaders.faa \
        --reorder $phylo_dir/polB_Needham_2019_aligned.fa \
        > $phylo_dir/tree03/polB_own_ESOM_Needham2019_aligned.fa \
        "

# trimming
conda activate trimal
trimal \
    -in $phylo_dir/tree03/polB_own_ESOM_Needham2019_aligned.fa \
    -out $phylo_dir/tree03/polB_own_Needham2019_aligned_trimmed.fa \
    -gappyout \
    -htmlout $phylo_dir/tree03/report_polB_own_Needham2019_aligned_trimmed.html

# treebuilding
conda activate iqtree
sbatch \
    --job-name=iqtree3 \
    --nodes=1 \
    --tasks-per-node=4 \
    --cpus-per-task=1 \
    --mem=150G \
    --time=12:00:00 \
    --output=$base_dir/reports/iqtree_3_out \
    --error=$base_dir/reports/iqtree_3_err \
    --mail-type=END,FAIL,TIME_LIMIT \
    --mail-user=theokrueger.infomail@gmail.com\
    --wrap="\
        iqtree \
        -s $phylo_dir/tree03/polB_own_Needham2019_aligned_trimmed.fa \
        -T AUTO \
        -B 1000 \
        "

### tree4
# more genomes
for file in $(ls $phylo_dir/additional_genomes/03_NCVOGs_annotated/); do
    awk -v file=$file '$0 ~ /^NCVOG0038/ {print file, $3, $5}' $phylo_dir/additional_genomes/03_NCVOGs_annotated/$file \
        >> $phylo_dir/additional_genomes/add_genomes_polB_stats.txt
done
# throwing out bad/ multiple hits

# extract polB
while read line; do 
    file=$(echo $line | cut -f1 -d ' ' | sed 's/.NCVOG_tblout//g')
    bin=$(echo $file | cut -f1 -d '.')
    contig=$(echo $line | cut -f2 -d ' ')
    grep -w "$contig" $phylo_dir/additional_genomes/02_predicted/$file | \
        sed 's/>//g' | \
        seqtk subseq $phylo_dir/additional_genomes/02_predicted/$file - | \
        awk -v pat="$contig" -v rep="${bin}-${contig}" '{gsub(pat, rep)} 1' \
        >> $phylo_dir/additional_genomes/04_polB/polB_additional_genomes.fa
done < $phylo_dir/additional_genomes/add_genomes_polB_stats.txt
# shorten headers
awk -F " " '{if ($0 ~ /^>/) { print $1;} else { print $0}}' $phylo_dir/additional_genomes/04_polB/polB_additional_genomes.fa \
    > $phylo_dir/additional_genomes/04_polB/polB_additional_genomes_shortheaders.fa

# add own + Needham alignment
conda activate mafft
sbatch \
    --job-name=align_mafft \
    --nodes=1 \
    --tasks-per-node=32 \
    --cpus-per-task=1 \
    --mem=10G \
    --time=0:30:00 \
    --output=$base_dir/reports/align_out \
    --error=$base_dir/reports/align_err \
    --mail-type=END,FAIL,TIME_LIMIT \
    --mail-user=theokrueger.infomail@gmail.com \
    --wrap="\
        mafft \
        --add $phylo_dir/additional_genomes/04_polB/polB_additional_genomes_shortheaders.fa \
        --reorder $phylo_dir/tree03/polB_own_ESOM_Needham2019_aligned.fa \
        > $phylo_dir/tree04/polB_own_Needham2019_additional_aligned.fa \
        "

# trimming
conda activate trimal
trimal \
    -in $phylo_dir/tree04/polB_own_Needham2019_additional_aligned.fa \
    -out $phylo_dir/tree04/polB_own_Needham2019_additional_aligned_trimmed.fa \
    -gappyout \
    -htmlout $phylo_dir/tree04/report_polB_own_Needham2019_additional_aligned.html

# treebuilding
conda activate iqtree
sbatch \
    --job-name=iqtree4 \
    --nodes=1 \
    --tasks-per-node=4 \
    --cpus-per-task=1 \
    --mem=150G \
    --time=12:00:00 \
    --output=$base_dir/reports/iqtree_4_out \
    --error=$base_dir/reports/iqtree_4_err \
    --mail-type=END,FAIL,TIME_LIMIT \
    --mail-user=theokrueger.infomail@gmail.com\
    --wrap="\
        iqtree \
        -s $phylo_dir/tree04/polB_own_Needham2019_additional_aligned_trimmed.fa \
        -T AUTO \
        -B 1000 \
        "
# download and change in figtree

### tree5
# remove LCMAC102, LCMAC103: can't find source, plot weird in the tree

# trimming
conda activate trimal
trimal \
    -in $phylo_dir/tree05/polB_own_Needham2019_additional_aligned.fa \
    -out $phylo_dir/tree05/polB_own_Needham2019_additional_aligned_trimmed.fa \
    -gappyout \
    -htmlout $phylo_dir/tree05/report_polB_own_Needham2019_additional_aligned.html

# treebuilding
conda activate iqtree
sbatch \
    --job-name=iqtree5 \
    --nodes=1 \
    --tasks-per-node=4 \
    --cpus-per-task=1 \
    --mem=150G \
    --time=12:00:00 \
    --output=$base_dir/reports/iqtree_5_out \
    --error=$base_dir/reports/iqtree_5_err \
    --mail-type=END,FAIL,TIME_LIMIT \
    --mail-user=theokrueger.infomail@gmail.com\
    --wrap="\
        iqtree \
        -s $phylo_dir/tree05/polB_own_Needham2019_additional_aligned_trimmed.fa \
        -T AUTO \
        -B 1000 \
        "