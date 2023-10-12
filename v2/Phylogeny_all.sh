#!/bin/bash

base_dir=$WORK/bird_bags/YuChens_pipeline
phylo_dir=$WORK/bird_bags/YuChens_pipeline/Phylogeny


mkdir $phylo_dir/allpolB_v1

# copy own polB extracted before
cp $WORK/bird_bags/YuChens_pipeline/Comparison_Schulz/v1/polB_contigs_YC_V1_shortheaders.faa \
    $phylo_dir/allpolB_v1

# add own to Needham alignment
conda activate mafft
sbatch \
    --job-name=align_mafft \
    --nodes=1 \
    --tasks-per-node=32 \
    --cpus-per-task=1 \
    --mem=10G \
    --time=0:30:00 \
    --output=$base_dir/reports/align_all_out \
    --error=$base_dir/reports/align_all_err \
    --mail-type=END,FAIL,TIME_LIMIT \
    --mail-user=theokrueger.infomail@gmail.com \
    --wrap="\
        mafft \
        --add $phylo_dir/allpolB_v1/polB_contigs_YC_V1_shortheaders.faa \
        --reorder $phylo_dir/polB_Needham_2019_aligned.fa \
        > $phylo_dir/allpolB_v1/polB_own_Needham2019_aligned.fa \
        "

# trimming
conda activate trimal
trimal \
    -in $phylo_dir/allpolB_v1/polB_own_Needham2019_aligned.fa \
    -out $phylo_dir/allpolB_v1/polB_own_Needham2019_aligned_trimmed.fa \
    -gappyout \
    -htmlout $phylo_dir/allpolB_v1/report_polB_own_Needham2019_aligned_trimmed.html

# treebuilding
conda activate iqtree
sbatch \
    --job-name=iqtree \
    --nodes=1 \
    --tasks-per-node=4 \
    --cpus-per-task=1 \
    --mem=150G \
    --time=24:00:00 \
    --output=$base_dir/reports/iqtree_all_out \
    --error=$base_dir/reports/iqtree_all_err \
    --mail-type=END,FAIL,TIME_LIMIT \
    --mail-user=theokrueger.infomail@gmail.com\
    --wrap="\
        iqtree \
        -s $phylo_dir/allpolB_v1/polB_own_Needham2019_aligned_trimmed.fa \
        -T AUTO \
        -B 1000 \
        "

### again with all found polB with reasonable quality
conda activate mafft
mafft \
    --add $phylo_dir/additional_genomes/04_polB/polB_additional_genomes_shortheaders.fa \
    --reorder $phylo_dir/polB_Needham_2019_aligned.fa \
    > $phylo_dir/polB_Needham2019_additional_aligned.fa 

# add own to Needham alignment
mafft \
    --add $phylo_dir/allpolB_v1_unique/polB_contigs_YC_V1_unique_shortheaders.faa \
    --reorder $phylo_dir/polB_Needham2019_additional_aligned.fa  \
    > $phylo_dir/allpolB_v1_unique/polB_own_unique_Needham2019_additional_aligned.fa


# trimming
conda activate trimal
trimal \
    -in $phylo_dir/allpolB_v1_unique/polB_own_unique_Needham2019_additional_aligned.fa \
    -out $phylo_dir/allpolB_v1_unique/polB_own_unique_Needham2019_additional_aligned_trimmed.fa \
    -gappyout \
    -htmlout $phylo_dir/allpolB_v1_unique/report_polB_own_unique_Needham2019_additional_aligned_trimmed.html

# tree building
conda activate iqtree
sbatch \
    --job-name=iqtree \
    --nodes=1 \
    --tasks-per-node=4 \
    --cpus-per-task=1 \
    --mem=150G \
    --time=24:00:00 \
    --output=$base_dir/reports/iqtree_all_out \
    --error=$base_dir/reports/iqtree_all_err \
    --mail-type=END,FAIL,TIME_LIMIT \
    --mail-user=theokrueger.infomail@gmail.com\
    --wrap="\
        iqtree \
        -s $phylo_dir/allpolB_v1_unique/polB_own_unique_Needham2019_additional_aligned_trimmed.fa \
        -T AUTO \
        -B 1000 \
        "


# stats
for seq in $(grep ">" $base_dir/Phylogeny/tree07_allpolB_v1_unique/polB_contigs_YC_V1_unique_shortheaders.faa); do
    contig=$(echo $seq | sed 's/>//g')

    # file=$(grep "$contig" $base_dir/01_DNA_sequences/v1_raw/*.fasta | \
    #     head -n1 | \
    #     cut -f1 -d ":")
    # bin=$(echo $file | cut -f9 -d "/")
    # echo $file $bin $contig >> $base_dir/Phylogeny/tree07_allpolB_v1_unique/own_polB_seqs.txt
    # echo $contig >> $base_dir/Phylogeny/tree07_allpolB_v1_unique/own_polB_seqs_contigs.txt

    grep "$contig" $base_dir/00_stats/v1/genomes_contig_lengths.txt | \
        head -n1 | cut -f2,3 -d $'\t' \
        >> $base_dir/Phylogeny/tree07_allpolB_v1_unique/own_polB_seqs_contiglengths.txt
done

#polB length
conda activate bioawk
bioawk -c fastx -t '{print $name, length($seq)}' $base_dir/Phylogeny/tree07_allpolB_v1_unique/polB_contigs_YC_V1_unique_shortheaders.faa \
    > $base_dir/Phylogeny/tree07_allpolB_v1_unique/own_polB_seqs_lengths.txt

# when were sequences removed
COUNTER=1
for contig in $(cat $base_dir/Phylogeny/tree07_allpolB_v1_unique/own_polB_seqs_contigs.txt); do
    
    echo "$COUNTER: Searching $contig in v2:"
    grep --exclude-dir=metabat2_only $contig $base_dir/01_DNA_sequences/v2/*
    echo " "

    echo "$COUNTER: Searching $contig in v1:"
    grep --exclude-dir=metabat2_only $contig $base_dir/01_DNA_sequences/v1_raw/*
    echo " "

    echo " "
    ((COUNTER++))
done
echo -e "contig\tbins_v2\tbins_v1" > $base_dir/Phylogeny/tree07_allpolB_v1_unique/own_polB_seqs_bins.txt
for contig in $(cat $base_dir/Phylogeny/tree07_allpolB_v1_unique/own_polB_seqs_contigs.txt); do
    v2=$(grep --exclude-dir=metabat2_only $contig $base_dir/01_DNA_sequences/v2/* | cut -f1 -d ":" | rev | cut -f1 -d "/"| rev | tr '\n' ';')
    v1=$(grep --exclude-dir=metabat2_only $contig $base_dir/01_DNA_sequences/v1_raw/* | cut -f1 -d ":" | rev | cut -f1 -d "/"| rev | tr '\n' ';')
    echo -e "$contig\t$v2\t$v1" >> $base_dir/Phylogeny/tree07_allpolB_v1_unique/own_polB_seqs_bins.txt
done

# file of protein ids
awk '{print $1}' $base_dir/Phylogeny/tree07_allpolB_v1_unique/own_polB_seqs_lengths.txt \
    > $base_dir/Phylogeny/tree07_allpolB_v1_unique/own_polB_seqs_protein_ids.txt

# get all proteins from polB contigs
mkdir $base_dir/Phylogeny/tree07_allpolB_v1_unique/polB_contigs_proteins/
conda activate seqtk
while read line; do
    contig=$(echo $line | cut -f1 -d $' ')
    file=$(echo $line | rev | cut -f2 -d ";" | rev | cut -f2 -d $' ')
    echo "Finding $file: $contig"
    sequences=$(grep "$contig" $base_dir/02_predicted/v1/$file.faa | cut -f1 -d " " | sed 's/>//g')
    echo $sequences | sed 's/ /\n/g' > read_file.temp
    seqtk subseq $base_dir/02_predicted/v1/$file.faa read_file.temp > $base_dir/Phylogeny/tree07_allpolB_v1_unique/polB_contigs_proteins/"$contig"_proteins.fa
done < $base_dir/Phylogeny/tree07_allpolB_v1_unique/own_polB_seqs_bins.txt

# get DNA seq (full contig)
conda activate seqtk
for file in $(ls $base_dir/01_DNA_sequences/v1_raw/*); do
    seqtk subseq $file $base_dir/Phylogeny/tree07_allpolB_v1_unique/own_polB_seqs_contigs.txt >> $base_dir/Phylogeny/tree07_allpolB_v1_unique/polB_contigsYC_V1.fa
done
# only unique
conda activate seqkit
cat $base_dir/Phylogeny/tree07_allpolB_v1_unique/polB_contigsYC_V1.fa | \
    seqkit rmdup -s -o $base_dir/Phylogeny/tree07_allpolB_v1_unique/polB_contigsYC_V1_unique.fa 

# get DNA seq (only polB)
conda activate seqkit
for seq_id in $(cat $base_dir/Phylogeny/tree07_allpolB_v1_unique/own_polB_seqs_protein_ids.txt); do
    output=$(grep -r "$seq_id" $WORK/bird_bags/YuChens_pipeline/02_predicted/v1/ | head -n1)
    
    file=$(echo $output | cut -f1 -d ":" | rev | cut -f1 -d "/" | rev | sed 's/.faa//g')
    bin=$(echo $output | cut -f2 -d ":" | cut -f1 -d "#" | rev | cut -f2- -d "_" | rev | sed 's/>//g')
    sstart=$(echo $output | cut -f2 -d ":" | cut -f2 -d "#" | xargs)
    send=$(echo $output | cut -f2 -d ":" | cut -f3 -d "#" | xargs)
    direction=$(echo $output | cut -f2 -d ":" | cut -f4 -d "#" | xargs)
    echo -e "$file\t$bin\t$sstart\t$send\t$direction\n"

    if [ $direction = "1" ]; then
        seq_nuc=$(seqkit faidx $base_dir/01_DNA_sequences/v1_raw/$file $bin:$sstart-$send)
    elif [ $direction = "-1" ]; then
        seq_nuc=$(seqkit faidx $base_dir/01_DNA_sequences/v1_raw/$file $bin:$send-$sstart)
    else
        echo "error"
    fi
    
    echo $seq_nuc
    echo $seq_nuc | tr ' ' '\n' \
        >> $base_dir/Phylogeny/tree07_allpolB_v1_unique/own_polB_seqs_nucleotide.fa
done

# n proteins on each polB contig
for seq in $(cat $base_dir/Phylogeny/tree07_allpolB_v1_unique/own_polB_seqs_protein_ids.txt); do
    contig=$(echo $seq | rev | cut -d "_" -f2- | rev)
    n_proteins=$(grep ">" $base_dir/Phylogeny/tree07_allpolB_v1_unique/polB_contigs_proteins/"$contig"_proteins.fa | wc -l)
    echo -e "$seq\t$n_proteins" >> $base_dir/Phylogeny/tree07_allpolB_v1_unique/own_polB_seqs_contigs_nproteins.txt
done


### new tree based on blasts / detailed look at alignment

mkdir $phylo_dir/tree08_allpolB_v1_unique_refined
cp $phylo_dir/tree07_allpolB_v1_unique/polB_contigs_YC_V1_unique_shortheaders.faa $phylo_dir/tree08_allpolB_v1_unique_refined/polB_contigs_YC_V1_unique_shortheaders.faa
# delete based on decisions -> lab book/ tables
# BB_2016_1040m_12C_4m_20_L002_000000000919_11
# BB_2016_1040m_12C_4m_20_L002_000000001389_12
# BB_2016_1040m_12C_4d_13_L00m_000000007167_4
# BB_2016_1040m_12C_4d_14_L00m_000000002372_9
# BB_2016_1040m_12C_4d_14_L00m_000000013819_5
# BB_2016_1040m_12C_4d_13_L00m_000000000113_29
# BB_2016_1040m_12C_4d_14_L00m_000000006182_12
# BB_2016_1040m_12C_4d_14_L00m_000000000391_31

# add
conda activate mafft
mafft \
    --add $phylo_dir/tree08_allpolB_v1_unique_refined/polB_contigs_YC_V1_unique_shortheaders.faa \
    --reorder $phylo_dir/polB_Needham2019_additional_aligned.fa  \
    > $phylo_dir/tree08_allpolB_v1_unique_refined/polB_own_unique_Needham2019_additional_aligned.fa
# trim
conda activate trimal
trimal \
    -in $phylo_dir/tree08_allpolB_v1_unique_refined/polB_own_unique_Needham2019_additional_aligned.fa \
    -out $phylo_dir/tree08_allpolB_v1_unique_refined/polB_own_unique_Needham2019_additional_aligned_trimmed.fa \
    -gappyout \
    -htmlout $phylo_dir/tree08_allpolB_v1_unique_refined/report_polB_own_unique_Needham2019_additional_aligned_trimmed.html
# tree building
# model selection based on previous model selections and waiting for the AIC for the most commonly used ones, then restarting the run
conda activate iqtree
sbatch \
    --job-name=iqtree \
    --nodes=1 \
    --tasks-per-node=32 \
    --cpus-per-task=1 \
    --mem=150G \
    --time=24:00:00 \
    --output=$base_dir/reports/iqtree_all_out \
    --error=$base_dir/reports/iqtree_all_err \
    --mail-type=END,FAIL,TIME_LIMIT \
    --mail-user=theokrueger.infomail@gmail.com\
    --wrap="\
        iqtree \
        -s $phylo_dir/tree08_allpolB_v1_unique_refined/polB_own_unique_Needham2019_additional_aligned_trimmed.fa \
        -m LG+F+R7 \
        -T AUTO \
        -B 1000 \
        "

### tree9

# add together
mkdir $WORK/bird_bags/YuChens_pipeline/Phylogeny/tree09_allpolB_new_samples
cat $WORK/bird_bags/YuChens_pipeline/Phylogeny/additional_genomes_new_sampling/polB_add_combined.fa \
    > $WORK/bird_bags/YuChens_pipeline/Phylogeny/tree09_allpolB_new_samples/polB_add_own_combined.fa
cat $WORK/bird_bags/YuChens_pipeline/Phylogeny/tree08_allpolB_v1_unique_refined/polB_contigs_YC_V1_unique_shortheaders.faa \
    >> $WORK/bird_bags/YuChens_pipeline/Phylogeny/tree09_allpolB_new_samples/polB_add_own_combined.fa
# remove duplicate seqs
conda activate seqkit
seqkit rmdup \
    -s \
    -o $WORK/bird_bags/YuChens_pipeline/Phylogeny/tree09_allpolB_new_samples/polB_add_own_combined_clean.fa \
    $WORK/bird_bags/YuChens_pipeline/Phylogeny/tree09_allpolB_new_samples/polB_add_own_combined.fa
# tree
conda activate mafft
mafft \
    --auto \
        $WORK/bird_bags/YuChens_pipeline/Phylogeny/tree09_allpolB_new_samples/polB_add_own_combined_clean.fa \
        > $WORK/bird_bags/YuChens_pipeline/Phylogeny/tree09_allpolB_new_samples/polB_add_own_combined_clean_aligned.fa
conda activate trimal
trimal \
    -in $WORK/bird_bags/YuChens_pipeline/Phylogeny/tree09_allpolB_new_samples/polB_add_own_combined_clean_aligned.fa \
    -out $WORK/bird_bags/YuChens_pipeline/Phylogeny/tree09_allpolB_new_samples/polB_add_own_combined_clean_aligned_trimmed.fa \
    -gappyout \
    -htmlout $WORK/bird_bags/YuChens_pipeline/Phylogeny/tree09_allpolB_new_samples/report_polB_add_own_combined_clean_aligned_trimmed.html
conda activate iqtree
sbatch \
    --job-name=iqtree \
    --nodes=1 \
    --tasks-per-node=32 \
    --cpus-per-task=1 \
    --mem=150G \
    --time=24:00:00 \
    --output=$WORK/bird_bags/YuChens_pipeline/reports/iqtree_all_out \
    --error=$WORK/bird_bags/YuChens_pipeline/reports/iqtree_all_err \
    --mail-type=END,FAIL,TIME_LIMIT \
    --mail-user=theokrueger.infomail@gmail.com\
    --wrap="\
        iqtree \
        -s $WORK/bird_bags/YuChens_pipeline/Phylogeny/tree09_allpolB_new_samples/polB_add_own_combined_clean_aligned_trimmed.fa \
        -m LG+F+R7 \
        -T AUTO \
        -B 1000 \
        "

### tree10
# add Needham back in
mkdir $WORK/bird_bags/YuChens_pipeline/Phylogeny/tree10_allpolB_new_samples
cp $WORK/bird_bags/YuChens_pipeline/Phylogeny/tree09_allpolB_new_samples/polB_add_own_combined_clean.fa \
    $WORK/bird_bags/YuChens_pipeline/Phylogeny/tree10_allpolB_new_samples/polB_add_own_combined.fa

# remove sequences from before, too short
# M-3300014204-43|Ga0172381_10000352_43
# M-3300017989-5|Ga0180432_10013138_6
# M-3300023184-184|Ga0214919_10020353_2
# S-1063923-109|1063923_contig_14058_1
# S-ERX556101-89|ERX556101_contig_183_1

conda activate mafft
mafft \
    --add $WORK/bird_bags/YuChens_pipeline/Phylogeny/tree10_allpolB_new_samples/polB_add_own_combined.fa \
    --reorder $WORK/bird_bags/YuChens_pipeline/Phylogeny/polB_Needham2019_additional_aligned.fa  \
    > $WORK/bird_bags/YuChens_pipeline/Phylogeny/tree10_allpolB_new_samples/polB_own_Needham2019_additional_aligned.fa
conda activate trimal
trimal \
    -in $WORK/bird_bags/YuChens_pipeline/Phylogeny/tree10_allpolB_new_samples/polB_own_Needham2019_additional_aligned.fa \
    -out $WORK/bird_bags/YuChens_pipeline/Phylogeny/tree10_allpolB_new_samples/polB_own_Needham2019_additional_aligned_trimmed.fa \
    -gappyout \
    -htmlout $WORK/bird_bags/YuChens_pipeline/Phylogeny/tree10_allpolB_new_samples/report_polB_own_Needham2019_additional_aligned_trimmed.html
conda activate iqtree
sbatch \
    --job-name=iqtree \
    --nodes=1 \
    --tasks-per-node=32 \
    --cpus-per-task=1 \
    --mem=150G \
    --time=24:00:00 \
    --output=$WORK/bird_bags/YuChens_pipeline/reports/iqtree_all_out \
    --error=$WORK/bird_bags/YuChens_pipeline/reports/iqtree_all_err \
    --mail-type=END,FAIL,TIME_LIMIT \
    --mail-user=theokrueger.infomail@gmail.com\
    --wrap="\
        iqtree \
        -s $WORK/bird_bags/YuChens_pipeline/Phylogeny/tree10_allpolB_new_samples/polB_own_Needham2019_additional_aligned_trimmed.fa \
        -m LG+F+R7 \
        -T AUTO \
        -B 1000 \
        "

# tree 11
mkdir $WORK/bird_bags/YuChens_pipeline/Phylogeny/tree11_allpolB_new_samples
cp $WORK/bird_bags/YuChens_pipeline/Phylogeny/tree10_allpolB_new_samples/polB_add_own_combined.fa \
    $WORK/bird_bags/YuChens_pipeline/Phylogeny/tree11_allpolB_new_samples/polB_add_own_combined.fa

# manually add Bachy 2021 sequences which are not already part of Needham 2019

# remove duplicates from same viruses
# NC_010191.2_224
# NC_014765.1_201

conda activate mafft
mafft \
    --add $WORK/bird_bags/YuChens_pipeline/Phylogeny/tree11_allpolB_new_samples/polB_add_own_combined.fa \
    --reorder $WORK/bird_bags/YuChens_pipeline/Phylogeny/polB_Needham2019_additional_aligned.fa  \
    > $WORK/bird_bags/YuChens_pipeline/Phylogeny/tree11_allpolB_new_samples/polB_add_own_combined_aligned.fa
conda activate trimal
trimal \
    -in $WORK/bird_bags/YuChens_pipeline/Phylogeny/tree11_allpolB_new_samples/polB_add_own_combined_aligned.fa \
    -out $WORK/bird_bags/YuChens_pipeline/Phylogeny/tree11_allpolB_new_samples/polB_add_own_combined_aligned_trimmed.fa \
    -gappyout \
    -htmlout $WORK/bird_bags/YuChens_pipeline/Phylogeny/tree11_allpolB_new_samples/report_polB_add_own_combined_aligned_trimmed.html
conda activate iqtree
sbatch \
    --job-name=iqtree \
    --nodes=1 \
    --tasks-per-node=1 \
    --cpus-per-task=32 \
    --mem=150G \
    --time=24:00:00 \
    --output=$WORK/bird_bags/YuChens_pipeline/reports/iqtree_all_out \
    --error=$WORK/bird_bags/YuChens_pipeline/reports/iqtree_all_err \
    --mail-type=END,FAIL,TIME_LIMIT \
    --mail-user=theokrueger.infomail@gmail.com\
    --wrap="\
        iqtree \
        -s $WORK/bird_bags/YuChens_pipeline/Phylogeny/tree11_allpolB_new_samples/polB_add_own_combined_aligned_trimmed.fa \
        -m LG+F+R7 \
        -T AUTO \
        -B 1000 \
        "

### tree 12
mkdir $WORK/bird_bags/YuChens_pipeline/Phylogeny/tree12_allpolB_new_samples

# rerun with 1000 bootstraps
cp $WORK/bird_bags/YuChens_pipeline/Phylogeny/tree11_allpolB_new_samples/polB_add_own_combined_aligned.fa \
    $WORK/bird_bags/YuChens_pipeline/Phylogeny/tree12_allpolB_new_samples/
cp $WORK/bird_bags/YuChens_pipeline/Phylogeny/tree11_allpolB_new_samples/polB_add_own_combined_aligned_trimmed.fa \
    $WORK/bird_bags/YuChens_pipeline/Phylogeny/tree12_allpolB_new_samples/

sbatch \
    --job-name=iqtree \
    --nodes=1 \
    --tasks-per-node=32 \
    --cpus-per-task=1 \
    --mem=150G \
    --time=24:00:00 \
    --output=$WORK/bird_bags/YuChens_pipeline/reports/iqtree_all_out \
    --error=$WORK/bird_bags/YuChens_pipeline/reports/iqtree_all_err \
    --mail-type=END,FAIL,TIME_LIMIT \
    --mail-user=theokrueger.infomail@gmail.com\
    --wrap="\
        iqtree \
        -s $WORK/bird_bags/YuChens_pipeline/Phylogeny/tree12_allpolB_new_samples/polB_add_own_combined_aligned_trimmed.fa \
        -m LG+F+R7 \
        -T AUTO \
        -B 10000 \
        "

# tree 13
mkdir $WORK/bird_bags/YuChens_pipeline/Phylogeny/tree13_allpolB_new_samples

cp $WORK/bird_bags/YuChens_pipeline/Phylogeny/tree11_allpolB_new_samples/polB_add_own_combined.fa \
    $WORK/bird_bags/YuChens_pipeline/Phylogeny/tree13_allpolB_new_samples/polB_add_own_combined.fa
# manually add Baudoux2015

conda activate mafft
mafft \
    --add $WORK/bird_bags/YuChens_pipeline/Phylogeny/tree13_allpolB_new_samples/polB_add_own_combined.fa \
    --reorder $WORK/bird_bags/YuChens_pipeline/Phylogeny/polB_Needham2019_additional_aligned.fa  \
    > $WORK/bird_bags/YuChens_pipeline/Phylogeny/tree13_allpolB_new_samples/polB_add_own_combined_aligned.fa
conda activate trimal
trimal \
    -in $WORK/bird_bags/YuChens_pipeline/Phylogeny/tree13_allpolB_new_samples/polB_add_own_combined_aligned.fa \
    -out $WORK/bird_bags/YuChens_pipeline/Phylogeny/tree13_allpolB_new_samples/polB_add_own_combined_aligned_trimmed.fa \
    -gappyout \
    -htmlout $WORK/bird_bags/YuChens_pipeline/Phylogeny/tree13_allpolB_new_samples/report_polB_add_own_combined_aligned_trimmed.html
conda activate iqtree
sbatch \
    --job-name=iqtree \
    --nodes=1 \
    --tasks-per-node=1 \
    --cpus-per-task=32 \
    --mem=150G \
    --time=24:00:00 \
    --output=$WORK/bird_bags/YuChens_pipeline/reports/iqtree_all_out \
    --error=$WORK/bird_bags/YuChens_pipeline/reports/iqtree_all_err \
    --mail-type=END,FAIL,TIME_LIMIT \
    --mail-user=theokrueger.infomail@gmail.com\
    --wrap="\
        iqtree \
        -s $WORK/bird_bags/YuChens_pipeline/Phylogeny/tree13_allpolB_new_samples/polB_add_own_combined_aligned_trimmed.fa \
        -m LG+F+R7 \
        -T AUTO \
        -B 10000 \
        "
