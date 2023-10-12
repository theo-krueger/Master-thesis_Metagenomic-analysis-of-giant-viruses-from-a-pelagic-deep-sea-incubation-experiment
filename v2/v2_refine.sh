#!/bin/sh

############ anvi-refine ######################
# ssh -L 8092:localhost:8092 -X smomw539@nesh-fe.rz.uni-kiel.de
# ssh -L 8093:localhost:8093 -X smomw539@nesh-fe.rz.uni-kiel.de

# working with bins from Bins_abundances_v1
conda activate anvio-7.1
anvio_dir=$WORK/bird_bags/YuChens_pipeline/05_anvio

anvi-refine \
    -c $anvio_dir/rename.5k.all.spades.contigs.fasta.db \
    -p $anvio_dir/merged_samples/PROFILE.db \
    -C collection_select_NCLDV_coAssem_concoct \
    --server-only -P 8092 \
    -b BB_coAssem_output_concoct_40_2

anvi-refine \
    -c $anvio_dir/rename.5k.all.spades.contigs.fasta.db \
    -p $anvio_dir/merged_samples/PROFILE.db \
    -C collection_select_NCLDV_ESOM \
    --server-only -P 8092 \
    -b BB_2016_1040m_12C_4d_14_L00m_ESOM_34_1

## refining based on bins that contain polB

# anvi-refine
mkdir $WORK/bird_bags/YuChens_pipeline/05_anvio/summaries/v2_polBbins_refined/
for collection in $(ls $WORK/bird_bags/YuChens_pipeline/05_anvio/collections/c* | cut -f9 -d '/'); do
    anvi-summarize \
        -c $WORK/bird_bags/YuChens_pipeline/05_anvio/rename.5k.all.spades.contigs.fasta.db \
        -p $WORK/bird_bags/YuChens_pipeline/05_anvio/merged_samples/PROFILE.db \
        -C $collection \
        -o $WORK/bird_bags/YuChens_pipeline/05_anvio/summaries/v2_polBbins_refined/$collection-output_v2_polB_bins_refined
done

## trying to find polB from Fabians assembly in the new one

# old polB seqs:
# 15_Bin.85 c_000001244625_8
# 15_Bin.85 c_000002621871_36
# 8_Bin.4 c_000000221285_9
# 8_Bin.4 c_000000407204_40
mkdir $WORK/bird_bags/comparison/prasinovirus_polB/
# copy to $WORK/bird_bags/comparison/prasinovirus_polB/possible_prasinovirus_polB_Fpipeline.faa
# c_000002621871_36 and c_000000407204_40 are 100% the same

# make YC_pipeline bins into db v1
mkdir $WORK/bird_bags/db/YC_pipeline/
mkdir $WORK/bird_bags/db/YC_pipeline/protein_seqs/
for file in $(find $WORK/bird_bags/YuChens_pipeline/02_predicted/v1/ -maxdepth 1 -type f); do
    cat $file >> $WORK/bird_bags/db/YC_pipeline/protein_seqs/v1_all_proteins.faa
done

conda activate diamond-2.1.4
diamond makedb --in $WORK/bird_bags/db/YC_pipeline/protein_seqs/v1_all_proteins.faa -d $WORK/bird_bags/db/YC_pipeline/protein_seqs/v1_all_proteins.dmnd

sbatch \
    --job-name=find_F_polB \
    --nodes=1 \
    --tasks-per-node=32 \
    --cpus-per-task=1 \
    --mem=120G \
    --time=24:00:00 \
    --output=$WORK/bird_bags/comparison/reports/find_F_polB.out \
    --error=$WORK/bird_bags/comparison/reports/find_F_polB.err \
    --mail-type=END,FAIL,TIME_LIMIT \
    --mail-user=theokrueger.infomail@geomar.de \
    --wrap="\
        diamond blastp \
            --query $WORK/bird_bags/comparison/prasinovirus_polB/possible_prasinovirus_polB_Fpipeline.faa \
            --db $WORK/bird_bags/db/YC_pipeline/protein_seqs/v1_all_proteins.dmnd \
            --out $WORK/bird_bags/comparison/prasinovirus_polB/F_prasinovirus_polB_vs_YC_v1.diamond \
            --threads 32 \
            --evalue 1e-04 \
            --outfmt 6 \
            --max-target-seqs 10
    "
cat $WORK/bird_bags/comparison/prasinovirus_polB/F_prasinovirus_polB_vs_YC_v1.diamond | sort | uniq \
    > $WORK/bird_bags/comparison/prasinovirus_polB/F_prasinovirus_polB_vs_YC_v1_uniq.diamond
awk '{print $2}' $WORK/bird_bags/comparison/prasinovirus_polB/F_prasinovirus_polB_vs_YC_v1_uniq.diamond | sort -u \
    > $WORK/bird_bags/comparison/prasinovirus_polB/F_prasinovirus_polB_vs_YC_v1_uniq_hits.txt

for seq in $(cat $WORK/bird_bags/comparison/prasinovirus_polB/F_prasinovirus_polB_vs_YC_v1_uniq_hits.txt); do
    echo "Searching $seq in v2:"
    grep $seq $WORK/bird_bags/YuChens_pipeline/02_predicted/v2/* | cut -f1 -d ":" | rev | cut -f1 -d "/"| rev
    echo -e '\n'
    echo "Searching $seq in v1:"
    grep $seq $WORK/bird_bags/YuChens_pipeline/02_predicted/v1/* | cut -f1 -d ":" | rev | cut -f1 -d "/"| rev
    echo -e '\n\n'
done

# make YC_pipeline bins into db v0
for file in $(find $WORK/bird_bags/YuChens_pipeline/02_predicted/v0/ -maxdepth 1 -type f); do
    cat $file >> $WORK/bird_bags/db/YC_pipeline/protein_seqs/v0_all_proteins.faa
done

conda activate diamond-2.1.4
diamond makedb --in $WORK/bird_bags/db/YC_pipeline/protein_seqs/v0_all_proteins.faa -d $WORK/bird_bags/db/YC_pipeline/protein_seqs/v0_all_proteins.dmnd

sbatch \
    --job-name=find_F_polB \
    --nodes=1 \
    --tasks-per-node=32 \
    --cpus-per-task=1 \
    --mem=120G \
    --time=24:00:00 \
    --output=$WORK/bird_bags/comparison/reports/find_F_polB.out \
    --error=$WORK/bird_bags/comparison/reports/find_F_polB.err \
    --mail-type=END,FAIL,TIME_LIMIT \
    --mail-user=theokrueger.infomail@geomar.de \
    --wrap="\
        diamond blastp \
            --query $WORK/bird_bags/comparison/prasinovirus_polB/possible_prasinovirus_polB_Fpipeline.faa \
            --db $WORK/bird_bags/db/YC_pipeline/protein_seqs/v0_all_proteins.dmnd \
            --out $WORK/bird_bags/comparison/prasinovirus_polB/F_prasinovirus_polB_vs_YC_v0.diamond \
            --threads 32 \
            --evalue 1e-04 \
            --outfmt 6 \
            --max-target-seqs 10
    "

cat $WORK/bird_bags/comparison/prasinovirus_polB/F_prasinovirus_polB_vs_YC_v0.diamond | sort | uniq \
    > $WORK/bird_bags/comparison/prasinovirus_polB/F_prasinovirus_polB_vs_YC_v0_uniq.diamond
awk '{print $2}' $WORK/bird_bags/comparison/prasinovirus_polB/F_prasinovirus_polB_vs_YC_v0_uniq.diamond | sort -u \
    > $WORK/bird_bags/comparison/prasinovirus_polB/F_prasinovirus_polB_vs_YC_v0_uniq_hits.txt
for seq in $(cat $WORK/bird_bags/comparison/prasinovirus_polB/F_prasinovirus_polB_vs_YC_v0_uniq_hits.txt); do
    echo "Searching $seq in v2:"
    grep $seq $WORK/bird_bags/YuChens_pipeline/02_predicted/v2/* | cut -f1 -d ":" | rev | cut -f1 -d "/"| rev
    echo -e '\n'
    echo "Searching $seq in v1:"
    grep $seq $WORK/bird_bags/YuChens_pipeline/02_predicted/v1/* | cut -f1 -d ":" | rev | cut -f1 -d "/"| rev
    echo -e '\n'
    echo "Searching $seq in v0:"
    grep $seq $WORK/bird_bags/YuChens_pipeline/02_predicted/v0/* | cut -f1 -d ":" | rev | cut -f1 -d "/"| rev
    echo -e '\n\n'
done