#!/bin/sh

#make new db
sbatch \
    --job-name=makedb \
    --nodes=1 \
    --tasks-per-node=32 \
    --cpus-per-task=1 \
    --mem=120G \
    --time=24:00:00 \
    --output=$WORK/bird_bags/db/reports/makedb_out \
    --error=$WORK/bird_bags/db/reports/makedb_err \
    --mail-type=END,FAIL,TIME_LIMIT \
    --mail-user=theokrueger.infomail@gmail.com \
    --wrap="\
        source $HOME/miniconda3/bin/activate diamond-2.1.4

        diamond makedb \
            --in $WORK/bird_bags/db/NCBI_vir_20220927/nr.gz \
            --db $WORK/bird_bags/db/NCBI_vir_20220927_WithTaxonomy/NCBInr_WTax.dmnd \
            --taxonmap $WORK/bird_bags/db/NCBI_vir_20220927_WithTaxonomy/prot.accession2taxid.FULL.gz \
            --taxonnodes $WORK/bird_bags/db/NCBI_vir_20220927_WithTaxonomy/nodes.dmp \
            --taxonnames $WORK/bird_bags/db/NCBI_vir_20220927_WithTaxonomy/names.dmp \
            "


# run blast
for bin in $(ls $WORK/bird_bags/YuChens_pipeline/02_predicted/v2_polB_bins/); do
    sbatch \
    --job-name=diamond_$bin \
    --nodes=1 \
    --tasks-per-node=32 \
    --cpus-per-task=1 \
    --mem=120G \
    --time=12:00:00 \
    --output=$WORK/bird_bags/YuChens_pipeline/reports/diamond_"$bin"_out \
    --error=$WORK/bird_bags/YuChens_pipeline/reports/diamond_"$bin"_err \
    --mail-type=END,FAIL,TIME_LIMIT \
    --mail-user=theokrueger.infomail@gmail.com \
    --wrap="\
        diamond blastp \
            --query $WORK/bird_bags/YuChens_pipeline/02_predicted/v2/$bin \
            --db $WORK/bird_bags/db/NCBI_vir_20220927_WithTaxonomy/NCBInr_WTax.dmnd \
            --out $WORK/bird_bags/YuChens_pipeline/03_annotated/v2_polB_bins/NR/$bin.vs.NRwTax.diamond \
            --threads 32 \
            --evalue 1e-04 \
            --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle staxids sscinames sskingdoms skingdoms sphylums \
            --max-target-seqs 10 \
            "
done