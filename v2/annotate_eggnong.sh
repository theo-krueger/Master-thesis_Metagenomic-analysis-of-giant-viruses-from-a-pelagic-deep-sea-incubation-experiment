########### eggnog ##############

for bin in $(ls $WORK/bird_bags/YuChens_pipeline/02_predicted/v2/); do
    sbatch --job-name=$bin.eggnog \
        --nodes=1 \
        --tasks-per-node=32 \
        --cpus-per-task=1 \
        --mem=120G \
        --time=24:00:00 \
        --output=$WORK/bird_bags/YuChens_pipeline/reports/eggnog_out \
        --error=$WORK/bird_bags/YuChens_pipeline/reports/eggnog_err \
        --wrap="\
            python $HOME/miniconda3/envs/eggnog/bin/emapper.py \
            --output "$bin"_eggnog.tbl \
            --output_dir $WORK/bird_bags/YuChens_pipeline/03_annotated/v2/eggnog/ \
            --cpu 32 \
            -i $WORK/bird_bags/YuChens_pipeline/02_predicted/v2/$bin \
            -m diamond \
            --data_dir $WORK/bird_bags/db/eggnog_20230331/ \
            --dmnd_db $WORK/bird_bags/db/eggnog_20230331/eggnog_proteins.dmnd"
done

mkdir 03_annotated/v2/eggnog/annotations/ && mv 03_annotated/v2/eggnog/*.annotations 03_annotated/v2/eggnog/annotations/
mkdir 03_annotated/v2/eggnog/hits/ && mv 03_annotated/v2/eggnog/*.hits 03_annotated/v2/eggnog/hits/
mkdir 03_annotated/v2/eggnog/seed_orthologs/ && mv 03_annotated/v2/eggnog/*logs 03_annotated/v2/eggnog/seed_orthologs/

for file in $(ls 03_annotated/v2/eggnog/annotations/*.annotations); do
    bin=$(echo $file | cut -f5 -d '/' | cut -f1 -d '-')
    tail -n +5 $file | head -n -4 | sed 's/#query/query/g' \
        >03_annotated/v2/eggnog/annotations/clean/$bin-contigs.fa.faa_eggnog_clean.tbl.emapper.annotations.tsv
done
