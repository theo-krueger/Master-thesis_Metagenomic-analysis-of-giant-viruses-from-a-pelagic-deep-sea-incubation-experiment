
base_dir=$WORK/bird_bags/YuChens_pipeline

############# NCVOGs #################
for bin in $(ls $base_dir/02_predicted/v2/); do
    sbatch \
        --job-name=$bin.hmm \
        --nodes=1 \
        --tasks-per-node=32 \
        --cpus-per-task=1 \
        --mem=120G \
        --time=00:10:00 \
        --output=$base_dir/reports/hmm_out \
        --error=$base_dir/reports/hmm_err \
        --wrap="\
            source $HOME/miniconda3/bin/activate hmmer

            hmmscan -E 1e-04 --tblout $base_dir/03_annotated/v2/NCVOGs/$bin.NCVOG_tblout \
            $WORK/bird_bags/db/Needham_2019/47_NCVOGs/47.proteins.megavirales_Needham_2019.hmms \
            $base_dir/02_predicted/v2/$bin &>/dev/null \
            "
done

# again for bins not increasing in abundance
for bin in $(ls $base_dir/02_predicted/v2.2/); do
    sbatch \
        --job-name=$bin.hmm \
        --nodes=1 \
        --tasks-per-node=32 \
        --cpus-per-task=1 \
        --mem=120G \
        --time=00:10:00 \
        --output=$base_dir/reports/hmm_out \
        --error=$base_dir/reports/hmm_err \
        --wrap="\
            source $HOME/miniconda3/bin/activate hmmer

            hmmscan -E 1e-04 --tblout $base_dir/03_annotated/v2.2/NCVOGs/$bin.NCVOG_tblout \
            $WORK/bird_bags/db/Needham_2019/47_NCVOGs/47.proteins.megavirales_Needham_2019.hmms \
            $base_dir/02_predicted/v2.2/$bin &>/dev/null \
            "
done

# again for rest of new polB bins
for bin in $(ls $base_dir/02_predicted/v2_polB_bins/); do
    sbatch \
        --job-name=$bin.hmm \
        --nodes=1 \
        --tasks-per-node=32 \
        --cpus-per-task=1 \
        --mem=120G \
        --time=00:10:00 \
        --output=$base_dir/reports/hmm_out \
        --error=$base_dir/reports/hmm_err \
        --wrap="\
            source $HOME/miniconda3/bin/activate hmmer

            hmmscan -E 1e-04 --tblout $base_dir/03_annotated/v2_polB_bins/NCVOGs/$bin.NCVOG_tblout \
            $WORK/bird_bags/db/Needham_2019/47_NCVOGs/47.proteins.megavirales_Needham_2019.hmms \
            $base_dir/02_predicted/v2_polB_bins/$bin &>/dev/null \
            "
done