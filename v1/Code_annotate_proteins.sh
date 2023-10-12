
# finding NCVOGs
for bin in $(ls 02_predicted); do
sbatch --job-name=$bin.hmm \
--nodes=1 \
--tasks-per-node=32 \
--cpus-per-task=1 \
--mem=120G \
--time=00:10:00 \
--output=reports/hmm_out \
--error=reports/hmm_err \
--wrap="\
hmmscan -E 1e-04 --tblout 03_annotated/NCVOGs/$bin.NCVOG_tblout \
$WORK/bird_bags/db/Needham_2019/47_NCVOGs/47.proteins.megavirales_Needham_2019.hmms \
02_predicted/$bin &>/dev/null"
done

### annotate
# make db
sbatch --job-name=makedb \
--nodes=1 \
--tasks-per-node=4 \
--cpus-per-task=1 \
--mem=120G \
--time=12:00:00 \
--output=reports/makedb_out \
--error=reports/makedb_err \
--mail-type=END,FAIL,TIME_LIMIT \
--mail-user=tkrueger@geomar.de \
--wrap="\
diamond makedb --in $WORK/bird_bags/db/NCBI_vir_20220927/nr.gz -d $WORK/bird_bags/db/NCBI_vir_20220927/nr.gz.dmnd"

#annotate
#only ran for 32 highest probability ones for now
# change name and which files to run, 4 in parallel

#!/bin/sh

#SBATCH --job-name=diamond_4
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --qos=long
#SBATCH --time=5-00:00:00
#SBATCH --output=reports/diamond_4_out
#SBATCH --error=reports/diamond_4_err
#SBATCH --partition=cluster
#SBATCH --mail-type=END,FAIL,TIME-LIMIT
#SBATCH --mail-user=theokrueger.infomail@gmail.com

source $HOME/miniconda3/bin/activate diamond

for bin in $(cat third_filter_above25uniqueNCVOGs.txt | sed -n '16,20p') ; do
diamond blastp \
--query 02_predicted/individual/$bin.fasta.faa \
--db $WORK/bird_bags/db/NCBI_vir_20220927/nr.gz.dmnd \
--out 03_annotated/NR_32_over25NCVOG/$bin.vs.NR.diamond \
--threads 32 \
--evalue 1e-04 \
--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
--max-target-seqs 10
done


