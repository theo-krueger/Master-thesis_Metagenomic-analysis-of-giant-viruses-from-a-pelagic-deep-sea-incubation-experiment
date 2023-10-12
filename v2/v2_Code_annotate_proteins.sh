
base_dir=$WORK/bird_bags/YuChens_pipeline

############## diamond ##############

#only ran for 32 highest probability ones for now
# change name and which files to run, 4 in parallel

# now works for all of them too if using fewer jobs

# make diamond db with taxonomy
sbatch --job-name=makedb \
    --nodes=1 \
    --tasks-per-node=4 \
    --cpus-per-task=1 \
    --mem=120G \
    --time=24:00:00 \
    --output=reports/makedb_out \
    --error=reports/makedb_err \
    --mail-type=END,FAIL,TIME_LIMIT \
    --mail-user=theokrueger.infomail@geomar.de \
    --wrap="\
        diamond makedb \
        --in $WORK/bird_bags/db/NCBI_vir_20220927/nr.gz \
        --db $WORK/bird_bags/db/NCBI_vir_20220927_WithTaxonomy/nr.gz_WTax.dmnd \
        --taxonmap $WORK/bird_bags/db/NCBI_vir_20220927_WithTaxonomy/prot.accession2taxid.FULL.gz \
        --taxonnodes $WORK/bird_bags/db/NCBI_vir_20220927_WithTaxonomy/nodes.dmp \
        --taxonnames $WORK/bird_bags/db/NCBI_vir_20220927_WithTaxonomy/names.dmp"

# also installed tax2lin

#!/bin/sh

#SBATCH --job-name=diamond
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --qos=long
#SBATCH --time=5-00:00:00
#SBATCH --output=reports/diamond_out
#SBATCH --error=reports/diamond_err
#SBATCH --partition=cluster
#SBATCH --mail-type=END,FAIL,TIME-LIMIT
#SBATCH --mail-user=theokrueger.infomail@gmail.com

source $HOME/miniconda3/bin/activate diamond-2.1.4

for bin in $(ls 02_predicted/v2/); do
    diamond blastp \
        --query 02_predicted/v2/$bin \
        --db $WORK/bird_bags/db/NCBI_vir_20220927_WithTaxonomy/nr.gz_WTax.dmnd \
        --out 03_annotated/v2/NR/$bin.vs.NR.diamond \
        --threads 32 \
        --evalue 1e-04 \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle staxids sscinames \
        --max-target-seqs 10
done

# add taxonomy

python3 $HOME/tools/taxid2lin.py \
-i $WORK/bird_bags/YuChens_pipeline/03_annotated/v2/NR/ \
-o $WORK/bird_bags/YuChens/03_annotated/v2/NR_WNames/ \
-d $WORK/bird_bags/db/NCBI_vir_20220927_WithTaxonomy/ncbi_lineages_2023-04-04.csv

# join for easy searching
for file in $(ls 03_annotated/v2/NR/); do
    bin=$(echo $file | cut -f1 -d '-')
    awk -v bin=$bin '{print bin, $0}' 03_annotated/v2/NR/$file >>03_annotated/v2/NR_all.tbl
done
# make into library

# interesting stuff, do not consider names of origin species
awk -F '[\t[]' '$13 ~ /cyanide|cyanase|thiocyanate|arsen|mercur|quorum|methan|quinoprotein|copper|secretion|rhodanese|chito|chiti|phenol|benzene|benzo/ {print $0}' \
    $base_dir/03_annotated/v2/NR_all.tbl \
    > $base_dir/03_annotated/v2/interesting_proteins.txt

## only results for very likely NCLDV
for bin in $(cat $base_dir/00_stats/v2/v2_bins_afterfilter_above_25_uniqNCVOG.txt); do
    grep "$bin" $base_dir/03_annotated/v2/interesting_proteins.txt \
        >> $base_dir/03_annotated/v2/interesting_proteins_32bins.txt
done
# make bin-contigs file
for dir in $(ls $base_dir/06_virsorter2/v2/02_output/); do
    awk '$8 == "NCLDV" {print $1}' $base_dir/06_virsorter2/v2/02_output/$dir/final-viral-score.tsv | \
        cut -f1 -d '|' \
        >> $base_dir/06_virsorter2/v2/virsorter2_NCLDV_contigs
done
cat $base_dir/06_virsorter2/v2/virsorter2_NCLDV_contigs | sort | uniq > $base_dir/06_virsorter2/v2/virsorter2_NCLDV_contigs_uniq
# choose only proteins on NCLDV contigs
for contig in $(cat $base_dir/06_virsorter2/v2/virsorter2_NCLDV_contigs_uniq); do
    grep "$contig" $base_dir/03_annotated/v2/interesting_proteins_32bins.txt \
        >> $base_dir/03_annotated/v2/interesting_proteins_32bins_onlyNCLDVcontigs.txt
done
sort \
    -k1 \
    -o $base_dir/03_annotated/v2/interesting_proteins_32bins_onlyNCLDVcontigs.txt \
    $base_dir/03_annotated/v2/interesting_proteins_32bins_onlyNCLDVcontigs.txt
# only ESOM
for contig in $(cat $base_dir/06_virsorter2/v2/virsorter2_NCLDV_contigs_uniq); do
    grep "$contig" $base_dir/03_annotated/v2/NR_all.tbl >> $base_dir/03_annotated/v2/NR_NCLDVcontig.txt
done
awk '/ESOM/ {print}' $base_dir/03_annotated/v2/NR_NCLDVcontig.txt > temp && mv temp $base_dir/03_annotated/v2/NR_NCLDVcontig.txt
sort \
    -k1 \
    -o $base_dir/03_annotated/v2/NR_NCLDVcontig.txt \
    $base_dir/03_annotated/v2/NR_NCLDVcontig.txt

# add lineage info with custom py script
sbatch \
    --job-name=dmnd_addlin \
    --nodes=1 \
    --tasks-per-node=1 \
    --cpus-per-task=1 \
    --mem=120G \
    --time=2:00:00 \
    --output=$WORK/bird_bags/YuChens_pipeline/reports/dmnd_addlin_out \
    --error=$WORK/bird_bags/YuChens_pipeline/reports/dmnd_addlin_err \
    --wrap="\
        python3.9 $HOME/tools/taxid2lin.py \
        -i $WORK/bird_bags/YuChens_pipeline/03_annotated/v2/NR/ \
        -o $WORK/bird_bags/YuChens_pipeline/03_annotated/v2/NR_WNames/ \
        -d $WORK/bird_bags/db/NCBI_vir_20220927_WithTaxonomy/ncbi_lineages_2023-04-04.csv"

# combine
for file in $(ls 03_annotated/v2/NR/); do
    bin=$(echo $file | cut -f1 -d '.')
    tail -n +4 03_annotated/v2/NCVOGs/$bin.fa.faa.NCVOG_tblout | head -n -10 | awk '{split($1,a,"[-]"); print $3 "\t" $1 "\t" $5 "\t" $6}' >03_annotated/v2/$bin.hmm.temp
    join -j 1 -a 1 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,2.2,2.3,2.4 \
        <(sort -k1 03_annotated/v2/NR/$file) <(sort -k1 03_annotated/v2/$bin.hmm.temp) >03_annotated/v2/combined/$bin.NCVOG.diamond
    rm 03_annotated/v2/$bin.hmm.temp
done

# add header
for file in $(ls 03_annotated/v2/combined/); do
    echo -e 'qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstitle\tstaxonids\tsscinames\tNCVOG\tNCVOG_evalue\tNCVOG_score' |
        cat - 03_annotated/v2/combined/$file >temp && mv temp 03_annotated/v2/combined/$file-wtitle.tbl
done
rm 03_annotated/v2/combined/*.diamond



