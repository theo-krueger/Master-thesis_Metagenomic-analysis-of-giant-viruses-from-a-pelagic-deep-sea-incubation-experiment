import sys

sys.dont_write_bytecode = True

import argparse
import shutil
from pathlib import Path
import os
from contextlib import contextmanager

from lib import *


def get_args(argv=None):
    argParser = argparse.ArgumentParser()

    argParser.add_argument("-ig", "--genome-file", help="FASTA genome file")
    argParser.add_argument("-ip", "--predicted-proteins-file", help="FASTA predicted proteins file")
    argParser.add_argument("-ia", "--annotation-file", help="Diamond blast annotation file")
    argParser.add_argument("-itdb", "--taxon-database", help="Taxon database for translation of annotated taxon IDs")
    argParser.add_argument("-ic", "--crispr-file", help="cmt-mod crispr annotation fasta file")
    argParser.add_argument("-in", "--ncvog-file", help="hmmscan NCVOG annotation table")
    argParser.add_argument("-n", "--name", help="genome name")
    argParser.add_argument("-o", "--output-path", help="Output path for plots", required=True)

    return argParser.parse_args(argv)
    

@contextmanager
def cd(new_dir):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(new_dir))
    try:
        yield
    except:
        os.chdir(prevdir)


def create_paths(args):
    genome_name = args.genome_file.split("/")[-1].split(".")[0]

    base_path = f"{args.output_path}/{genome_name}"
    data_path = f"{base_path}/data"
    etc_path = f"{base_path}/config"


    Path(data_path).mkdir(parents=True, exist_ok=True)
    Path(etc_path).mkdir(exist_ok=True)

    return base_path, data_path, etc_path


def create_data_files(data_path, args):
    if not os.path.isfile(f"{data_path}/karyotype.txt"):
        print("Creating karyotype file...")
        gen_karyotype.create_karyotype_bed(file = args.genome_file, out_path = data_path)
    if not os.path.isfile(f"{data_path}/gc.txt"):   
        print("Creating gc file...")
        gen_gc.create_gc_bed(file = args.genome_file, out_path = data_path)
    if not os.path.isfile(f"{data_path}/proteins.txt"):   
        print("Creating proteins file...")
        gen_proteins.create_proteins_bed(file = args.predicted_proteins_file, out_path = data_path)
    if not os.path.isfile(f"{data_path}/protein_origin_tax.txt"):   
        print("Creating protein origins file:")
        gen_protein_origins.create_protein_origins_bed(annotation_file = args.annotation_file, 
                                                       protein_file = args.predicted_proteins_file,
                                                       taxon_database = args.taxon_database,
                                                       out_path = data_path)
    if not os.path.isfile(f"{data_path}/crispr.txt"):   
        print("Creating crispr file...")
        gen_crispr.create_crispr_bed(file = args.crispr_file, out_path = data_path)
    if not os.path.isfile(f"{data_path}/ncvogs.txt"):   
        print("Creating NCVOG file...")
        gen_NCVOGs.create_NCVOG_bed(file = args.ncvog_file, out_path = data_path)

def create_config(data_path, config_path):
    if not os.path.isfile(f"{config_path}/circos_final.conf"): 
        gen_config.create_config_file(data_path, config_path)

def copy_config(config_path):
    print("Checking and adding config files...")
    
    current_dir = os.path.dirname(__file__)
    all_files = os.listdir(f"{current_dir}/config/")
    for file in all_files:
        if not os.path.isfile(f"{config_path}/{file}"):
            shutil.copy(f"{current_dir}/config/{file}", config_path)


def create_label(data_path, args):
    if not os.path.isfile(f"{data_path}/label.txt"):
        print("Creating label file...")
        gen_label.create_label_file(args.name, data_path)

def create_circos(target_path, **args):
    print(f"Creating plot at: {target_path}")
    with cd(target_path):
        os.system(f"conda run -n circos circos -conf config/circos.conf")

def main(args):
    
    print("\n CUSTOM CIRCOS")
    print(f"\nCreating for: {args.name}")

    # paths
    base_path, data_path, etc_path = create_paths(args=args)
    # data
    create_data_files(data_path=data_path, args=args)
    # config
    copy_config(config_path=etc_path)
    # label
    create_label(data_path=data_path, args=args)
    # plot
    create_circos(target_path=base_path)



if __name__ == "__main__":
    args = get_args()
    main(args=args)


# python3 $WORK/bird_bags/YuChens_pipeline/circos/circos_terminal/circos_main.py -ig $WORK/bird_bags/YuChens_pipeline/01_DNA_sequences/v2_5MAGs/BB_2016_1040m_12C_4m_19_L002_ESOM_4_1-contigs.fa -ip $WORK/bird_bags/YuChens_pipeline/02_predicted/v2_polB_bins/BB_2016_1040m_12C_4m_19_L002_ESOM_4_1-contigs.fa.faa -ia $WORK/bird_bags/YuChens_pipeline/03_annotated/v2_polB_bins/NR/BB_2016_1040m_12C_4m_19_L002_ESOM_4_1-contigs.fa.faa.vs.NRwTax.diamond -itdb $WORK/bird_bags/db/NCBI_vir_20220927_WithTaxonomy/ncbi_lineages_2023-04-04.csv -ic $WORK/bird_bags/YuChens_pipeline/03_annotated/v2_polB_bins/crispr/BB_2016_1040m_12C_4m_19_L002_ESOM_4_1-contigs_CRISPRS.fasta -in $WORK/bird_bags/YuChens_pipeline/03_annotated/v2_polB_bins/NCVOGs/BB_2016_1040m_12C_4m_19_L002_ESOM_4_1-contigs.fa.faa.NCVOG_tblout -n BB_4_1 -o $WORK/bird_bags/YuChens_pipeline/circos/
