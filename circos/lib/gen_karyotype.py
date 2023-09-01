from Bio import SeqIO
import pandas as pd

def base_stats(fasta_file):
    refseq = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    list_base = []
    for seq_id, (seq_name, seq) in enumerate(refseq.items()):
        seq_length = len(seq)

        record = {
            "chr": "chr",
            "spacer": "-",
            "id": seq_name,
            "label": seq_id + 1,
            "start": 1,
            "end": seq_length,
            "color": "vdgey"

        }
        list_base.append(record)

    df_base = pd.DataFrame.from_records(list_base)
    return df_base

def create_karyotype_bed(file, out_path):
    kt = base_stats(file)
    with open(f"{out_path}/karyotype.txt", 'w') as file:
            dfAsString = kt.to_string(header = False, index = False)
            file.write(dfAsString)