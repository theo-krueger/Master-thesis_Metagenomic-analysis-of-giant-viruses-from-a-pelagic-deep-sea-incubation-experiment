from Bio import SeqIO
import pandas as pd
import re

def crispr_positions(crispr_file):
    try:
        seqs = SeqIO.parse(crispr_file, "fasta")

        crisprs = []
        for seq in seqs:
            seq_name = seq.description
            name = seq_name.split(";")[0]
            numbers = re.findall(r'\d+', seq_name.split(";")[1])
            
            record = {
                "name": name,
                "start": numbers[1],
                "end": numbers[2],
                "value": "value"
            }

            crisprs.append(record)
        
        df_crispr = pd.DataFrame.from_records(crisprs)

    except FileNotFoundError:
         print("\tNo CRISPRs found")
         df_crispr = pd.DataFrame()

    return df_crispr

def create_crispr_bed(file, out_path):
    df = crispr_positions(file)
    with open(f"{out_path}/crispr.txt", 'w') as file:
            if not df.empty:
                dfAsString = df.to_string(header = False, index = False)
                file.write(dfAsString)
            else:
                 pass

if __name__ == "__main__":
    crispr_positions("/Users/tkrueger/Library/CloudStorage/GoogleDrive-theokrueger.marbio@gmail.com/My Drive/AG_Worden/BirdBags/YuChens_pipeline/v2_ESOM/crispr/BB_2016_1040m_12C_4m_19_L002_ESOM_4_2-contigs_CRISPRS.fasta")