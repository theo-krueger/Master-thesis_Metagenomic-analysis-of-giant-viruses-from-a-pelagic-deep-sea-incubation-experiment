from Bio import SeqIO
import pandas as pd

def extract_protein_positions(file):
     refseq = SeqIO.parse(file, "fasta")

     records = []
     records_full = []
     for sequence in refseq:
          description_list = sequence.description.split(" # ")

          protein = description_list[0]
          name = protein.rpartition("_")[0]
          start = description_list[1]
          end = description_list[2]
          if description_list[3] == "1":
               direction = "forward"
          elif description_list[3] == "-1":
               direction = "reverse"
          else:
               print("BIG ERROR HAPPENING")
          
          record = {
             "name": name,
             "start": start,
             "end": end,
             "direction":direction,
          }
          records.append(record)

          record_full = {
             "protein": protein,
             "name": name,
             "start": start,
             "end": end,
             "direction":direction,
          }
          records_full.append(record_full)

        # print(description_list)
        # print(record)       
     
     df_full = pd.DataFrame.from_records(records_full)
     df = pd.DataFrame.from_records(records)
     df_forward = df.loc[df["direction"] == "forward"]
     df_reverse = df.loc[df["direction"] == "reverse"]

     return df_full, df, df_forward, df_reverse


def create_proteins_bed(file, out_path):
    prot_full, prot, prot_f, prot_r = extract_protein_positions(file)
    with open(f"{out_path}/proteins_lookup.txt", 'w') as file:
            dfAsString = prot_full.to_string(header = False, index = False)
            file.write(dfAsString)
    with open(f"{out_path}/proteins.txt", 'w') as file:
            dfAsString = prot.to_string(header = False, index = False)
            file.write(dfAsString)
    with open(f"{out_path}/proteins_forward.txt", 'w') as file:
            dfAsString = prot_f.to_string(header = False, index = False)
            file.write(dfAsString)
    with open(f"{out_path}/proteins_reverse.txt", 'w') as file:
            dfAsString = prot_r.to_string(header = False, index = False)
            file.write(dfAsString)

if __name__ == "__main__":
     extract_protein_positions("/Users/tkrueger/Library/CloudStorage/GoogleDrive-theokrueger.marbio@gmail.com/My Drive/AG_Worden/BirdBags/YuChens_pipeline/v2_ESOM/predicted/BB_2016_1040m_12C_4m_19_L002_ESOM_4_1-contigs.fa.faa")