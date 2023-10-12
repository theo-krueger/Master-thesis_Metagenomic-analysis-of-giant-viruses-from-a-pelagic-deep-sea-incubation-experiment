import pandas as pd

def read_NCVOG_table(annotation_file):
    clean_file = []
    with open(annotation_file, "r") as file:
        for line in file.readlines():
            if not line.startswith("#"):
                line_content = line.split()
                clean_file.append(line_content)
    # print(clean_file)
    df_column_headers = [
        "target name", "accession", "query name", "accession", "E-value", "score", 
        "bias", "E-value", "score", "bias", "exp", "reg", "clu", "ov", "env", "dom", 
        "rep", "inc", "description of target"
                         ]
    annotation_df = pd.DataFrame(clean_file, columns=df_column_headers)
    
    return annotation_df


def extract_NCVOG_positions(annotation_file, protein_lookup_file):
    annotation_df = read_NCVOG_table(annotation_file)
    protein_df = pd.read_csv(protein_lookup_file, sep=r"\s+", names=["protein", "name", "start", "end", "direction"])

    records = []
    ncvogs = []
    for idx, row in annotation_df.iterrows():
        protein = row['query name']
        ncvog_name = row['target name']

        protein_stats = protein_df.loc[protein_df["protein"] == protein].values.flatten().tolist()
        protein_name = protein_stats[1]
        start = protein_stats[2]
        end = protein_stats[3]

        record = {
                "name": protein_name,
                "start": start,
                "end": end,
                "value": "value"
            }
        records.append(record)

        ncvog = {
            "ncvog": ncvog_name,
            "name": protein_name,
            "start": start,
            "end": end,
            }
        ncvogs.append(ncvog)       

    df_ncvogs = pd.DataFrame.from_records(records)
    df_ncvogs_full = pd.DataFrame.from_records(ncvogs)
    
    return df_ncvogs, df_ncvogs_full


def create_NCVOG_bed(file, out_path):
    df, df_full = extract_NCVOG_positions(annotation_file=file,
                                 protein_lookup_file=f"{out_path}/proteins_lookup.txt")
    with open(f"{out_path}/ncvogs.txt", 'w') as file:
        dfAsString = df.to_string(header = False, index = False)
        file.write(dfAsString)
    with open(f"{out_path}/ncvogs_list.txt", 'w') as file:
        dfAsString = df_full.to_string(header = False, index = False)
        file.write(dfAsString)


if __name__ == "__main__":
    read_NCVOG_table("/Users/tkrueger/Library/CloudStorage/GoogleDrive-theokrueger.marbio@gmail.com/My Drive/AG_Worden/BirdBags/YuChens_pipeline/v2_ESOM/NCVOGs/hmmscan/BB_2016_1040m_12C_4m_19_L002_ESOM_4_1-contigs.fa.faa.NCVOG_tblout")
    