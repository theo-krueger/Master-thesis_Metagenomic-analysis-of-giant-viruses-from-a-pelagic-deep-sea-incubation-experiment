import pandas as pd
from Bio import SeqIO


def taxid2lin(annotation_file, database, header=False):

    # print("\t Loading database...")
    dictionary_df = pd.read_csv(database, low_memory=False)

    if header:
        annotation_df = pd.read_csv(annotation_file, sep="\t", index_col=False)
    else:
        col_names = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend",
                     "sstart", "send", "evalue", "bitscore", "stitle", "staxonids", 
                     "sscinames", "sskingdoms", "skingdoms", "sphylums" 
                     ]
        annotation_df = pd.read_csv(annotation_file, sep="\t", index_col=False, names=col_names)

    # print("\t Finding taxonomy...")
    records = []
    counters = {
        "multi_group": 0,
        "faulty_id": 0}
    for idx, row in annotation_df.iterrows():
        # print(idx)

        seq = row["qseqid"]

        taxon_group = row["sskingdoms"]
        if pd.isnull(taxon_group):
            continue
        
        elif ";" not in taxon_group:
            try:
                taxonomy = dictionary_df.loc[dictionary_df["tax_id"] == int(taxon_group), "superkingdom"].values[0]
            except IndexError:
                # print(f"\t\tImproper sskingdom id found: Row nr {idx}")
                taxonomy = "nan"
                counters["faulty_id"] += 1
        else:
            # print(f"\t\tMultiple superkingdoms found. Row nr {idx}")
            taxonomy = "nan"
            counters["multi_group"] += 1
            
            # list_taxids = taxon_group.split(";")
        
        record = {
                "qseqid": seq,
                "origin": taxonomy
            }
        records.append(record)
            
    df = pd.DataFrame.from_records(records)  
    print(f"\tErrors encountered: {counters['multi_group']} multiple IDs, {counters['faulty_id']} faulty IDs. Values replaced with 'nan'")
    return df


def extract_protein_origins(annotation_file, taxon_database, protein_lookup_file, **kwargs):
    print("\tAdding taxon names to id...")
    annotation_table = taxid2lin(annotation_file=annotation_file, database=taxon_database)
    
    print("\tCollecting protein stats...")
    protein_df = pd.read_csv(protein_lookup_file, sep=r"\s+", names=["protein", "contig", "start", "end", "direction"])

    origin_colors = {
        "Bacteria":     "80beaf",
        "Archaea":      "b3ddd1",
        "Eukaryota":    "6b3074",
        "Viruses":      "ee9c6c",
        "Unknown":      "d1dce2"
    }

    records = []
    for idx, row in annotation_table.iterrows():
        qseqid = row["qseqid"]
        origin = row['origin']
        #print(origin)

        if origin == "Bacteria":
            color = origin_colors["Bacteria"]
        elif origin == "Archaea":
            color = origin_colors["Archaea"]
        elif origin == "Eukaryota":
            color = origin_colors["Eukaryota"]
        elif origin == "Viruses":
            color = origin_colors["Viruses"]
        else: # Unknown
            color = origin_colors["Unknown"]

        protein_stats = protein_df.loc[protein_df["protein"] == qseqid].values.flatten().tolist()

        name = protein_stats[1]
        start = protein_stats[2]
        end = protein_stats[3]

        record = {
             "name": name,
             "start": start,
             "end": end,
             "color": f"color={color}"
        }
        records.append(record)

        # print(description_list)
        # print(record)       
        
    df = pd.DataFrame.from_records(records)

    return df


def create_protein_origins_bed(annotation_file, protein_file, taxon_database, out_path):
    origin = extract_protein_origins(
        annotation_file=annotation_file, 
        protein_file=protein_file, 
        taxon_database=taxon_database,
        protein_lookup_file=f"{out_path}/proteins_lookup.txt")
    
    with open(f"{out_path}/protein_origin_tax.txt", 'w') as file:
            dfAsString = origin.to_string(header = False, index = False)
            file.write(dfAsString)


if __name__ == "__main__":
    out = extract_protein_origins(
        protein_file="/Users/tkrueger/Library/CloudStorage/GoogleDrive-theokrueger.marbio@gmail.com/My Drive/AG_Worden/BirdBags/YuChens_pipeline/v2_ESOM/predicted/BB_2016_1040m_12C_4m_19_L002_ESOM_4_1-contigs.fa.faa",
        annotation_file="/Users/tkrueger/Library/CloudStorage/GoogleDrive-theokrueger.marbio@gmail.com/My Drive/AG_Worden/BirdBags/YuChens_pipeline/v2_ESOM/annot/NCBI_NR/BB_2016_1040m_12C_4m_19_L002_ESOM_4_1-contigs.fa.faa.vs.NRwTax.diamond",
        taxon_database="/Users/tkrueger/Desktop/ncbi_lineages_2023-04-04.csv")
    print(out.to_string())
     