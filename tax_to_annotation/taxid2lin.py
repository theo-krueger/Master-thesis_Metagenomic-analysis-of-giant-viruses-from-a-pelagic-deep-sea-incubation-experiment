import pandas as pd
from pathlib import Path
import numpy as np
import os
import warnings
import argparse
from colorama import Fore


argParser = argparse.ArgumentParser()
required = argParser.add_argument_group('required arguments')
optional = argParser.add_argument_group('optional arguments')

required.add_argument("-i", "--input-path", help="Input path to the diamond annotation files, folder or single file",
                      required=True)
required.add_argument("-o", "--output-dir", help="Output directory to save the output file to", required=True)
required.add_argument("-d", "--db-path", help="Path to the ncbi_lineages file as .csv", required=True)

optional.add_argument("--header", help="Input files have headers. Default: False", default=False)
optional.add_argument("--col-names", help="List of column names of the input files", nargs='+',
                      default=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend",
                               "sstart", "send", "evalue", "bitscore", "stitle", "staxonids", "sscinames"])
optional.add_argument("--input-extension", help="Specifiying the input file extension. Default: diamond",
                      default="diamond")

args = argParser.parse_args()


def fxn():
    warnings.warn("pending deprecation", PendingDeprecationWarning)
    warnings.warn("future", FutureWarning)


with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    fxn()

    path_input_annot = args.input_path
    path_input_annot = Path(path_input_annot)
    if os.path.isdir(path_input_annot):
        files_input_list = [fn for fn in os.listdir(path_input_annot) if fn.endswith("." + args.input_extension)]
        print(f"Number of files found: {len(files_input_list)}.")
    elif os.path.isfile(path_input_annot):
        files_input_list = [str(path_input_annot)]
    else:
        print("Input path not found.")

    path_output_annot = args.output_dir

    path_input_dic = args.db_path
    path_input_dic = Path(path_input_dic)
    print("Loading database...")
    dictionary_df = pd.read_csv(path_input_dic)

    lineage_columns = list(dictionary_df.columns)[1:8]

    print("Starting annotation.")
    for filename_annotation in files_input_list:

        print(Fore.LIGHTYELLOW_EX + f'Processing: {filename_annotation}' + Fore.RESET)

        if args.header:
            annotation_df = pd.read_csv(path_input_annot/filename_annotation, sep="\t")
        else:
            annotation_df = pd.read_csv(path_input_annot / filename_annotation, sep="\t", names=args.col_names)

        annotation_df[lineage_columns] = np.nan

        annotation_df_new = annotation_df.copy(deep=True)

        lineage_columns_idx = [annotation_df_new.columns.get_loc(a) for a in lineage_columns]
        annotation_df_new[lineage_columns] = annotation_df[lineage_columns].astype('object')

        for idx, row in annotation_df.iterrows():
            # print(idx)
            if pd.isnull(row["staxonids"]):
                print(f'ROW {idx}: No taxon id.')
                continue

            elif ";" not in row["staxonids"]:
                taxid_ann = row["staxonids"]
                names_dic = dictionary_df.loc[dictionary_df["tax_id"] == int(taxid_ann), lineage_columns].values.flatten().tolist()

                if len(names_dic):
                    annotation_df_new.loc[annotation_df_new.index[idx], lineage_columns] = names_dic
                else:
                    print(f'ROW {idx}: No lineage found for taxon id {taxid_ann}.')
                    continue
            else:
                list_taxids = row["staxonids"].split(";")
                # print(f'ROW {idx}: Multiple taxon ids found in line: {list_taxids}.')

                taxid_dic = [[] for i in range(len(lineage_columns))]

                for taxid in list_taxids:
                    names_dic = dictionary_df.loc[
                        dictionary_df["tax_id"] == int(taxid), lineage_columns].values.flatten().tolist()
                    for index, values in enumerate(names_dic):
                        taxid_dic[index].append(names_dic[index])

                new_taxid_dic = {}
                for ind, name in enumerate(lineage_columns):
                    new_taxid_dic[name] = taxid_dic[ind]

                if len(taxid_dic):
                    for name in new_taxid_dic:
                        annotation_df_new.at[annotation_df_new.index[idx], name] = new_taxid_dic[name]

                else:
                    print(f"ROW {idx}: No lineages found.")

        filename_output = os.path.basename(os.path.normpath(filename_annotation)) + "_names.csv"
        filename_output = Path(filename_output)
        annotation_df_new.to_csv(path_output_annot/filename_output, index=False)
        print(Fore.CYAN + f"File saved as {filename_output}" + Fore.RESET)
