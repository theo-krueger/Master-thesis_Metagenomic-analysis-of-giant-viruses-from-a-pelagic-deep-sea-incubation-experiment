import pandas as pd

def mean_gc(file_path):
    gc_df = pd.read_csv(f"{file_path}/gc.txt", sep=r"\s+", names=["name", "start", "end", "gc"])
    gc_mean = gc_df["gc"].mean()

    return gc_mean

def count_proteins(file_path):
    n_proteins = 0
    with open(f"{file_path}/proteins_lookup.txt") as file:
        for line in file:
            if line.strip():
                n_proteins += 1
    
    return n_proteins
    
def count_ncvogs(file_path):
    n_ncvogs = 0
    with open(f"{file_path}/ncvogs_list.txt") as file:
        ncvog_list = []
        for line in file:
            if line.strip():
                ncvog = line.split()[0]
                ncvog_list.append(ncvog)
                n_ncvogs += 1
        n_ncvogs_uniq = len(set(ncvog_list))
    
    return n_ncvogs, n_ncvogs_uniq

def generate_label(name, out_files_path):
    
    n_proteins = count_proteins(out_files_path)
    n_ncvogs, n_ncvogs_uniq = count_ncvogs(out_files_path)
    avg_gc = mean_gc(out_files_path) * 100


    label_text = f"\
        Name: {name}\n\
        Mean GC%: {avg_gc}\n\
        Proteins: {n_proteins}\n\
        NCVOGs: {n_ncvogs} ({n_ncvogs_uniq})\
    "
    return label_text


def create_label_file(name, out_path):
    label_text = generate_label(name, out_path)

    with open(f"{out_path}/label.txt", "w") as file:
        file.write(label_text)
