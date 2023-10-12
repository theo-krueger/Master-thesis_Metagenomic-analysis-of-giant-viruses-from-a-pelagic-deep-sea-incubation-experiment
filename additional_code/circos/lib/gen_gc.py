from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import pandas as pd

def sliding_gc(fasta_file, window_size=1000, steps=None):
    if steps is None:
        steps = window_size

    refseq = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    records = []
    gc_extremes = []
    for seq_name, seq in refseq.items():
        seq_length = len(seq)

        mod_window_size = min(seq_length, window_size)
        starts = list(range(0, seq_length - mod_window_size + 1, steps))

        max_gc = 0
        min_gc = 100
        
        for chunk_start in starts:
            chunk_end = chunk_start + mod_window_size
            chunk = seq[chunk_start:chunk_end]
            chunk_GC = gc_fraction(chunk)

            record = {
                "name": seq_name,
                "start": chunk_start + 1,
                "end": chunk_end,
                "gc": chunk_GC
            }
            records.append(record)

            if chunk_GC > max_gc:
                max_gc = chunk_GC
            if chunk_GC < min_gc:
                min_gc = chunk_GC 

        gc_extreme = {
            "name": seq_name,
            "min_gc": min_gc,
            "max_gc": max_gc
        }   
        gc_extremes.append(gc_extreme) 

    df = pd.DataFrame.from_records(records)
    stats = pd.DataFrame.from_records(gc_extremes)

    min_gc_all = pd.to_numeric(stats["min_gc"]).min(axis=0).astype(str)
    max_gc_all =  pd.to_numeric(stats["max_gc"]).max(axis=0).astype(str)
    gc_extreme_all = {
        "name": "all_contigs",
        "min_gc": min_gc_all,
        "max_gc": max_gc_all
    }
    gc_extreme_all_df = pd.DataFrame.from_records(gc_extreme_all, index=[0])
    stats_final = pd.concat([stats, gc_extreme_all_df], ignore_index = True)

    return df, stats_final

def create_gc_bed(file, out_path):
    gc, stats = sliding_gc(file)
    with open(f"{out_path}/gc.txt", 'w') as file:
        dfAsString = gc.to_string(header = False, index = False)
        file.write(dfAsString)
    with open(f"{out_path}/gc_stats.txt", 'w') as file:
        dfAsString = stats.to_string(header = False, index = False)
        file.write(dfAsString)


if __name__ == "__main__":
    gc, stats = sliding_gc("/Users/tkrueger/Library/CloudStorage/GoogleDrive-theokrueger.marbio@gmail.com/My Drive/AG_Worden/BirdBags/YuChens_pipeline/v2_ESOM/DNA_seq/BB_2016_1040m_12C_4m_19_L002_ESOM_4_1-contigs.fa")
    print(stats)