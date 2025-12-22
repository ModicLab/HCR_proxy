from localcider.sequenceParameters import SequenceParameters
import pandas as pd
import multiprocessing as mp

all_nuclear = pd.read_csv('data/curated_data_for_modeling/localcider_features.tsv', sep='\t')
all_nuclear = all_nuclear[["Protein ID", "Sequence"]]

def run_sequence_analysis(sequence, num_scrambles=100000, random_seed=None):
    try:
        print(f"Running analysis for {sequence}")
        SeqObj = SequenceParameters(sequence)
        SeqObj.save_zscoresAndPlots(num_scrambles=num_scrambles, random_seed=random_seed)
        return "Success"
    except Exception as e:
        return f"Error: {e}"

def run_on_dataframe(df, num_scrambles=100000, random_seed=None, n_jobs=None):
    
    if n_jobs is None:
        n_jobs = max(mp.cpu_count(), len(df))

    sequences = df['Sequence'].tolist()

    print(f"Number of sequences to process: {len(sequences)}")
    print(n_jobs)
    
    with mp.Pool(processes=n_jobs) as pool:
        pool.starmap(
            run_sequence_analysis,
            [(seq, num_scrambles, random_seed) for seq in sequences]
        )

run_on_dataframe(all_nuclear, num_scrambles=100000, random_seed=42, n_jobs=None)
