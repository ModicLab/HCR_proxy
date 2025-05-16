import os
import torch
import numpy as np
import pandas as pd
import esm
import deepspeed

configs = [
    ("esm2_t33_650M_UR50D", 33,   "650M"),
    ("esm2_t36_3B_UR50D",  36,   "3B"),
    ("esm2_t48_15B_UR50D", 48,   "15B"),
]

data_csv   = "data/for_training/sequence_and_groups_for_esm2_15B_representations.tsv"
out_dir    = "data/curated_data_for_modeling/esm2_representations"

df   = pd.read_csv(data_csv, sep="\t")
data = list(zip(df["Protein ID"], df["Sequence"]))

device = torch.device("cuda", int(os.environ.get("LOCAL_RANK", 0)))

for fn_name, layer, suffix in configs:
    print(f"\n--- Running {fn_name} (repr_layers=[{layer}]) ---")
    
    model_fn = getattr(esm.pretrained, fn_name)
    model, vocab = model_fn()
    batch_converter = vocab.get_batch_converter()
    model.eval()

    model = deepspeed.init_inference(
        model,
        mp_size=4,
        dtype=torch.float16,
        replace_method='auto',
        replace_with_kernel_inject=True
    )

    all_reps = []
    for pid, seq in data:
        labels, strs, tokens = batch_converter([(pid, seq)])
        tokens = tokens.to(device)
        with torch.no_grad():
            out = model(tokens, repr_layers=[layer], return_contacts=False)
        reps = out["representations"][layer]     
        mean_rep = reps.cpu().numpy().mean(axis=1)  
        all_reps.append(mean_rep)
    all_reps = np.concatenate(all_reps, axis=0) 

    out_path = os.path.join(out_dir, f"representations_45S_47S_esm2_{suffix}.npy")
    np.save(out_path, all_reps)
    print(f"Saved {all_reps.shape} array to {out_path}")
