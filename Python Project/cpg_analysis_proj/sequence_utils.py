
def read_sequence(seq_or_path: str) -> str:
    try:
        if seq_or_path.endswith((".fasta", ".fa")): # FASTA file path
            with open(seq_or_path, "r") as f:
                lines = f.readlines()
                return "".join(line.strip() for line in lines if not line.startswith(">"))
        elif seq_or_path.strip().startswith(">"): # FASTA formatted string
            lines = seq_or_path.strip().splitlines()
            return "".join(line.strip() for line in lines if not line.startswith(">"))
        else: # just raw DNA seq
            return seq_or_path.strip()

    except FileNotFoundError:
        raise ValueError("File not found.")

def is_valid_dna(seq):
    return all(base in "ACTGactg" for base in seq)

def global_cpg(seq: str) -> float:
    seq = seq.upper()
    if len(seq) < 2:
        return 0.0
    count = sum(1 for i in range(len(seq)-1) if seq[i:i+2] == 'CG')
    return round((count / (len(seq) -1)) * 100, 2) if len(seq) > 1 else 0

def sliding_window_cpg(seq: str, window: int, step: int):
    seq = seq.upper()
    results = []

    for i in range(0, len(seq) - window + 1, step):
        window_seq = seq[i:i+window]
        cpg_count = sum(1 for j in range(window - 1) if window_seq[j:j+2] == 'CG')
        percent = (cpg_count/(window - 1)) * 100 if window > 1 else 0
        results.append({'position':i, 
                        'cpg_percent': round(percent,2)})
    return results

def detect_cpg_islands(seq, threshold=60.0, window=200, step=50):
    seq = seq.upper()
    islands = []

    for i in range(0, len(seq) - window +1, step):
        window_seq = seq[i:i+window]
        cpg_count = sum(1 for j in range(window - 1) if window_seq[j:j+2] == 'CG')
        cpg_percent = (cpg_count/(window -1)) * 100 if window > 1 else 0

        if cpg_percent > threshold:
            islands.append({
                'start':i,
                'end':i+window,
                'cpg_percent': round(cpg_percent, 2)
            })
    return islands

