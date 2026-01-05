from sequence_utils import global_cpg, sliding_window_cpg, detect_cpg_islands
from biological_info import annotate_islands

def run_analysis(task, sequence, window=None, step=None, threshold=None):
    task = task.lower()

    #  global CpG content
    if task == 'global':
        return {"cpg_percentage": global_cpg(sequence)}

    # sliding window content
    elif task == 'sliding':
        if window is None or step is None:
            raise ValueError("Window size and step size must be provided for sliding window analysis.")
        window = int(window)
        step = int(step)        
        window_results = sliding_window_cpg(sequence, window, step)
        return {"window": window, "step": step, "sliding_results": window_results}
    
    # CpG island detection
    elif task == 'island':
        if threshold is None or window is None or step is None:
            raise ValueError("Window size, step size and threshold must be provided for island detection.")
        
        step = int(step)
        window = int(window)
        threshold = float(threshold)
        islands = detect_cpg_islands(sequence, threshold=threshold, window=window, step=step)
        annotated_islands = annotate_islands(islands)
        return {"window": window, "step": step ,"threshold": threshold,"cpg_islands": annotated_islands}
    else:
        raise ValueError("Unknown task")

