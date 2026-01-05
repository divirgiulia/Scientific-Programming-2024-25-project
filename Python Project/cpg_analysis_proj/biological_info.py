import csv
import os

Gene_annotation_file = os.path.join("test_data", "gene_annotations.csv")
Promoter_annotation_file = os.path.join("test_data", "promoter_annotations.csv")

def load_annotations(file_path):
    annotations = []
    with open(file_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            annotations.append({
                "name": row.get("name", ""),
                "start": int(row["start"]),
                "end": int(row["end"])
            })
    return annotations

if not os.path.exists(Gene_annotation_file):
    print(f"Warning: Gene annotation file {Gene_annotation_file} not found.")
    GENE_ANNOTATIONS = []

if not os.path.exists(Promoter_annotation_file):
    print(f"Warning: Promoter annotation file {Promoter_annotation_file} not found.")
    PROMOTER_ANNOTATIONS = []

Gene_annotations = load_annotations(Gene_annotation_file) if os.path.exists(Gene_annotation_file) else []
if not Gene_annotations:
    print(f"Warning: Gene annotation file {Gene_annotation_file} not found.")

Promoter_annotations = load_annotations(Promoter_annotation_file) if os.path.exists(Promoter_annotation_file) else []
if not Promoter_annotations:
    print(f"Warning: Promoter annotation file {Promoter_annotation_file} not found.")


def find_gene(start, end):
    for gene in Gene_annotations:
        if start <= gene["end"] and end >= gene["start"]:
            return gene["name"]
    return None

def is_promoter(start, end):
    for promoter in Promoter_annotations:
        if start <= promoter["end"] and end >= promoter["start"]:
            return True
    return False

def annotate_islands(islands):
    annotated = []
    for island in islands:
        start, end = island["start"], island["end"]
        annotated.append({
            "start": start,
            "end": end,
            "cpg_percent": island.get("gpg_percent"),
            "gene": find_gene(start, end),
            "promoter": is_promoter(start, end)    
        })
    return annotated

if __name__ == "__main__":
    test_islands = [(90, 150), (810, 900), (2100, 2300)]
    results = annotate_islands(test_islands)
    for r in results:
        print(r)
