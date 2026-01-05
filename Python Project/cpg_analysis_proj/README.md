# CpG Analysis Client-Server Project

## Short description

This project implements a client-server application in Python to analyze the CpG content in DNA sequences. 
The server accepts DNA sequences as input (raw or FASTA format) and performs the following analyses:
- **Global CpG content:** Percentage of continuout CG nucleotides in the entire sequence.
- **Sliding window CpG profile:** Calculates CpG percentage in windows sliding over the sequence.
- **CpG island detection:** Identifies regions with CpG content above a given threshold and annotates them.
Whenever CpG islands are detected, they can be annotated with biological information.
The client sends requests to the server with analyss parameters and receives JSON-formatted results.

---

## Features

- Supports raw DNA and FASTA inputs;
- Configurable server host and port via 'config.ini';
- Efficient handling of multiple analysis tasks;
- Biological annotation of detected CpG islands;
- Simple command-line client interface.

---

## Usage
1. Start the server
```bash
python server.py --host 127.0.0.1 --port 5050
```
2. Use the client
```bash
python client.py --host 127.0.0.1 --port 5050 --task global --input ATCGCGCG
```
Other examples:
```bash
python client.py --host 127.0.0.1 --port 5050 --task sliding --input ATCGCGCG --window 4 --step 2
python client.py --host 127.0.0.1 --port 5050 --task island --input ATCGCGCG --window 6 --step 2 --threshold 60
```

3. Running tests
Unit tests are located in the tests/ folder. Run tests using:
```bash
python -m unittest discover tests/
```

---

## Requirements

- Python 3.7+
- Flask

Install dependencies via pip:
```bash
pip install -r requirements.txt
```
