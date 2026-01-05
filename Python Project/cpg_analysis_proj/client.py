import argparse
import requests
import json
import os

def load_input(input_str):
    sequence = ""
    if os.path.isfile(input_str):
        with open(input_str, 'r') as file:
            lines = file.readlines()
        sequence = ''.join([line.strip() for line in lines if not line.startswith(">")])
    else:
        sequence = input_str.strip()

    sequence = sequence.replace(" ", "").upper()
    sequence = ''.join([base for base in sequence if base in "ACGT"])
    return sequence

def main():
    parser = argparse.ArgumentParser(description='Client for CpG analysis server')

    parser.add_argument("--host", type=str, required=True, help='Server host')
    parser.add_argument("--port", type=int, required=True, help="Server port")
    parser.add_argument('--task', type=str, choices=['global', 'sliding', 'island'], required=True, help='Task to perform')
    parser.add_argument('--input', type=str, required=True, help='DNA sequence or path to FASTA file')

    parser.add_argument('--window', type=int, help='Window size')
    parser.add_argument('--step', type=int, help='Step size')
    parser.add_argument('--threshold', type=float, help='CpG % threshold')

    args = parser.parse_args()

    if args.task == "sliding":
        if args.window is None or args.step is None:
            print("Error: --window and --step must be provided for sliding window analysis.")
            print("Example: python client.py --host 127.0.0.1 --port 5000 --task sliding --input ATCG --window 100 --step 20")
            return
        
    if args.task == "island":
        if args.window is None or args.step is None or args.threshold is None:
            print("Error: Error: --window, --step, and --threshold must all be provided for island detection.")
            print("Example: python client.py --host 127.0.0.1 --port 5000 --task island --input ATCG --window 200 --step 50 --threshold 60")
            return

    try:
        sequence = load_input(args.input)
    except Exception as e:
        print(f"Error loading input: {e}")
        return
    
    payload = {"sequence": sequence,
               "task": args.task
    }

    if args.task in ['sliding', 'island']:
        payload['window'] = args.window
        payload['step'] = args.step
        if args.task == 'island':
            payload['threshold'] = args.threshold

    url = f"http://{args.host}:{args.port}/analyze"

    try:
        response = requests.post(url, json=payload)
        response.raise_for_status()
        result = response.json()
        print(json.dumps(result, indent=2))
    except requests.exceptions.RequestException as e:
        print(f"Request failed: {e}")
    except ValueError:
        print("Server returned non-JSON response.")

if __name__ == '__main__':
    main()