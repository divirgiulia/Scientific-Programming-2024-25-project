from sequence_utils import read_sequence, is_valid_dna
from flask import Flask, request, jsonify
from flask_cors import CORS
import configparser
from analysis import run_analysis

config = configparser.ConfigParser()
config.read('config.ini')
HOST = config['server']['host'] # in this way i dont hard code the host and port variables
PORT = int(config['server']['port'])

app = Flask(__name__)
CORS(app)

@app.route("/analyze", methods=["POST"])
def analyze():
    try:
        data = request.get_json()
        task = data.get('task')
        sequence_input = data.get('sequence')
        window = data.get('window')
        step = data.get('step')
        threshold = data.get('threshold', 60)

        sequence = read_sequence(sequence_input)

        if not is_valid_dna(sequence):
            return jsonify({'error': 'Invalid DNA sequence.'}), 400
        
        results = run_analysis(task, sequence, window, step, threshold)

        return jsonify({"task": task, "results": results})
    
    except Exception as e:
        print(f"Error during analysis: {e}")
        return jsonify({'error': str(e)}), 500
    
if __name__ == '__main__':
        print(f"Starting server on http://{HOST}:{PORT}")
        app.run(host=HOST, port=PORT, debug=True)
