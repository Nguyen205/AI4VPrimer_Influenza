"""
Web API for editing primer_pipeline.py CONFIG and PARAMETERS.
Run: python app.py
Then open http://localhost:5000
"""
from flask import Flask, jsonify, request, send_file
import re, os

app = Flask(__name__)
PIPELINE = os.path.join(os.path.dirname(__file__), 'primer_pipeline.py')

# Define which variables belong to CONFIG vs PARAMETERS
CONFIG_VARS = ['GENE_NAME', 'ALN_FILE', 'CROSS_PANELS', 'SPECIFICITY_FILE', 'OUT_FILE', 'PRIMERS_PER_REGION']
PARAM_VARS = [
    'MAX_PRIMER_LEN',
    'MIN_TM', 'MAX_LOW_IDENTITY_BRIDGE',
    'MIN_AMP_SEQ', 'MAX_AMP_SEQ', 'MIN_AMP_QPCR', 'MAX_AMP_QPCR',
    'MIN_PRIMER_SENS', 'SENS_TARGET', 'OVERLAP_LIMITS',
    'COMBO_DEGEN_PREFER'
]


def read_pipeline():
    with open(PIPELINE, 'r') as f:
        return f.read()


def parse_value(raw):
    """Parse a raw string value from the source code."""
    raw = raw.strip()
    # Strip inline comments (but not inside strings or complex structures)
    if not raw.startswith(("'", '"', '{', '[')):
        raw = raw.split('#')[0].strip()
    if raw in ('None', ''):
        return None
    if raw.startswith("'") or raw.startswith('"'):
        return raw.strip("'\"")
    if raw.startswith('{'):
        return eval(raw)
    if raw.startswith('['):
        return eval(raw)
    try:
        if '.' in raw:
            return float(raw)
        return int(raw)
    except ValueError:
        return raw


def extract_params(source, var_names):
    """Extract variable assignments from source code."""
    params = {}
    for var in var_names:
        # Match single-line assignments and multi-line dicts/lists
        pattern = rf'^{var}\s*=\s*(.+?)(?:\n(?=[A-Z_]|\n|#\s*[═]))'
        m = re.search(pattern, source, re.MULTILINE | re.DOTALL)
        if not m:
            # Try single line
            pattern2 = rf'^{var}\s*=\s*(.+)$'
            m = re.search(pattern2, source, re.MULTILINE)
        if m:
            raw = m.group(1).rstrip()
            # For multi-line values, keep collecting until braces/brackets close
            if raw.count('{') > raw.count('}') or raw.count('[') > raw.count(']'):
                start = m.start(1)
                depth_b, depth_c = 0, 0
                end = start
                for i, ch in enumerate(source[start:], start):
                    if ch == '{': depth_c += 1
                    elif ch == '}': depth_c -= 1
                    elif ch == '[': depth_b += 1
                    elif ch == ']': depth_b -= 1
                    if depth_b <= 0 and depth_c <= 0 and i > start:
                        end = i + 1
                        break
                raw = source[start:end].rstrip()
            params[var] = parse_value(raw)
    return params


def update_pipeline(var, value):
    """Update a single variable assignment in the pipeline source."""
    source = read_pipeline()
    # Format value for Python source
    if isinstance(value, str):
        formatted = f"'{value}'"
    elif value is None:
        formatted = 'None'
    elif isinstance(value, dict):
        # Format dict nicely
        if not value:
            formatted = '{}'
        else:
            items = ',\n    '.join(f"'{k}': '{v}'" for k, v in value.items())
            formatted = '{\n    ' + items + ',\n}'
    elif isinstance(value, list):
        formatted = repr(value)
    else:
        formatted = repr(value)

    # Replace the assignment line(s)
    # Handle multi-line values (dicts, lists spanning lines)
    pattern = rf'^({var}\s*=\s*)(.+?)(\n(?=[A-Z_]|\n|#\s*[═]))'
    m = re.search(pattern, source, re.MULTILINE | re.DOTALL)
    if not m:
        pattern = rf'^({var}\s*=\s*)(.+)$'
        m = re.search(pattern, source, re.MULTILINE)

    if m:
        # For multi-line original values, find the full extent
        raw = m.group(2).rstrip()
        if raw.count('{') > raw.count('}') or raw.count('[') > raw.count(']'):
            start = m.start(2)
            depth_b, depth_c = 0, 0
            end = start
            for i, ch in enumerate(source[start:], start):
                if ch == '{': depth_c += 1
                elif ch == '}': depth_c -= 1
                elif ch == '[': depth_b += 1
                elif ch == ']': depth_b -= 1
                if depth_b <= 0 and depth_c <= 0 and i > start:
                    end = i + 1
                    break
            source = source[:m.start(2)] + formatted + source[end:]
        else:
            source = source[:m.start(2)] + formatted + source[m.end(2):]

        with open(PIPELINE, 'w') as f:
            f.write(source)
        return True
    return False


@app.route('/')
def index():
    return send_file('index.html')


@app.route('/api/params', methods=['GET'])
def get_params():
    source = read_pipeline()
    config = extract_params(source, CONFIG_VARS)
    params = extract_params(source, PARAM_VARS)
    # Show config fields as blank so users fill their own
    for key in config:
        config[key] = None
    return jsonify({
        'config': config,
        'parameters': params
    })


@app.route('/api/params', methods=['POST'])
def set_params():
    data = request.json
    errors = []
    for var, value in data.items():
        if var not in CONFIG_VARS + PARAM_VARS:
            errors.append(f"Unknown variable: {var}")
            continue
        if not update_pipeline(var, value):
            errors.append(f"Failed to update: {var}")
    if errors:
        return jsonify({'status': 'partial', 'errors': errors}), 207
    return jsonify({'status': 'ok'})


@app.route('/api/run', methods=['POST'])
def run_pipeline():
    import subprocess
    result = subprocess.run(['python', PIPELINE], capture_output=True, text=True, cwd=os.path.dirname(PIPELINE))
    return jsonify({'stdout': result.stdout, 'stderr': result.stderr, 'returncode': result.returncode})


if __name__ == '__main__':
    import webbrowser, threading
    threading.Timer(1, lambda: webbrowser.open('http://127.0.0.1:5000')).start()
    app.run(host='127.0.0.1', port=5000)
