"""
Flask API for Splicify
Main application file
"""

from flask import Flask, request, jsonify
from flask_cors import CORS

# Import the plannotate blueprint
from plannotate_endpoint import plannotate_bp

# Create Flask app
app = Flask(__name__)

# Enable CORS (allows frontend to call API)
CORS(app)

# Register blueprints
app.register_blueprint(plannotate_bp, url_prefix='/plannotate')

# Basic health check for main app
@app.route('/health', methods=['GET'])
def health():
    return jsonify({
        "ok": True,
        "service": "Splicify API",
        "version": "1.0"
    })

# Example route (you can add your Gibson routes here)
@app.route('/', methods=['GET'])
def index():
    return jsonify({
        "message": "Splicify API",
        "endpoints": {
            "main_health": "/health",
            "plannotate_health": "/plannotate/health",
            "plannotate_annotate": "/plannotate/annotate_genbank"
        }
    })

# Run the app
if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=True)
