#!/usr/bin/env python3
'''
Test Gateway cloning with feature lookup
'''

import sys
import re
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

# Test feature extraction from message
message = "Please design primers for insertion of GFP into this vector via Gateway cloning"

# Common patterns for feature references
feature_patterns = [
    r"insertion of (\w+)",
    r"clone (\w+) into",
    r"insert (\w+)",
    r"design primers for (\w+)",
]

potential_features = []
for pattern in feature_patterns:
    matches = re.findall(pattern, message, re.IGNORECASE)
    potential_features.extend(matches)

common_words = {"this", "that", "the", "into", "via", "using", "with", "vector", "plasmid", "primers"}
potential_features = [f for f in potential_features if f.lower() not in common_words]

print(f"Message: {message}")
print(f"\nExtracted features: {potential_features}")

# Test feature lookup
from splicify_api.chat import lookup_genetic_features

if potential_features:
    feature_results = lookup_genetic_features(potential_features)
    print(f"\nFeature lookup results:")
    for feature in feature_results:
        print(f"  - {feature['name']}: {len(feature['sequence'])} bp from {feature.get('source', 'unknown')}")
else:
    print("\nNo features found")
