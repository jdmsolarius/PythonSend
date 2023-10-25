#!/bin/sh

# Set the environment variable
export PYTHONPATH=/app/ProteinRepo

# Run Flask
flask run --host=0.0.0.0 --port=5555