#!/bin/bash

if [ -d ".venv" ]; then
  source .venv/bin/activate
fi

pip install -r requirements.txt

streamlit run streamlit_pipeline_app.py --server.port=8501 --server.enableCORS=false
