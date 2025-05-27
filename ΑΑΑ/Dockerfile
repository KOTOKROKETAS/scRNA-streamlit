# Dockerfile για Streamlit εφαρμογή ανάλυσης scRNA-seq

FROM python:3.10-slim

ENV PYTHONUNBUFFERED=1

WORKDIR /app

COPY requirements.txt ./

RUN apt-get update &&     apt-get install -y build-essential python3-dev libglib2.0-0 libsm6 libxrender1 libxext6 &&     pip install --no-cache-dir --upgrade pip &&     pip install --no-cache-dir -r requirements.txt

COPY . .

CMD ["streamlit", "run", "streamlit_pipeline_app.py", "--server.port=8501", "--server.enableCORS=false"]
