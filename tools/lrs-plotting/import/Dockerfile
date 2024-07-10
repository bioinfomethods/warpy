FROM python:3.11.9-slim-bookworm

WORKDIR /opt/lrs-plotting

COPY prepare-alignments.py requirements.txt ./

RUN pip install --no-cache-dir -r requirements.txt
