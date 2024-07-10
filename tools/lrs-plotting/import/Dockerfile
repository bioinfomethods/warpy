FROM python:3.9.19-slim-bookworm

RUN pip install --no-cache-dir --upgrade pip

WORKDIR /opt/lrs-plotting

COPY prepare-alignments.py requirements.txt ./

RUN pip install --no-cache-dir -r requirements.txt
