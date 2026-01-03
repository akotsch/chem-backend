FROM continuumio/miniconda3

WORKDIR /app

# Install RDKit
RUN conda install -c conda-forge rdkit python=3.11 -y

# Install Python deps
RUN pip install fastapi uvicorn

COPY . .

EXPOSE 10000

CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "10000"]
