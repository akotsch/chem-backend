from fastapi import FastAPI
from pydantic import BaseModel
from rdkit import Chem

app = FastAPI()

# ---------- Request format ----------
class MoleculeRequest(BaseModel):
    smiles: str

# ---------- Root endpoint ----------
@app.get("/")
def root():
    return {"status": "Chem backend running"}

# ---------- Molecule analysis ----------
@app.post("/analyze")
def analyze_molecule(data: MoleculeRequest):
    mol = Chem.MolFromSmiles(data.smiles)
    if mol is None:
        return {"error": "Invalid SMILES"}

    atoms = []
    for atom in mol.GetAtoms():
        atoms.append({
            "index": atom.GetIdx(),
            "symbol": atom.GetSymbol(),
            "degree": atom.GetDegree(),
            "charge": atom.GetFormalCharge()
        })

    return {
        "num_atoms": mol.GetNumAtoms(),
        "atoms": atoms
    }
