from fastapi import FastAPI
from pydantic import BaseModel
from rdkit import Chem
from fastapi.middleware.cors import CORSMiddleware

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # you can restrict to your Cloudflare domain later
    allow_methods=["*"],
    allow_headers=["*"],
)


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
from fastapi.responses import JSONResponse
from rdkit.Chem import Draw
from io import BytesIO
import base64

# ---------- Molecule render ----------
@app.post("/render")
def render_molecule(data: MoleculeRequest):
    mol = Chem.MolFromSmiles(data.smiles)
    if mol is None:
        return {"error": "Invalid SMILES"}

    # Generate image
    img = Draw.MolToImage(mol, size=(300,300))

    # Convert to bytes
    buf = BytesIO()
    img.save(buf, format="PNG")
    byte_data = buf.getvalue()

    # Encode as base64
    b64_img = base64.b64encode(byte_data).decode("utf-8")

    # Prepare atom info (optional for frontend)
    atoms = [{"index": a.GetIdx(), "symbol": a.GetSymbol()} for a in mol.GetAtoms()]

    return JSONResponse({
        "image": b64_img,
        "atoms": atoms,
        "bonds": []  # leave empty for now
    })

@app.post("/suggest")
def suggest(data: dict):
    # Dummy response
    return {"arrows": []}
