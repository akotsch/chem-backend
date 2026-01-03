from fastapi import FastAPI
from pydantic import BaseModel
from rdkit import Chem
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from rdkit.Chem import Draw
from io import BytesIO
import base64

# ---------- App ----------
app = FastAPI()

# ---------- Enable CORS ----------
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

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
    if not mol:
        return {"error": "Invalid SMILES"}

    atoms = [{"index": a.GetIdx(), "symbol": a.GetSymbol(),
              "degree": a.GetDegree(), "charge": a.GetFormalCharge()} for a in mol.GetAtoms()]

    return {"num_atoms": mol.GetNumAtoms(), "atoms": atoms}

# ---------- Molecule render ----------
@app.post("/render")
def render_molecule(data: MoleculeRequest):
    mol = Chem.MolFromSmiles(data.smiles)
    if not mol:
        return {"error": "Invalid SMILES"}

    img = Draw.MolToImage(mol, size=(300,300))
    buf = BytesIO()
    img.save(buf, format="PNG")
    b64_img = base64.b64encode(buf.getvalue()).decode("utf-8")

    atoms = [{"index": a.GetIdx(), "symbol": a.GetSymbol()} for a in mol.GetAtoms()]

    return JSONResponse({"image": b64_img, "atoms": atoms, "bonds": []})

# ---------- Suggest mechanism ----------
@app.post("/suggest")
def suggest(data: dict):
    # Dummy response for now
    return {"arrows": []}
