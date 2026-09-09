"""
sources.py — Fuentes de datos intercambiables para el cribado
=============================================================
LigandNexus no depende de una sola base de datos. Aquí vive una abstracción de
«fuente»: cada una sabe buscar derivados (por subestructura) y devolver un
DataFrame NORMALIZADO con el mismo esquema, para que el resto de la app (filtros,
clasificación, geometrías) sea agnóstico a la procedencia.

Tres tipos de fuente:
  · EN VIVO por API con búsqueda de subestructura: PubChem, ChEMBL, COCONUT.
  · ARCHIVO LOCAL: el usuario sube un SDF/SMILES/CSV y la búsqueda de
    subestructura la hace RDKit localmente — esto habilita CUALQUIER base de
    datos que publique un volcado descargable (LOTUS, NPASS, ZINC, HMDB, …).
  · REGISTRO: catálogo de ~20 bases conocidas; las que no tienen API de
    subestructura se usan por la vía de descarga + archivo local.

Esquema normalizado de salida (columnas):
  ID, SMILES, MolecularFormula, MolecularWeight, IUPACName, Name,
  HBondDonorCount, HBondAcceptorCount, XLogP, TPSA, Source
Los campos que la fuente no provea se calculan con RDKit.
"""

from __future__ import annotations

import io
import csv
import json
import time
import urllib.request
import urllib.parse
import urllib.error

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Lipinski

from . import pubchem, chem

# Se reutiliza el contexto TLS de pubchem (verificacion activa salvo que se
# pida lo contrario con LIGANDNEXUS_TLS_INSEGURO=1). Un solo sitio donde
# cambiarlo evita que las dos mitades del programa queden desalineadas.
_SSL = pubchem._SSL_CTX

_COLS = ["ID", "SMILES", "MolecularFormula", "MolecularWeight", "IUPACName",
         "Name", "HBondDonorCount", "HBondAcceptorCount", "XLogP", "TPSA", "Source"]


# ── Normalización / enriquecimiento con RDKit ─────────────────────────────────

def _enriquecer(registros, source_id, log=print) -> pd.DataFrame:
    """
    Recibe una lista de dicts (cada uno con al menos 'SMILES') y devuelve un
    DataFrame con el esquema normalizado, calculando con RDKit lo que falte
    (fórmula, peso, donores/aceptores de H, logP, TPSA). Descarta SMILES inválidos.
    """
    filas = []
    n_malos = 0
    for r in registros:
        smi = str(r.get("SMILES", "") or "")
        mol = Chem.MolFromSmiles(smi) if smi else None
        if mol is None:
            n_malos += 1
            continue
        fila = {c: r.get(c, "") for c in _COLS}
        fila["SMILES"] = smi
        fila["Source"] = source_id
        # Rellenar con RDKit lo que falte.
        if not fila.get("MolecularFormula"):
            fila["MolecularFormula"] = rdMolDescriptors.CalcMolFormula(mol)
        if fila.get("MolecularWeight") in ("", None):
            fila["MolecularWeight"] = round(Descriptors.MolWt(mol), 2)
        if fila.get("HBondDonorCount") in ("", None):
            fila["HBondDonorCount"] = Lipinski.NumHDonors(mol)
        if fila.get("HBondAcceptorCount") in ("", None):
            fila["HBondAcceptorCount"] = Lipinski.NumHAcceptors(mol)
        if fila.get("XLogP") in ("", None):
            try:
                fila["XLogP"] = round(Descriptors.MolLogP(mol), 2)
            except Exception:
                fila["XLogP"] = ""
        if fila.get("TPSA") in ("", None):
            fila["TPSA"] = round(rdMolDescriptors.CalcTPSA(mol), 1)
        filas.append(fila)
    if n_malos:
        log(f"  ({n_malos} SMILES no válidos descartados)")
    df = pd.DataFrame(filas, columns=_COLS)
    # 'CID' como alias del ID genérico, para compatibilidad con el resto de la app.
    df["CID"] = df["ID"]
    return df


def _get(url, timeout=45, headers=None):
    h = {"User-Agent": "LigandNexus/0.2 (research)", "Accept": "application/json"}
    if headers:
        h.update(headers)
    req = urllib.request.Request(url, headers=h)
    with urllib.request.urlopen(req, timeout=timeout, context=_SSL) as r:
        return r.read()


def _post_json(url, payload, timeout=45):
    data = json.dumps(payload).encode()
    req = urllib.request.Request(url, data=data, method="POST", headers={
        "Content-Type": "application/json", "Accept": "application/json",
        "User-Agent": "LigandNexus/0.2 (research)"})
    with urllib.request.urlopen(req, timeout=timeout, context=_SSL) as r:
        return json.loads(r.read())


# ── Fuente base ───────────────────────────────────────────────────────────────

class Source:
    id = ""
    name = ""
    def buscar(self, smiles_core, max_records=200, log=print) -> pd.DataFrame:
        raise NotImplementedError


# ── PubChem (envuelve el cliente robusto ya existente) ────────────────────────

class PubChemSource(Source):
    id, name = "pubchem", "PubChem"

    def buscar(self, smiles_core, max_records=200, log=print) -> pd.DataFrame:
        cids = pubchem.buscar_derivados(smiles_core, max_records=max_records, log=log)
        if not cids:
            return _enriquecer([], self.id, log)
        df = pubchem.descargar_propiedades(cids, log=log)
        registros = []
        for _, row in df.iterrows():
            registros.append({
                "ID": row.get("CID", ""), "SMILES": pubchem.extraer_smiles(row.to_dict()),
                "MolecularFormula": row.get("MolecularFormula", ""),
                "MolecularWeight": row.get("MolecularWeight", ""),
                "IUPACName": row.get("IUPACName", ""), "Name": "",
                "HBondDonorCount": row.get("HBondDonorCount", ""),
                "HBondAcceptorCount": row.get("HBondAcceptorCount", ""),
                "XLogP": row.get("XLogP", ""), "TPSA": row.get("TPSA", ""),
            })
        return _enriquecer(registros, self.id, log)


# ── ChEMBL (API REST de subestructura del EBI) ────────────────────────────────

class ChEMBLSource(Source):
    id, name = "chembl", "ChEMBL"
    BASE = "https://www.ebi.ac.uk/chembl/api/data"

    def buscar(self, smiles_core, max_records=200, log=print) -> pd.DataFrame:
        limit = min(200, max_records)
        url = (f"{self.BASE}/substructure/{urllib.parse.quote(smiles_core)}.json"
               f"?limit={limit}")
        registros = []
        log(f"ChEMBL: subestructura {smiles_core}")
        while url and len(registros) < max_records:
            try:
                d = json.loads(_get(url))
            except urllib.error.HTTPError as e:
                log(f"  ChEMBL HTTP {e.code}"); break
            except Exception as e:
                log(f"  ChEMBL error: {e}"); break
            total = d.get("page_meta", {}).get("total_count", "?")
            for m in d.get("molecules", []):
                st = m.get("molecule_structures") or {}
                pr = m.get("molecule_properties") or {}
                registros.append({
                    "ID": m.get("molecule_chembl_id", ""),
                    "SMILES": st.get("canonical_smiles", ""),
                    "MolecularFormula": pr.get("full_molformula", ""),
                    "MolecularWeight": pr.get("full_mwt", ""),
                    "IUPACName": "", "Name": m.get("pref_name") or "",
                    "HBondDonorCount": pr.get("hbd", ""),
                    "HBondAcceptorCount": pr.get("hba", ""),
                    "XLogP": pr.get("alogp", ""), "TPSA": pr.get("psa", ""),
                })
                if len(registros) >= max_records:
                    break
            log(f"  ChEMBL: {len(registros):,}/{total}")
            nxt = d.get("page_meta", {}).get("next")
            url = ("https://www.ebi.ac.uk" + nxt) if nxt else None
            time.sleep(0.2)
        return _enriquecer(registros, self.id, log)


# ── COCONUT (productos naturales; POST /api/search) ───────────────────────────

class COCONUTSource(Source):
    id, name = "coconut", "COCONUT"
    URL = "https://coconut.naturalproducts.net/api/search"

    def buscar(self, smiles_core, max_records=200, log=print) -> pd.DataFrame:
        per = min(100, max_records)
        registros, page = [], 1
        log(f"COCONUT: subestructura {smiles_core}")
        while len(registros) < max_records:
            try:
                d = _post_json(self.URL, {"query": smiles_core, "type": "substructure",
                                          "limit": per, "page": page})
            except urllib.error.HTTPError as e:
                log(f"  COCONUT HTTP {e.code}"); break
            except Exception as e:
                log(f"  COCONUT error: {e}"); break
            data = d.get("data", {})
            lote = data.get("data", []) or []
            if not lote:
                break
            for r in lote:
                registros.append({
                    "ID": r.get("identifier", ""), "SMILES": r.get("canonical_smiles", ""),
                    "MolecularFormula": r.get("molecular_formula") or "",
                    "MolecularWeight": r.get("molecular_weight") or "",
                    "IUPACName": r.get("iupac_name") or "", "Name": r.get("name") or "",
                    "HBondDonorCount": "", "HBondAcceptorCount": "", "XLogP": "", "TPSA": "",
                })
                if len(registros) >= max_records:
                    break
            total = data.get("total", "?")
            log(f"  COCONUT: {len(registros):,}/{total}")
            if page >= data.get("last_page", page):
                break
            page += 1
            time.sleep(0.3)
        return _enriquecer(registros, self.id, log)


# ── Archivo local (SDF / SMILES / CSV) + subestructura con RDKit ──────────────

class LocalFileSource(Source):
    id, name = "local", "Archivo local"

    def _iter_mols(self, contenido: bytes, nombre: str):
        """Itera (mol, id) de un archivo SDF, SMILES (.smi/.txt) o CSV con SMILES."""
        low = nombre.lower()
        if low.endswith((".sdf", ".sd", ".mol")):
            supl = Chem.ForwardSDMolSupplier(io.BytesIO(contenido))
            for i, m in enumerate(supl):
                if m is not None:
                    idv = m.GetProp("_Name") if m.HasProp("_Name") else f"row{i}"
                    yield m, (idv or f"row{i}")
        elif low.endswith(".csv"):
            texto = contenido.decode("utf-8", errors="replace")
            rd = csv.DictReader(io.StringIO(texto))
            scol = next((c for c in (rd.fieldnames or []) if "smiles" in c.lower()), None)
            idcol = next((c for c in (rd.fieldnames or [])
                          if c.lower() in ("id", "cid", "identifier", "name")), None)
            for i, row in enumerate(rd):
                smi = row.get(scol, "") if scol else ""
                m = Chem.MolFromSmiles(str(smi)) if smi else None
                if m is not None:
                    yield m, (row.get(idcol) if idcol else f"row{i}")
        else:   # .smi / .txt: 'SMILES [id]' por línea
            texto = contenido.decode("utf-8", errors="replace")
            for i, linea in enumerate(texto.splitlines()):
                linea = linea.strip()
                if not linea or linea.lower().startswith("smiles"):
                    continue
                partes = linea.split()
                m = Chem.MolFromSmiles(partes[0])
                if m is not None:
                    yield m, (partes[1] if len(partes) > 1 else f"row{i}")

    def buscar_local(self, smiles_core, contenido, nombre, max_records=200,
                     log=print) -> pd.DataFrame:
        core = chem.query_nucleo(smiles_core) if smiles_core else None
        registros, n_leidos = [], 0
        log(f"Archivo local «{nombre}»: buscando subestructura…")
        for mol, idv in self._iter_mols(contenido, nombre):
            n_leidos += 1
            if core is not None and not mol.HasSubstructMatch(core):
                continue
            registros.append({"ID": idv, "SMILES": Chem.MolToSmiles(mol),
                              "MolecularFormula": "", "MolecularWeight": "",
                              "IUPACName": "", "Name": str(idv),
                              "HBondDonorCount": "", "HBondAcceptorCount": "",
                              "XLogP": "", "TPSA": ""})
            if len(registros) >= max_records:
                break
            if n_leidos % 2000 == 0:
                log(f"  leídas {n_leidos:,} moléculas, {len(registros)} coinciden…")
        log(f"  {len(registros)} coincidencias de {n_leidos:,} moléculas leídas.")
        return _enriquecer(registros, self.id, log)


# Instancias reutilizables.
PUBCHEM = PubChemSource()
CHEMBL = ChEMBLSource()
COCONUT = COCONUTSource()
LOCAL = LocalFileSource()


# ── Registro de bases de datos (~20) ──────────────────────────────────────────
# 'modo': "api" = búsqueda en vivo; "descarga" = bajar el volcado y usar la
# pestaña de archivo local; los de descarga traen enlace y formato.

REGISTRO = [
    # --- En vivo (API con subestructura) ---
    {"id": "pubchem", "nombre": "PubChem", "cat": "General", "modo": "api",
     "src": PUBCHEM, "desc": "≈119 M de compuestos. La fuente por defecto, la más amplia.",
     "url": "https://pubchem.ncbi.nlm.nih.gov/"},
    {"id": "chembl", "nombre": "ChEMBL", "cat": "Bioactividad / fármacos", "modo": "api",
     "src": CHEMBL, "desc": "≈2,4 M de moléculas con datos de bioactividad. API REST del EBI.",
     "url": "https://www.ebi.ac.uk/chembl/"},
    {"id": "coconut", "nombre": "COCONUT", "cat": "Productos naturales", "modo": "api",
     "src": COCONUT, "desc": "≈700 k productos naturales abiertos. Rico en motivos que unen metales.",
     "url": "https://coconut.naturalproducts.net/"},
    # --- Descarga + búsqueda local ---
    {"id": "lotus", "nombre": "LOTUS", "cat": "Productos naturales", "modo": "descarga",
     "desc": "≈750 k pares estructura–organismo. Volcado en Zenodo (SMILES/SDF).",
     "url": "https://lotus.naturalproducts.net/"},
    {"id": "npass", "nombre": "NPASS", "cat": "Productos naturales", "modo": "descarga",
     "desc": "Productos naturales con actividad biológica cuantitativa. Descarga TSV/SMILES.",
     "url": "https://bidd.group/NPASS/"},
    {"id": "npatlas", "nombre": "NP Atlas", "cat": "Productos naturales", "modo": "descarga",
     "desc": "Productos naturales microbianos curados. Descarga TSV/SDF.",
     "url": "https://www.npatlas.org/download"},
    {"id": "supernatural3", "nombre": "SuperNatural 3", "cat": "Productos naturales", "modo": "descarga",
     "desc": "≈560 k productos naturales. Descarga de estructuras.",
     "url": "https://bioinf-applied.charite.de/supernatural_3/"},
    {"id": "cmaup", "nombre": "CMAUP", "cat": "Plantas medicinales", "modo": "descarga",
     "desc": "Ingredientes activos de plantas medicinales. Descarga con SMILES.",
     "url": "https://bidd.group/CMAUP/"},
    {"id": "imppat", "nombre": "IMPPAT", "cat": "Plantas medicinales", "modo": "descarga",
     "desc": "Fitoquímicos de plantas medicinales indias. Descarga con SMILES.",
     "url": "https://cb.imsc.res.in/imppat/"},
    {"id": "chebi", "nombre": "ChEBI", "cat": "Interés biológico", "modo": "descarga",
     "desc": "Entidades químicas de interés biológico (ontología). Descarga SDF.",
     "url": "https://www.ebi.ac.uk/chebi/downloadsForward.do"},
    {"id": "drugbank", "nombre": "DrugBank", "cat": "Fármacos", "modo": "descarga",
     "desc": "Fármacos y dianas. Volcado con estructuras (requiere cuenta gratuita).",
     "url": "https://go.drugbank.com/releases/latest#structures"},
    {"id": "drugcentral", "nombre": "DrugCentral", "cat": "Fármacos", "modo": "descarga",
     "desc": "Fármacos aprobados y en desarrollo. Descarga SMILES.",
     "url": "https://drugcentral.org/download"},
    {"id": "bindingdb", "nombre": "BindingDB", "cat": "Bioactividad", "modo": "descarga",
     "desc": "Afinidades de unión proteína–ligando. Descarga SDF/TSV.",
     "url": "https://www.bindingdb.org/bind/chemsearch/marvin/Download.jsp"},
    {"id": "surechembl", "nombre": "SureChEMBL", "cat": "Patentes", "modo": "descarga",
     "desc": "Química extraída de patentes. Volcados descargables.",
     "url": "https://www.surechembl.org/"},
    {"id": "hmdb", "nombre": "HMDB", "cat": "Metabolómica", "modo": "descarga",
     "desc": "Metaboloma humano. Descarga SDF de estructuras.",
     "url": "https://hmdb.ca/downloads"},
    {"id": "foodb", "nombre": "FooDB", "cat": "Metabolómica", "modo": "descarga",
     "desc": "Compuestos de alimentos. Descarga con estructuras.",
     "url": "https://foodb.ca/downloads"},
    {"id": "swisslipids", "nombre": "SwissLipids", "cat": "Lípidos", "modo": "descarga",
     "desc": "Lípidos conocidos y predichos. Descarga con estructuras.",
     "url": "https://www.swisslipids.org/#/downloads"},
    {"id": "zinc", "nombre": "ZINC20", "cat": "Cribado / comprables", "modo": "descarga",
     "desc": "Compuestos comprables para cribado virtual. Subconjuntos en SMILES.",
     "url": "https://zinc20.docking.org/"},
    {"id": "cod", "nombre": "Crystallography Open DB", "cat": "Estructuras cristalinas", "modo": "descarga",
     "desc": "Estructuras cristalinas abiertas (útil para complejos metálicos). Descarga CIF.",
     "url": "https://www.crystallography.net/cod/"},
    # --- Tu propia biblioteca ---
    {"id": "local", "nombre": "Archivo propio", "cat": "Tu biblioteca", "modo": "local",
     "src": LOCAL, "desc": "Sube tu SDF/SMILES/CSV y busca la subestructura localmente.",
     "url": ""},
]


def por_id(sid):
    return next((e for e in REGISTRO if e["id"] == sid), None)
