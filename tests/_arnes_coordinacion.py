# -*- coding: utf-8 -*-
"""
_arnes_coordinacion.py — arnes de AppTest para views/coordinacion
=================================================================
Monta la vista de coordinacion con las fuentes SIMULADAS (sin red, sin PubChem)
para poder conducirla en las pruebas. Dos lotes de moleculas distintos, para
comprobar que cambiar de molecula no arrastra resultados del lote anterior.

No es un test: lo ejecuta AppTest desde tests/verificar_motores.py.
"""
import os
import sys

_RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _RAIZ)

import pandas as pd
import streamlit as st

st.set_page_config(page_title="arnes", layout="wide")
from views import common, coordinacion

LOTES = {
    "1318": dict(
        info={"CID": 1318, "SMILES": "c1cnc2c(c1)ccc1cccnc21",
              "Title": "1,10-phenanthroline", "IUPACName": "1,10-phenanthroline",
              "MolecularFormula": "C12H8N2", "MolecularWeight": "180.21"},
        filas=[
            {"CID": 1318, "SMILES": "c1cnc2c(c1)ccc1cccnc21", "IUPACName": "phen",
             "MolecularWeight": 180.21, "XLogP": 2.0, "TPSA": 25.8},
            {"CID": 70848, "SMILES": "Cc1cnc2c(c1)ccc1cccnc21", "IUPACName": "metil-phen",
             "MolecularWeight": 194.24, "XLogP": 2.4, "TPSA": 25.8},
            {"CID": 111111, "SMILES": "O=[N+]([O-])c1cnc2c(c1)ccc1cccnc21",
             "IUPACName": "nitro-phen", "MolecularWeight": 225.20, "XLogP": 1.9,
             "TPSA": 71.6},
            # Pesada y polar: cae con el tope de peso y es la unica «no atraviesa»,
            # asi que sirve para probar filtros que dejan un subconjunto propio.
            {"CID": 333333, "SMILES": "OC(=O)c1cnc2c(c1)ccc1cc(C(=O)O)cnc21",
             "IUPACName": "phen-dicarboxi", "MolecularWeight": 268.22, "XLogP": 1.1,
             "TPSA": 100.4},
        ],
    ),
    "1474": dict(
        info={"CID": 1474, "SMILES": "c1ccc(-c2ccccn2)nc1",
              "Title": "2,2'-bipyridine", "IUPACName": "2,2'-bipyridine",
              "MolecularFormula": "C10H8N2", "MolecularWeight": "156.18"},
        filas=[
            {"CID": 1474, "SMILES": "c1ccc(-c2ccccn2)nc1", "IUPACName": "bpy",
             "MolecularWeight": 156.18, "XLogP": 1.8, "TPSA": 25.8},
            {"CID": 222222, "SMILES": "Cc1ccnc(-c2ccccn2)c1", "IUPACName": "metil-bpy",
             "MolecularWeight": 170.21, "XLogP": 2.2, "TPSA": 25.8},
        ],
    ),
}


def _resolver(texto, modo):
    t = str(texto).strip()
    if t in LOTES:
        return LOTES[t]["info"], None
    return None, f"No encontrado: {t}"


def _buscar(entry, smiles_core, max_records, log, subido=None):
    for k, v in LOTES.items():
        if v["info"]["SMILES"] == smiles_core:
            log(f"  simulado: {len(v['filas'])} filas del lote {k}")
            return pd.DataFrame(v["filas"])
    return pd.DataFrame()


common.resolver_entrada = _resolver
common.buscar_derivados_fuente = _buscar
common.selector_fuente = lambda key: ({"nombre": "SIMULADO"}, None)

coordinacion.render()
