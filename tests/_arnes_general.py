# -*- coding: utf-8 -*-
"""
_arnes_general.py — arnes de AppTest para views/general
=======================================================
Igual que _arnes_coordinacion, pero para el cribado general: hace falta para
comprobar que mover un filtro por propiedad NO tira el perfil ADMET ya
calculado.

No es un test: lo ejecuta AppTest desde tests/verificar_motores.py.
"""
import os
import sys

_RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _RAIZ)

import pandas as pd
import streamlit as st

st.set_page_config(page_title="arnes-gen", layout="wide")
from views import common, general

INFO = {"CID": 1318, "SMILES": "c1cnc2c(c1)ccc1cccnc21", "Title": "1,10-phenanthroline",
        "IUPACName": "1,10-phenanthroline", "MolecularFormula": "C12H8N2",
        "MolecularWeight": "180.21"}

FILAS = [
    {"CID": 1318, "SMILES": "c1cnc2c(c1)ccc1cccnc21", "IUPACName": "phen",
     "MolecularWeight": 180.21, "XLogP": 2.0, "TPSA": 25.8},
    {"CID": 70848, "SMILES": "Cc1cnc2c(c1)ccc1cccnc21", "IUPACName": "metil-phen",
     "MolecularWeight": 194.24, "XLogP": 2.4, "TPSA": 25.8},
    {"CID": 111111, "SMILES": "O=[N+]([O-])c1cnc2c(c1)ccc1cccnc21",
     "IUPACName": "nitro-phen", "MolecularWeight": 225.20, "XLogP": 1.9, "TPSA": 71.6},
    {"CID": 333333, "SMILES": "OC(=O)c1cnc2c(c1)ccc1cc(C(=O)O)cnc21",
     "IUPACName": "phen-dicarboxi", "MolecularWeight": 268.22, "XLogP": 1.1,
     "TPSA": 100.4},
]

common.resolver_entrada = lambda t, m: ((INFO, None) if str(t).strip() == "1318"
                                        else (None, "no encontrado"))
common.buscar_derivados_fuente = lambda e, s, n, log, sub=None: pd.DataFrame(FILAS)
common.selector_fuente = lambda key: ({"nombre": "SIMULADO"}, None)

general.render()
