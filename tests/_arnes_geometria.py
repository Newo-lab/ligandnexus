# -*- coding: utf-8 -*-
"""
_arnes_geometria.py — arnes de AppTest para views/geometria
===========================================================
Solo monta la vista: los conjuntos se siembran desde la prueba con
`at.session_state[...]`, que es justo lo que hay que comprobar (que ofrezca los
conjuntos ya filtrados por ADMET y no los vacios).

No es un test: lo ejecuta AppTest desde tests/verificar_motores.py.
"""
import os
import sys

_RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _RAIZ)

import streamlit as st

st.set_page_config(page_title="arnes-geo", layout="wide")
from views import geometria

geometria.render()
