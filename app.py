"""
LigandNexus — aplicación web (Streamlit)
=============================================
Cribado virtual de derivados moleculares a partir de PubChem, con cinco modos:
cribado general, análisis de coordinación, identificador de átomos coordinantes,
geometrías 3D + inputs de química cuántica, y complejos de metales con polifenoles.

Ejecutar:  streamlit run app.py
"""

import os
import sys

import streamlit as st

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

st.set_page_config(page_title="LigandNexus", page_icon="🔗", layout="wide")

from views import common, inicio, general, coordinacion, donadores, hierro, geometria, acerca  # noqa: E402
import ligandnexus as lf  # noqa: E402

common.aplicar_tema()

with st.sidebar:
    st.markdown("### 🔗 LigandNexus")
    st.caption(f"v{lf.__version__} · cribado desde PubChem")
    st.divider()

pg = st.navigation({
    "Inicio": [
        # la portada abre en el indice de modos, no en un modo concreto: la
        # herramienta es multiproposito y arrancar dentro de uno solo da una
        # primera impresion equivocada a quien la recibe sin contexto.
        st.Page(inicio.render, title="Inicio / otros modos", icon="🏠",
                url_path="inicio", default=True),
        st.Page(hierro.render, title="Metales + polifenoles", icon="🩸",
                url_path="hierro-polifenoles"),
    ],
    "Cribado": [
        st.Page(general.render, title="Cribado general", icon="🔎",
                url_path="cribado-general"),
        st.Page(coordinacion.render, title="Análisis de coordinación", icon="🧭",
                url_path="coordinacion"),
    ],
    "Herramientas": [
        st.Page(donadores.render, title="Átomos coordinantes", icon="🎯",
                url_path="atomos-coordinantes"),
        st.Page(geometria.render, title="Geometrías 3D", icon="🧬",
                url_path="geometrias-3d"),
    ],
    "Ayuda": [
        st.Page(acerca.render, title="Acerca de / Referencias", icon="📚",
                url_path="acerca"),
    ],
})

with st.sidebar:
    st.divider()
    if st.button("🔄 Reiniciar sesión", width="stretch"):
        st.session_state.clear()
        st.rerun()
    st.caption("Cada base de datos tiene sus propios términos de uso "
               "(ver «Acerca»). Software libre · MIT.")

pg.run()
