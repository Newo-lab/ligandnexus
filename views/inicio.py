"""
Vista: Inicio
=============
"""

import streamlit as st

import ligandnexus as lf
from . import common


def render():
    st.markdown("<div style='text-align:right;color:#9aa5b1;font-size:.8rem'>O. Molina</div>",
                unsafe_allow_html=True)
    st.title("🔗 LigandNexus")
    st.caption(f"Cribado virtual de derivados moleculares a partir de PubChem · v{lf.__version__}")

    st.write(
        "LigandNexus recupera los **derivados** de una molécula base desde **PubChem**, "
        "los depura y los prepara para estudios computacionales. Elija el flujo que "
        "necesite:")

    st.markdown(
        """<div class="lf-card">
        <h3>🩸 Metales + polifenoles</h3>
        <p>Coordinación de <b>Fe(II)/Fe(III)</b> y <b>Cu(I)/Cu(II)</b> por polifenoles:
        <b>taninos</b> (ácido gálico, ácido tánico) y <b>flavonoides</b> (quercetina,
        luteolina). Detecta los sitios catecol, galoílo, salicilato e
        <b>hidroxi-cetona</b>, estima la especiación vs pH y construye el complejo 3D
        <b>mono o polinuclear</b> como semilla para DFT.</p>
        </div>""", unsafe_allow_html=True)

    c1, c2 = st.columns(2)
    with c1:
        st.markdown(
            f"""<div class="lf-card">
            <h3>🔎 Cribado general</h3>
            <p>Solo necesita los <b>datos de las moléculas</b>: derivados de PubChem con
            sus propiedades, filtrado por rangos y exportación a CSV/Excel. Sin suponer
            química de coordinación.</p>
            </div>""", unsafe_allow_html=True)
        st.markdown(
            f"""<div class="lf-card">
            <h3>🎯 Átomos coordinantes</h3>
            <p>Identifica y <b>numera los átomos donadores</b> (N, O, S) de una molécula,
            coloreados por fuerza. Base para numerar posiciones.</p>
            </div>""", unsafe_allow_html=True)
    with c2:
        st.markdown(
            f"""<div class="lf-card">
            <h3>🧭 Análisis de coordinación</h3>
            <p>Pipeline completo para <b>ligandos</b>: depura, clasifica sustituyentes
            por grupo funcional y posición respecto al donador, y selecciona un conjunto
            representativo con balance donador/aceptor.</p>
            </div>""", unsafe_allow_html=True)
        st.markdown(
            f"""<div class="lf-card">
            <h3>🧬 Geometrías 3D</h3>
            <p>Pre-optimiza geometrías con campo de fuerzas, las muestra en un
            <b>visor 3D interactivo</b> y genera inputs de Gaussian / ORCA / Psi4.</p>
            </div>""", unsafe_allow_html=True)

    st.info("Use el menú de la izquierda para navegar entre modos. Los conjuntos que "
            "genere en un modo quedan disponibles en **Geometrías 3D**.")

    with st.expander("🚀 Inicio rápido"):
        st.markdown(
            "1. **Cribado general** → escriba `1474` (bipiridina) → *Verificar* → "
            "*Buscar derivados* → explore y exporte.\n"
            "2. **Átomos coordinantes** → escriba `glycine` o un SMILES → *Analizar*.\n"
            "3. **Análisis de coordinación** → recorra las pestañas ①→⑤.\n"
            "4. **Geometrías 3D** → una molécula o un conjunto → *Generar*.")
