"""
Vista: Acerca de / Referencias
==============================
"""

import pandas as pd
import streamlit as st
from ligandnexus import sources as _sources
from . import common


def _tabla_fuentes():
    st.subheader("Fuentes de datos disponibles")
    st.caption("🟢 búsqueda en vivo por API · ⬇️ descarga el volcado y búscalo localmente · "
               "📁 tu propio archivo. En todas, la búsqueda es por subestructura.")
    filas = []
    for e in _sources.REGISTRO:
        icono = {"api": "🟢", "local": "📁", "descarga": "⬇️"}.get(e["modo"], "•")
        filas.append({"": icono, "Base de datos": e["nombre"], "Categoría": e["cat"],
                      "Descripción": e["desc"]})
    st.dataframe(pd.DataFrame(filas), hide_index=True, width="stretch",
                 column_config={"": st.column_config.TextColumn(width="small")})
    st.caption("¿Falta alguna? Si publica un volcado en SDF/SMILES/CSV, ya es usable con "
               "la opción **Archivo propio**.")


def render():
    st.header("📚 Acerca de LigandNexus y cómo citarlo")
    st.write("LigandNexus es una herramienta de **software libre** (licencia MIT) para "
             "el cribado virtual de derivados moleculares. Si la utiliza en un trabajo "
             "publicable, cite el software que la compone, la(s) base(s) de datos que "
             "consulte y el programa de química cuántica que ejecute.")

    _tabla_fuentes()
    st.divider()

    estilo = st.radio("Estilo de cita", ["ACS", "APA"], horizontal=True, key="ac_estilo")

    st.subheader("Software de la aplicación")
    for r in common.REFERENCIAS:
        cita = r["acs"] if estilo == "ACS" else r["apa"]
        st.markdown(f"**{r['tool']}** · [documentación]({r['doc']})  \n"
                    f"<span style='font-size:.85rem'>{cita}</span>", unsafe_allow_html=True)

    st.subheader("Programas de química cuántica (cite el que ejecute)")
    for r in common.REFERENCIAS_QM:
        cita = r["acs"] if estilo == "ACS" else r["apa"]
        st.markdown(f"**{r['tool']}** · [sitio]({r['doc']})  \n"
                    f"<span style='font-size:.85rem'>{cita}</span>", unsafe_allow_html=True)

    st.subheader("Reglas y descriptores del panel ADMET")
    st.caption("Son criterios de la literatura, no software: si usa el tamiz "
               "farmacológico en un trabajo publicable, cite la regla concreta que "
               "aplique.")
    for r in common.REFERENCIAS_ADMET:
        cita = r["acs"] if estilo == "ACS" else r["apa"]
        st.markdown(f"**{r['tool']}** · [DOI]({r['doc']})  \n"
                    f"<span style='font-size:.85rem'>{cita}</span>", unsafe_allow_html=True)

    st.divider()
    st.subheader("Cómo citar LigandNexus")
    st.code("O. Molina. LigandNexus: cribado virtual de derivados moleculares a partir "
            "de PubChem, v0.2, 2026.", language=None)

    st.subheader("Notas técnicas")
    st.markdown(
        "- Las **geometrías 3D** son una pre-optimización con campo de fuerzas "
        "(MMFF94/UFF), **no** una optimización DFT.\n"
        "- Solo `.gjf` (Gaussian), `.inp` (ORCA) y `.in` (Psi4) se ejecutan en su "
        "programa; `.xyz/.sdf/.mol` son para visualización.\n"
        "- La identificación de átomos coordinantes es una **heurística** química; "
        "para casos límite conviene el criterio del especialista.\n"
        "- Los datos provienen de **PubChem** (dominio público); respete sus términos "
        "de uso y el límite de ~5 peticiones/segundo (la app lo controla sola).")
