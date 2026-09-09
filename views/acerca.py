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
                      "Descripción": e["desc"], "Términos": e.get("url", "")})
    st.dataframe(pd.DataFrame(filas), hide_index=True, width="stretch",
                 column_config={"": st.column_config.TextColumn(width="small"),
                                "Términos": st.column_config.LinkColumn(
                                    display_text="consultar", width="small")})
    st.warning("**Las licencias no son iguales.** PubChem es de dominio público, "
               "pero ChEMBL se distribuye bajo CC BY-SA 3.0 (atribución y "
               "compartir-igual) y otras restringen el uso comercial o exigen "
               "cuenta. Revise los términos de la fuente que use antes de publicar "
               "o redistribuir lo que exporte.", icon="⚖️")
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
        "- **Cada base de datos tiene su propia licencia y usted es responsable de "
        "respetarla.** PubChem es de dominio público, pero no todas lo son: ChEMBL "
        "se distribuye bajo **CC BY-SA 3.0** (exige atribución y compartir-igual) y "
        "varias restringen el uso comercial. Consulte los términos en el enlace de "
        "cada fuente antes de publicar o redistribuir lo que exporte.\n"
        "- La app respeta por sí sola el límite de ~5 peticiones/segundo de PubChem.")
