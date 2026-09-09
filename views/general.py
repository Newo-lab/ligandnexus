"""
Vista: Cribado general
======================
Para quien solo necesita, a partir de una molécula base, la lista de sus
derivados con sus propiedades — SIN suponer química de coordinación ni cobre.
Recupera de la base de datos elegida, explora/filtra por propiedades y exporta.
"""

import io
import pandas as pd
import streamlit as st

from ligandnexus import excel
from . import common


def render():
    # El conjunto perfilado que se publica para Geometrías 3D es CONSECUENCIA de
    # lo que quede en pantalla: se borra al entrar y solo se repone si al final
    # hay un perfil calculado. Así no puede sobrevivir a una búsqueda nueva ni a
    # una salida temprana de esta función.
    st.session_state.pop("gen_perfil", None)

    st.header("🔎 Cribado general de derivados")
    st.write("A partir de **una molécula base**, recupera sus derivados por "
             "subestructura en la **base de datos que elija** (PubChem, ChEMBL, "
             "COCONUT… o tu propio archivo). Pensado para exploración general: no "
             "asume química de coordinación. Para el análisis de cobre use la pestaña "
             "**Análisis de coordinación**.")

    # 1) Molécula base
    st.subheader("1 · Molécula base")
    modo = st.radio("Buscar por", ["CID de PubChem", "Nombre"], horizontal=True,
                    key="gen_modo")
    modo_key = "CID" if modo.startswith("CID") else "Nombre"
    c1, c2 = st.columns([3, 1])
    with c1:
        ph = "Ej.: 1474" if modo_key == "CID" else "Ej.: 2,2'-bipyridine"
        entrada = st.text_input("Molécula base", key="gen_entrada", placeholder=ph,
                                label_visibility="collapsed")
    with c2:
        buscar = st.button("🔎 Verificar", type="primary", width="stretch")

    with st.expander("¿No conoce el CID? Ejemplos"):
        st.dataframe(pd.DataFrame([{"Molécula": k, "CID": v}
                                   for k, v in common.CABEZAS.items()]),
                     hide_index=True, width="stretch")

    if buscar:
        info, err = common.resolver_entrada(entrada, modo_key)
        if err:
            st.error(err)
            st.session_state.pop("gen_info", None)
        else:
            st.session_state["gen_info"] = info
            st.session_state.pop("gen_df", None)

    info = st.session_state.get("gen_info")
    if not info:
        return

    with st.container(border=True):
        v1, v2 = st.columns([1, 2])
        with v1:
            im = common.img_mol(info["SMILES"], (260, 200))
            if im:
                st.image(im, width="stretch")
        with v2:
            st.markdown(f"### {info.get('Title') or info.get('IUPACName') or '—'}")
            st.write(f"**Fórmula:** {info.get('MolecularFormula','—')}  ·  "
                     f"**MW:** {info.get('MolecularWeight','—')} g/mol  ·  "
                     f"**CID:** {info['CID']}")
            st.code(info["SMILES"], language=None)

    # 2) Buscar derivados (en la base de datos elegida)
    st.subheader("2 · Recuperar derivados")
    entry, subido = common.selector_fuente("gen")
    max_rec = st.number_input("Máx. de derivados", 10, 100000, 200, 50, key="gen_maxrec",
                              help="200 para explorar rápido; súbelo para una corrida grande.")

    if st.button("📥 Buscar derivados", type="primary"):
        try:
            with st.status(f"Buscando en {entry['nombre']}…", expanded=True) as s:
                df = common.buscar_derivados_fuente(entry, info["SMILES"], int(max_rec),
                                                    lambda m: s.write(m), subido)
                if df is None or df.empty:
                    s.update(label="Sin resultados", state="error")
                    st.warning("No se encontraron derivados con esta subestructura.")
                    st.session_state.pop("gen_df", None)
                else:
                    st.session_state["gen_df"] = df
                    s.update(label=f"Completado · {len(df)} derivados de {entry['nombre']}",
                             state="complete")
        except Exception as e:
            st.error(f"Error al buscar en {entry['nombre']}: {e}")

    df = st.session_state.get("gen_df")
    if df is None or df.empty:
        return

    # 3) Explorar + filtrar por propiedades numéricas
    st.subheader("3 · Explorar y filtrar")
    st.metric("Derivados recuperados", len(df))

    dfx = df.copy()
    num_cols = [c for c in ("MolecularWeight", "XLogP", "TPSA",
                            "HBondDonorCount", "HBondAcceptorCount") if c in dfx.columns]
    with st.expander("⚙️ Filtros por propiedad (opcional)"):
        for c in num_cols:
            serie = pd.to_numeric(dfx[c], errors="coerce")
            if serie.notna().sum() == 0:
                continue
            lo, hi = float(serie.min()), float(serie.max())
            if lo == hi:
                continue
            r = st.slider(c, lo, hi, (lo, hi), key=f"gen_flt_{c}")
            dfx = dfx[(pd.to_numeric(dfx[c], errors="coerce").between(*r))]
        dfx = dfx.reset_index(drop=True)

    st.caption(f"Mostrando {len(dfx)} de {len(df)} derivados.")
    st.dataframe(dfx, width="stretch", height=340)

    # 4) Perfil ADMET — opcional, sobre lo que haya filtrado en pantalla. Lo que
    # devuelva el panel es lo que se exporta después, así que si se filtra por
    # Lipinski o por barrera, la descarga ya sale filtrada.
    st.subheader("4 · Farmacología")
    _scol = next((c for c in dfx.columns if "smiles" in c.lower()), None)
    dfx = common.panel_admet(dfx, _scol, key="gen_admet")
    if "Lipinski" in getattr(dfx, "columns", []):
        # Disponible en Geometrías 3D: lo que se ve aquí, ya filtrado.
        st.session_state["gen_perfil"] = dfx

    # 5) Exportar — sale SIEMPRE lo que hay en pantalla: si el panel ADMET filtró,
    # la descarga ya va filtrada.
    st.subheader("5 · Exportar")
    scol = next((c for c in dfx.columns if "smiles" in c.lower()), None)
    st.caption(f"{len(dfx):,} molécula(s) · {len(dfx.columns)} columnas.")
    if len(dfx) == 0:
        st.warning("No queda ninguna molécula tras los filtros: la descarga saldrá vacía.")

    e1, e2 = st.columns(2)
    with e1:
        # `data` como invocable: el CSV se arma al pulsar y no en cada rerun.
        st.download_button("⬇️ Descargar CSV", data=lambda: dfx.to_csv(index=False),
                           file_name="derivados.csv", mime="text/csv",
                           width="stretch", on_click="ignore")
    with e2:
        con_img = st.checkbox("Excel con estructuras 2D (más lento)", value=False)
        if con_img and not scol:
            st.warning("Este conjunto no trae columna de SMILES: el Excel saldrá "
                       "sin las estructuras 2D.")
        # El .xlsx se guarda en sesión y su botón se pinta FUERA del `if`: antes
        # vivía dentro y desaparecía en el primer rerun, incluido el que
        # provocaba su propio clic. La firma lo ata al conjunto con el que se
        # generó, para que nunca se descargue un Excel de otra cosa.
        firma_x = (common.firma_conjunto(dfx, scol), bool(con_img))
        if st.button("📊 Generar Excel", width="stretch"):
            buf = io.BytesIO()
            with st.status("Generando Excel…") as s:
                excel.exportar_excel(dfx, scol, buf, con_imagenes=con_img,
                                     log=lambda m: s.write(m))
                s.update(label="Listo", state="complete")
            st.session_state["gen_xlsx"] = (firma_x, buf.getvalue())
        guardado = st.session_state.get("gen_xlsx")
        if guardado and guardado[0] == firma_x:
            st.download_button("⬇️ Descargar Excel", data=guardado[1],
                               file_name="derivados.xlsx",
                               mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                               width="stretch", on_click="ignore")
        elif guardado:
            st.caption("El Excel guardado es de otro conjunto u otras opciones: "
                       "vuelva a generarlo.")
