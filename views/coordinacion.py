"""
Vista: Análisis de coordinación (pipeline completo)
===================================================
El flujo orientado a química de coordinación (p. ej. quelantes de cobre):
identificar la molécula → depurar → clasificar sustituyentes por grupo funcional
y posición respecto al donador → seleccionar un conjunto representativo → exportar.
"""

import io
import pandas as pd
import streamlit as st

from ligandnexus import screening, excel, chem
from . import common


def _hecho(k):
    return st.session_state.get(k) is not None


def render():
    st.header("🧭 Análisis de coordinación")
    st.write("Pipeline completo para preparar familias de **ligandos**: depura la "
             "biblioteca, clasifica los sustituyentes por **grupo funcional** y "
             "**posición** respecto al átomo coordinante, y selecciona un conjunto "
             "representativo con balance donador/aceptor.")

    tabs = st.tabs(["① Molécula", "② Depurar", "③ Clasificar",
                    "④ Seleccionar", "💊 ADMET", "⑤ Exportar"])

    # ── ① Molécula ────────────────────────────────────────────────────────────
    with tabs[0]:
        modo = st.radio("Buscar por", ["CID de PubChem", "Nombre"], horizontal=True,
                        key="coo_modo")
        modo_key = "CID" if modo.startswith("CID") else "Nombre"
        c1, c2 = st.columns([3, 1])
        with c1:
            ph = "Ej.: 1474" if modo_key == "CID" else "Ej.: 1,10-phenanthroline"
            entrada = st.text_input("Molécula", key="coo_entrada", placeholder=ph,
                                    label_visibility="collapsed")
        with c2:
            buscar = st.button("🔎 Verificar", type="primary", width="stretch")
        with st.expander("Ejemplos (cabezas de serie)"):
            st.dataframe(pd.DataFrame([{"Molécula": k, "CID": v}
                                       for k, v in common.CABEZAS.items()]),
                         hide_index=True, width="stretch")

        if buscar:
            info, err = common.resolver_entrada(entrada, modo_key)
            if err:
                st.error(err)
            else:
                st.session_state["coo_info"] = info
                for k in ("coo_df", "coo_filt", "coo_clas", "coo_sel"):
                    st.session_state.pop(k, None)

        info = st.session_state.get("coo_info")
        if info:
            with st.container(border=True):
                v1, v2 = st.columns([1, 2])
                with v1:
                    img, dons = chem.imagen_donadores(info["SMILES"], (300, 230),
                                                      dark=common.es_oscuro())
                    if img:
                        st.image(img, width="stretch")
                    st.caption("Donadores resaltados (verde=fuerte, ámbar=moderado).")
                with v2:
                    st.markdown(f"### {info.get('Title') or info.get('IUPACName') or '—'}")
                    st.write(f"**Fórmula:** {info.get('MolecularFormula','—')}  ·  "
                             f"**MW:** {info.get('MolecularWeight','—')}  ·  "
                             f"**CID:** {info['CID']}")
                    st.code(info["SMILES"], language=None)
                    if dons:
                        st.write("**Átomos coordinantes:** " +
                                 ", ".join(f"{d['elemento']}({d['nivel']})" for d in dons))

            st.subheader("Recuperar derivados")
            entry, subido = common.selector_fuente("coo")
            max_rec = st.number_input("Máx. de derivados", 10, 100000, 200, 50,
                                      key="coo_maxrec")
            if st.button("📥 Buscar derivados", type="primary", key="coo_buscar_der"):
                try:
                    with st.status(f"Buscando en {entry['nombre']}…", expanded=True) as s:
                        df = common.buscar_derivados_fuente(
                            entry, info["SMILES"], int(max_rec), lambda m: s.write(m), subido)
                        if df is None or df.empty:
                            s.update(label="Sin resultados", state="error")
                        else:
                            st.session_state.update(coo_df=df, coo_scol="SMILES",
                                                    coo_nucleo=info["SMILES"],
                                                    coo_mwpadre=float(info["MolecularWeight"] or 0))
                            for k in ("coo_filt", "coo_clas", "coo_sel"):
                                st.session_state.pop(k, None)
                            s.update(label=f"Completado · {len(df)} derivados de {entry['nombre']}",
                                     state="complete")
                except Exception as e:
                    st.error(f"Error: {e}")

            if _hecho("coo_df"):
                st.metric("Derivados", len(st.session_state["coo_df"]))
                st.dataframe(st.session_state["coo_df"].head(12), width="stretch")
                st.success("Continúe en **② Depurar**.")

    # ── ② Depurar ─────────────────────────────────────────────────────────────
    with tabs[1]:
        if not _hecho("coo_df"):
            st.info("Complete el paso ①.")
        else:
            f_met = st.checkbox("Descartar metales/semimetales", True)
            f_sal = st.checkbox("Descartar sales y mezclas", True)
            f_car = st.checkbox("Descartar especies con carga", True)
            f_nuc = st.checkbox("Descartar los que no contienen el núcleo", True)
            f_mw = st.checkbox("Descartar los demasiado pesados", True)
            mw_rel = st.number_input("Tope = peso del núcleo + …", 0.0, 2000.0, 350.0, 50.0)
            if st.button("🧹 Aplicar filtros", type="primary"):
                with st.status("Depurando…", expanded=True) as s:
                    dff, ct = screening.aplicar_filtros(
                        st.session_state["coo_df"], st.session_state["coo_scol"],
                        st.session_state["coo_nucleo"],
                        filtrar_metales=f_met, filtrar_sales=f_sal, filtrar_cargas=f_car,
                        filtrar_nucleo=f_nuc,
                        mw_relativo=(mw_rel if f_mw else None),
                        mw_padre=(st.session_state["coo_mwpadre"] if f_mw else None),
                        log=lambda m: s.write(m))
                    st.session_state["coo_filt"] = dff
                    st.session_state["coo_conteos"] = ct
                    for k in ("coo_clas", "coo_sel"):
                        st.session_state.pop(k, None)
                    s.update(label=f"Quedan {len(dff)}", state="complete")
            if _hecho("coo_filt"):
                ct = st.session_state["coo_conteos"]
                a, b, c = st.columns(3)
                a.metric("Antes", ct.get("inicial", 0))
                b.metric("Después", ct.get("final", 0))
                c.metric("Descartados", ct.get("inicial", 0) - ct.get("final", 0))
                desc = {k: v for k, v in ct.items() if k not in ("inicial", "final")}
                if desc:
                    st.bar_chart(pd.Series(desc))
                st.success("Continúe en **③ Clasificar**.")

    # ── ③ Clasificar ──────────────────────────────────────────────────────────
    with tabs[2]:
        if not _hecho("coo_filt"):
            st.info("Complete el paso ②.")
        else:
            min_at = st.number_input("Mínimo de átomos por sustituyente", 1, 20, 1)
            if st.button("🏷️ Clasificar", type="primary"):
                with st.status("Clasificando…", expanded=True) as s:
                    dfc, stt = screening.clasificar_dataframe(
                        st.session_state["coo_filt"], st.session_state["coo_scol"],
                        st.session_state["coo_nucleo"], min_atomos_sub=int(min_at),
                        log=lambda m: s.write(m))
                    st.session_state["coo_clas"] = dfc
                    st.session_state["coo_statclas"] = stt
                    st.session_state.pop("coo_sel", None)
                    s.update(label=f"{len(dfc)} clasificados", state="complete")
            if _hecho("coo_clas"):
                stt = st.session_state["coo_statclas"]
                a, b = st.columns(2)
                a.metric("Clasificados", stt["n_total"])
                b.metric("Máx. sustituyentes", stt["max_subs"])
                if stt["gf_counts"]:
                    st.caption("Grupos funcionales más frecuentes:")
                    st.bar_chart(pd.Series(dict(list(stt["gf_counts"].items())[:12])))
                cols = [c for c in st.session_state["coo_clas"].columns
                        if c in ("CID", "n_sust") or c.startswith("sub_1_")]
                st.dataframe(st.session_state["coo_clas"][cols].head(12), width="stretch")
                st.success("Continúe en **④ Seleccionar**.")

    # ── ④ Seleccionar ─────────────────────────────────────────────────────────
    with tabs[3]:
        if not _hecho("coo_clas"):
            st.info("Complete el paso ③.")
        else:
            max_mol = st.slider("Máx. de candidatos", 10, 300, 150, 10)
            with st.expander("⚙️ Opciones avanzadas"):
                c1, c2 = st.columns(2)
                min_gf = c1.number_input("Mínimo por grupo×posición", 1, 50, 2)
                max_gf = c2.number_input("Máximo por grupo×posición", 1, 50, 8)
                balance = st.checkbox("Balance donador/aceptor (ED/EW)", True)
                pct = st.slider("Mínimo de cada tipo", 0.0, 0.5, 0.20, 0.05)
            if st.button("🎯 Seleccionar", type="primary"):
                with st.status("Seleccionando…", expanded=True) as s:
                    dfs, st4 = screening.seleccionar_curada(
                        st.session_state["coo_clas"], max_por_mol=int(max_mol),
                        min_por_gf=int(min_gf), max_por_gf=int(max_gf),
                        balance_ed_ew=balance, pct_min_ed_ew=pct,
                        log=lambda m: s.write(m))
                    st.session_state["coo_sel"] = dfs
                    st.session_state["coo_statsel"] = st4
                    s.update(label=f"{st4['n_sel']} candidatos", state="complete")
            if _hecho("coo_sel"):
                st4 = st.session_state["coo_statsel"]
                a, b, c, d = st.columns(4)
                a.metric("Seleccionados", st4["n_sel"])
                b.metric("Donadores (ED)", st4["ed"])
                c.metric("Aceptores (EW)", st4["ew"])
                d.metric("Neutros", st4["neutro"])
                dfsel = st.session_state["coo_sel"]
                cols = [c for c in dfsel.columns
                        if c in ("CID", "n_sust", "ED_EW") or c.startswith("sub_1_")]
                st.dataframe(dfsel[cols].head(20), width="stretch")
                st.success("Continúe en **⑤ Exportar** o pase a **Geometrías 3D**.")

    # ── 💊 ADMET ──────────────────────────────────────────────────────────────
    # Va en pestaña propia y con selector de etapa a propósito: así el perfil
    # farmacológico puede mirarse en CUALQUIER punto del cribado —sobre los
    # derivados en bruto, sobre la lista depurada o sobre los candidatos— y no
    # solo al final.
    with tabs[4]:
        # Lo que la pestaña ⑤ exporta es CONSECUENCIA de lo que hay aquí en
        # pantalla, no una copia paralela con vida propia: se rehace en cada
        # pasada y se borra en cuanto deja de haber un perfil calculado. Antes
        # se guardaba y no se borraba nunca, de modo que al cambiar de molécula
        # Exportar seguía ofreciendo —por defecto, y a veces como única opción—
        # el perfil de la molécula anterior.
        st.session_state.pop("coo_admet", None)
        st.session_state.pop("coo_admet_de", None)

        etapas = [("Candidatos seleccionados (④)", "coo_sel"),
                  ("Clasificados (③)", "coo_clas"),
                  ("Lista depurada (②)", "coo_filt"),
                  ("Todos los derivados (①)", "coo_df")]
        hechas = [(e, k) for e, k in etapas if _hecho(k)]
        if not hechas:
            st.info("Complete al menos el paso ① para poder perfilar las moléculas.")
        else:
            et_a = st.selectbox("¿Sobre qué etapa del cribado?",
                                [e for e, _ in hechas], key="coo_admet_etapa")
            clave_a = dict((e, k) for e, k in hechas)[et_a]
            df_a = st.session_state[clave_a]
            scol_a = next((c for c in df_a.columns if "smiles" in c.lower()),
                          st.session_state.get("coo_scol"))
            if not scol_a or scol_a not in df_a.columns:
                st.warning("Esta etapa no conserva la columna de SMILES; use otra.")
            else:
                res_a = common.panel_admet(df_a, scol_a, key=f"coo_admet_{clave_a}")
                # Sin comprobar `len(res_a)`: que el filtro no deje ninguna
                # molécula es un resultado legítimo, y descartarlo dejaba en pie
                # el perfil anterior mientras la pantalla decía «0 de N».
                if res_a is not None and "Lipinski" in getattr(res_a, "columns", []):
                    st.session_state["coo_admet"] = res_a
                    st.session_state["coo_admet_de"] = et_a

    # ── ⑤ Exportar ────────────────────────────────────────────────────────────
    with tabs[5]:
        # El rótulo del perfil dice de qué etapa sale, para que no haya duda de
        # qué se está descargando.
        de = st.session_state.get("coo_admet_de")
        fuentes = [(f"Perfil ADMET 💊 · {de}" if de else "Perfil ADMET (💊)", "coo_admet"),
                   ("Candidatos seleccionados (④)", "coo_sel"),
                   ("Clasificados (③)", "coo_clas"),
                   ("Lista depurada (②)", "coo_filt"),
                   ("Todos los derivados (①)", "coo_df")]
        disp = [(e, k) for e, k in fuentes if _hecho(k)]
        if not disp:
            st.info("Complete al menos el paso ①.")
        else:
            et = st.selectbox("¿Qué exportar?", [e for e, _ in disp])
            clave = dict((e, k) for e, k in disp)[et]
            con_img = st.checkbox("Incluir estructuras 2D (más lento)", False)
            df = st.session_state[clave]
            scol = next((c for c in df.columns if "smiles" in c.lower()),
                        st.session_state.get("coo_scol"))
            st.caption(f"{len(df):,} molécula(s) · {len(df.columns)} columnas.")
            if len(df) == 0:
                st.warning("El conjunto elegido está vacío: la descarga saldrá sin filas.")
            if con_img and (not scol or scol not in df.columns):
                st.warning("Este conjunto no trae columna de SMILES: el Excel saldrá "
                           "sin las estructuras 2D.")

            cA, cB = st.columns(2)
            # `data` como invocable: el CSV se arma al pulsar y no en cada rerun.
            cA.download_button("⬇️ CSV", data=lambda: df.to_csv(index=False),
                               file_name="coordinacion.csv", mime="text/csv",
                               width="stretch", on_click="ignore")
            with cB:
                # El .xlsx se guarda en sesión y su botón se pinta FUERA del
                # `if`: antes vivía dentro y desaparecía en el primer rerun,
                # incluido el que provocaba su propio clic. La firma ata el
                # archivo guardado a la selección con la que se generó, para
                # que nunca se descargue un Excel de otra cosa.
                firma_x = (clave, bool(con_img), common.firma_conjunto(df, scol))
                if st.button("📊 Generar Excel", width="stretch"):
                    buf = io.BytesIO()
                    with st.status("Generando…") as s:
                        excel.exportar_excel(df, scol, buf, con_imagenes=con_img,
                                             log=lambda m: s.write(m))
                        s.update(label="Listo", state="complete")
                    st.session_state["coo_xlsx"] = (firma_x, buf.getvalue())
                guardado = st.session_state.get("coo_xlsx")
                if guardado and guardado[0] == firma_x:
                    st.download_button("⬇️ Descargar Excel", data=guardado[1],
                                       file_name="coordinacion.xlsx",
                                       mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                                       width="stretch", on_click="ignore")
                elif guardado:
                    st.caption("El Excel guardado es de otra selección u otras "
                               "opciones: vuelva a generarlo.")
