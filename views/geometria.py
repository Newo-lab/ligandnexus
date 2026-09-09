"""
Vista: Geometrías 3D y archivos de cálculo
==========================================
Genera geometrías 3D (pre-optimización con campo de fuerzas) para una molécula
o un conjunto, las muestra en un visor 3D interactivo y produce los inputs de
Gaussian / ORCA / Psi4.
"""

import pandas as pd
import streamlit as st

from ligandnexus import pubchem, geometry, chem
from . import common


def _origen_conjunto():
    """
    DataFrame de moléculas a procesar según lo que haya en sesión.

    Incluye los conjuntos ya filtrados por ADMET: quien acote por Lipinski o por
    barrera hematoencefálica quiere llevarse ESOS candidatos a 3D, no la lista
    entera. Los rótulos dicen de qué etapa sale cada uno.
    """
    de_coo = st.session_state.get("coo_admet_de")
    fuentes = [(f"Perfil ADMET 💊 de coordinación · {de_coo}" if de_coo
                else "Perfil ADMET 💊 de coordinación", "coo_admet"),
               ("Candidatos de coordinación (④)", "coo_sel"),
               ("Clasificados de coordinación (③)", "coo_clas"),
               ("Perfil ADMET 💊 del cribado general", "gen_perfil"),
               ("Derivados del cribado general", "gen_df")]
    # Un conjunto vacío (p. ej. un filtro ADMET que no dejó nada) no sirve para
    # generar geometrías: se ofrece solo lo que tiene filas.
    disp = []
    for e, k in fuentes:
        df = st.session_state.get(k)
        if df is not None and len(df):
            disp.append((e, k))
    if not disp:
        return None, None
    et = st.selectbox("Conjunto", [e for e, _ in disp], key="geo_fuente")
    return st.session_state[dict((e, k) for e, k in disp)[et]], et


def render():
    st.header("🧬 Geometrías 3D y archivos de cálculo")
    st.write("Genera geometrías 3D **pre-optimizadas con campo de fuerzas** (MMFF94, "
             "o UFF si falta parametrización) como punto de partida — **no** una "
             "optimización DFT — y arma los inputs de Gaussian, ORCA o Psi4.")

    origen = st.radio("¿Qué procesar?", ["Una sola molécula", "Un conjunto del cribado"],
                      horizontal=True, key="geo_origen")

    df_geo = None
    if origen == "Una sola molécula":
        modo = st.radio("Entrada", ["CID/Nombre", "SMILES"], horizontal=True, key="geo_modo1")
        c1, c2 = st.columns([3, 1])
        with c1:
            entrada = st.text_input("Molécula", key="geo_uno",
                                    placeholder="Ej.: 1474  ·  o  ·  c1ccncc1",
                                    label_visibility="collapsed")
        with c2:
            cargar = st.button("Cargar", type="primary", width="stretch")
        if cargar and entrada.strip():
            s = entrada.strip()
            try:
                if modo == "SMILES" or (not s.isdigit() and chem.mol_from_smiles(s)):
                    if chem.mol_from_smiles(s) is None:
                        st.error("SMILES no válido.")
                    else:
                        st.session_state["geo_df"] = pd.DataFrame(
                            [{"CID": "manual", "smiles": s, "abrev": "mol"}])
                        st.session_state["geo_nombre"] = "(SMILES)"
                else:
                    info = common.c_info_molecula(int(s))
                    st.session_state["geo_df"] = pd.DataFrame(
                        [{"CID": int(s), "smiles": info["SMILES"], "abrev": "mol"}])
                    st.session_state["geo_nombre"] = info.get("Title") or f"CID {s}"
            except Exception as e:
                st.error(f"No se pudo cargar: {e}")
        df_geo = st.session_state.get("geo_df")
        if df_geo is not None:
            smi0 = str(df_geo.iloc[0]["smiles"])
            u1, u2 = st.columns([1, 3])
            im = common.img_mol(smi0, (220, 170))
            if im:
                u1.image(im)
            u2.markdown(f"**Molécula:** {st.session_state.get('geo_nombre','')}")
            u2.code(smi0, language=None)
    else:
        df_geo, _ = _origen_conjunto()
        if df_geo is None:
            st.info("No hay conjuntos cargados. Use el cribado general o el análisis "
                    "de coordinación, o elija «Una sola molécula».")

    if df_geo is None:
        return

    st.divider()
    gA, gB, gC = st.columns(3)
    with gA:
        modo3d = st.radio("Calidad", ["Rápido (1 confórmero)", "Mejor energía (varios)"],
                          key="geo_calidad")
    with gB:
        limite = st.number_input("Máx. de moléculas", 1, int(len(df_geo)),
                                 min(int(len(df_geo)), 50), key="geo_lim")
    with gC:
        mult = st.number_input("Multiplicidad", 1, 10, 1, key="geo_mult",
                               help="La carga se calcula automáticamente.")
    with st.expander("⚙️ Opciones avanzadas"):
        n_confs = st.number_input("Confórmeros a explorar (modo mejor energía)", 1, 50, 10,
                                  key="geo_nconf")
        semilla = st.number_input("Semilla aleatoria", 0, 99999, 42, key="geo_seed")

    formatos = st.multiselect(
        "Formatos", [".xyz", ".sdf", ".mol", "Gaussian (.gjf)", "ORCA (.inp)", "Psi4 (.in)"],
        default=[".xyz", "Gaussian (.gjf)"], key="geo_fmt")

    # Plantillas por motor
    tpl = {}
    if "Gaussian (.gjf)" in formatos:
        with st.expander("📝 Plantilla Gaussian"):
            a, b, c = st.columns(3)
            g_mem = a.text_input("%mem", "8GB", key="g_mem")
            g_np = b.number_input("%nprocshared", 1, 256, 8, key="g_np")
            g_chk = c.checkbox("Incluir %chk", True, key="g_chk")
            g_met = st.text_input("Método", "M06-2X", key="g_met")
            g_bas = st.text_input("Base", "6-311+G(d,p)", key="g_bas")
            g_kw = st.text_input("Palabras clave", "", key="g_kw",
                                 help="P. ej. 'opt freq' o 'SCRF=(SMD,solvent=water)'.")
            g_ex = st.text_area("Secciones extra", "", key="g_ex")
            tpl["gaussian"] = dict(mem=g_mem, nprocshared=g_np, usar_chk=g_chk,
                                   metodo=g_met, base=g_bas, keywords=g_kw, secciones_extra=g_ex)
    if "ORCA (.inp)" in formatos:
        with st.expander("📝 Plantilla ORCA"):
            a, b, c = st.columns(3)
            o_np = a.number_input("nprocs", 1, 256, 8, key="o_np")
            o_mc = b.number_input("maxcore (MB)", 256, 64000, 3000, key="o_mc")
            o_kw = c.text_input("Palabras clave (!)", "TightSCF", key="o_kw")
            o_met = st.text_input("Método", "M06-2X", key="o_met")
            o_bas = st.text_input("Base", "def2-TZVP", key="o_bas")
            o_ex = st.text_area("Bloques extra", "", key="o_ex")
            tpl["orca"] = dict(metodo=o_met, base=o_bas, keywords=o_kw, nprocs=o_np,
                               maxcore=o_mc, bloques_extra=o_ex)
    if "Psi4 (.in)" in formatos:
        with st.expander("📝 Plantilla Psi4"):
            a, b = st.columns(2)
            p_mem = a.text_input("memory", "8 GB", key="p_mem")
            p_met = b.text_input("Método", "m06-2x", key="p_met")
            p_bas = st.text_input("Base", "6-311+G(d,p)", key="p_bas")
            p_ex = st.text_area("Opciones extra", "", key="p_ex")
            tpl["psi4"] = dict(memoria=p_mem, metodo=p_met, base=p_bas, extra=p_ex)

    if st.button("🧬 Generar geometrías", type="primary", key="geo_go"):
        if not formatos:
            st.warning("Seleccione al menos un formato.")
        else:
            scol = next((c for c in df_geo.columns if "smiles" in c.lower()), None)
            nconf = int(n_confs) if modo3d.startswith("Mejor") else 1
            sub = df_geo.head(int(limite))
            archivos, manifiesto, primer_pdb = {}, [], None
            n_ok = n_fail = 0
            total = len(sub)
            prog = st.progress(0.0, text="Generando…")
            for i, (_, row) in enumerate(sub.iterrows(), 1):
                smi = str(row.get(scol, ""))
                cid = row.get("CID", i)
                abrev = str(row.get("abrev") or row.get("cabeza_serie") or "mol")
                nombre = f"{abrev}_{cid}"
                mol, conf, metodo, ener, ncf = geometry.generar_3d(
                    smi, n_confs=nconf, semilla=int(semilla))
                if mol is None:
                    n_fail += 1
                    manifiesto.append({"CID": cid, "nombre": nombre, "estado": metodo})
                    prog.progress(i / total, text=f"{i}/{total}")
                    continue
                carga = geometry.carga_formal(mol)
                xyzl = geometry.xyz_lineas(mol, conf)
                if primer_pdb is None:
                    primer_pdb = geometry.pdb_block(mol, conf)
                    st.session_state["geo_primer_nombre"] = nombre
                if ".xyz" in formatos:
                    archivos[f"xyz/{nombre}.xyz"] = geometry.mol_a_xyz(mol, conf, f"{nombre} SMILES={smi}")
                if ".sdf" in formatos:
                    archivos[f"sdf/{nombre}.sdf"] = geometry.mol_a_sdf(mol, conf)
                if ".mol" in formatos:
                    archivos[f"mol/{nombre}.mol"] = geometry.mol_a_molblock(mol, conf)
                if "gaussian" in tpl:
                    archivos[f"gaussian/{nombre}.gjf"] = geometry.construir_input_gaussian(
                        nombre, carga, int(mult), xyzl, **tpl["gaussian"])
                if "orca" in tpl:
                    archivos[f"orca/{nombre}.inp"] = geometry.construir_input_orca(
                        nombre, carga, int(mult), xyzl, **tpl["orca"])
                if "psi4" in tpl:
                    archivos[f"psi4/{nombre}.in"] = geometry.construir_input_psi4(
                        nombre, carga, int(mult), xyzl, **tpl["psi4"])
                manifiesto.append({"CID": cid, "nombre": nombre, "estado": "OK",
                                   "n_atomos": mol.GetNumAtoms(), "carga": carga,
                                   "multiplicidad": int(mult), "metodo": metodo,
                                   "E_kcalmol": round(ener, 3) if ener is not None else "",
                                   "conformeros": ncf})
                n_ok += 1
                prog.progress(i / total, text=f"{i}/{total}")
            prog.empty()

            carpetas = {}
            for ruta, cont in archivos.items():
                carp, nom = ruta.split("/", 1)
                carpetas.setdefault(carp, {})[nom] = cont
            st.session_state["geo_zips"] = {c: common.zip_de(fs) for c, fs in carpetas.items()}
            st.session_state["geo_manif"] = pd.DataFrame(manifiesto)
            st.session_state["geo_res"] = {"ok": n_ok, "fail": n_fail}
            st.session_state["geo_pdb"] = primer_pdb

    # Resultados
    if "geo_zips" in st.session_state:
        r = st.session_state["geo_res"]
        st.success(f"Geometrías generadas: {r['ok']} · fallidas: {r['fail']}")

        if st.session_state.get("geo_pdb"):
            st.markdown(f"#### Vista 3D — {st.session_state.get('geo_primer_nombre','')}")
            estilo = st.radio("Estilo", ["stick", "ball"], horizontal=True, key="geo_estilo")
            common.mostrar_3d(st.session_state["geo_pdb"], height=420, estilo=estilo)

        st.dataframe(st.session_state["geo_manif"], width="stretch")
        st.download_button("⬇️ Manifiesto (CSV)",
                           data=st.session_state["geo_manif"].to_csv(index=False),
                           file_name="manifiesto_3D.csv", mime="text/csv", key="geo_dlman")
        if st.session_state["geo_zips"]:
            st.caption("Descargue cada formato:")
            etq = {"xyz": "⬇️ XYZ", "sdf": "⬇️ SDF", "mol": "⬇️ MOL",
                   "gaussian": "⬇️ Gaussian", "orca": "⬇️ ORCA", "psi4": "⬇️ Psi4"}
            zips = st.session_state["geo_zips"]
            cols = st.columns(min(3, len(zips)))
            for j, (carp, data) in enumerate(zips.items()):
                cols[j % len(cols)].download_button(
                    etq.get(carp, carp), data=data, file_name=f"{carp}.zip",
                    mime="application/zip", key=f"geo_dl_{carp}")
