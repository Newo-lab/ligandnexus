"""
common.py — Utilidades compartidas por las vistas de la interfaz
================================================================
Tema visual, envoltorios con caché para PubChem (evita repetir descargas al
re-ejecutarse la app), render 2D/3D y datos de citación.
"""

from __future__ import annotations

import io
import re
import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw

from ligandnexus import pubchem, chem, geometry, admet
from ligandnexus import sources as _sources

# ── Paleta ────────────────────────────────────────────────────────────────────
# El tema (claro/oscuro) lo fija Streamlit desde .streamlit/config.toml y con el
# menú ☰ → Settings → Theme. Los colores de aquí son ADAPTABLES: azul de acento
# que contrasta en ambos fondos y tarjetas translúcidas (gris con transparencia)
# que se ven bien tanto en oscuro como en claro.
AZUL = "#4a9fd4"        # acento (headings) — legible en claro y oscuro
COBRE = "#d98a3d"
VERDE = "#2E8B57"
AMBAR = "#D98C1A"
ROJO = "#C24A4A"

CSS = f"""
<style>
.block-container {{ padding-top: 2rem; max-width: 1180px; }}
h1, h2, h3 {{ color: {AZUL}; }}
[data-testid="stMetricValue"] {{ font-size: 1.5rem; color: {AZUL}; }}
.lf-card {{
    border: 1px solid rgba(127,127,127,0.22); border-radius: 12px;
    padding: 1rem 1.2rem; background: rgba(127,127,127,0.10);
    margin-bottom: .6rem;
}}
.lf-card h3 {{ margin-top: 0; }}
.lf-pill {{
    display:inline-block; padding: 2px 10px; border-radius: 999px;
    font-size: .78rem; font-weight: 600; color: white;
}}
.stButton>button {{ border-radius: 8px; }}
</style>
"""

def es_oscuro():
    """True si el tema activo de Streamlit es oscuro (por defecto asume oscuro)."""
    try:
        tema = getattr(getattr(st, "context", None), "theme", None)
        return getattr(tema, "type", "dark") != "light"
    except Exception:
        return True


def _fondo_3d():
    return "0x0e1117" if es_oscuro() else "0xffffff"


def aplicar_tema():
    st.markdown(CSS, unsafe_allow_html=True)


# ── PubChem con caché (10 min) ────────────────────────────────────────────────

@st.cache_data(show_spinner=False, ttl=600)
def c_info_molecula(cid):
    return pubchem.info_molecula(int(cid))


@st.cache_data(show_spinner=False, ttl=600)
def c_sinonimos(cid, n=6):
    return pubchem.sinonimos(int(cid), n=n)


@st.cache_data(show_spinner=False, ttl=600)
def c_buscar_cid_por_nombre(nombre):
    return pubchem.buscar_cid_por_nombre(nombre)


# ── Render 2D / 3D ────────────────────────────────────────────────────────────

def img_mol(smiles, size=(320, 240)):
    """Estructura 2D adaptada al tema (fondo transparente, trazos claros en oscuro)."""
    from io import BytesIO
    from PIL import Image
    from rdkit.Chem.Draw import rdMolDraw2D

    mol = Chem.MolFromSmiles(str(smiles)) if smiles else None
    if mol is None:
        return None
    d = rdMolDraw2D.MolDraw2DCairo(*size)
    if es_oscuro():
        rdMolDraw2D.SetDarkMode(d)
    d.drawOptions().setBackgroundColour((0, 0, 0, 0))   # transparente
    rdMolDraw2D.PrepareAndDrawMolecule(d, mol)
    d.FinishDrawing()
    return Image.open(BytesIO(d.GetDrawingText()))


# Solo se dejan pasar los registros de coordenadas, que RDKit escribe en
# columnas de ancho fijo. El titulo COMPND, en cambio, se copia literalmente
# desde la propiedad _Name, y en un SDF esa propiedad la escribe QUIEN HIZO
# EL ARCHIVO: es texto arbitrario que py3Dmol interpola sin escapar dentro
# del HTML del visor. Un titulo con <script> se ejecutaria en la pagina.
_PDB_COORD = re.compile(r"^(ATOM  |HETATM|CONECT|MODEL |ENDMDL|TER   |END)")


def _limpiar_pdb(bloque: str) -> str:
    """Descarta todo lo que no sean coordenadas de formato fijo."""
    return "\n".join(l for l in bloque.splitlines() if _PDB_COORD.match(l))


def mostrar_3d(pdb_block, height=420, estilo="stick"):
    """Visor 3D interactivo (py3Dmol) embebido en la página."""
    import py3Dmol
    view = py3Dmol.view(width="100%", height=height)
    view.addModel(_limpiar_pdb(pdb_block), "pdb")
    if estilo == "stick":
        view.setStyle({"stick": {"radius": 0.14}, "sphere": {"scale": 0.22}})
    elif estilo == "ball":
        view.setStyle({"stick": {"radius": 0.12}, "sphere": {"scale": 0.30}})
    else:
        view.setStyle({estilo: {}})
    view.setBackgroundColor(_fondo_3d())
    view.zoomTo()
    # st.components.v1.html quedó obsoleto (retirado tras 2026-06-01); st.iframe
    # acepta HTML crudo y por defecto ocupa todo el ancho, igual que antes.
    st.iframe(view._make_html(), height=height + 15)


# ── Verificador de molécula (compartido por varias vistas) ────────────────────

def resolver_entrada(texto, modo):
    """
    Resuelve una entrada del usuario a un dict de info de molécula.
    modo: 'CID' o 'Nombre'. Devuelve (info|None, mensaje_error|None).
    """
    try:
        if modo == "Nombre":
            if not (texto or "").strip():
                return None, "Escriba un nombre para buscar."
            cids = c_buscar_cid_por_nombre(texto.strip())
            if not cids:
                return None, "PubChem no encontró ningún compuesto con ese nombre."
            cid = cids[0]
        else:
            cid = int(texto)
        info = c_info_molecula(cid)
        info["sinonimos"] = c_sinonimos(cid, 6)
        if not info.get("SMILES"):
            return None, "PubChem no devolvió estructura para ese compuesto."
        return info, None
    except Exception as e:
        return None, f"No se pudo completar la consulta: {e}"


TIPOS_ARCHIVO = ["sdf", "sd", "mol", "smi", "txt", "csv"]


def selector_fuente(key):
    """
    Selector de base de datos (registro de ~20). Devuelve (entry, archivo_subido).
    Muestra el uploader/enlace de descarga según el modo de la fuente.
    """
    reg = _sources.REGISTRO
    # Agrupar etiquetas por categoría para el desplegable.
    etiquetas = {}
    for e in reg:
        icono = {"api": "🟢", "local": "📁", "descarga": "⬇️"}.get(e["modo"], "•")
        etiquetas[f"{icono} {e['nombre']} — {e['cat']}"] = e["id"]
    labels = list(etiquetas.keys())
    sel = st.selectbox("Base de datos", labels, index=0, key=f"{key}_src",
                       help="🟢 búsqueda en vivo · ⬇️ descarga su volcado y búscalo aquí · "
                            "📁 tu propio archivo")
    entry = _sources.por_id(etiquetas[sel])
    up = None
    if entry["modo"] == "api":
        st.caption(f"🟢 Búsqueda en vivo. {entry['desc']}")
    elif entry["modo"] == "local":
        st.caption(entry["desc"])
        up = st.file_uploader("Sube tu archivo (SDF · SMILES · CSV)",
                              type=TIPOS_ARCHIVO, key=f"{key}_up")
    else:  # descarga
        st.info(f"**{entry['nombre']}** no ofrece búsqueda de subestructura por API. "
                f"Descarga su volcado (una sola vez) y súbelo aquí: la búsqueda de "
                f"subestructura se hace localmente con RDKit.\n\n_{entry['desc']}_")
        if entry.get("url"):
            st.markdown(f"[⬇️ Ir a la descarga de {entry['nombre']}]({entry['url']})")
        up = st.file_uploader("Sube el volcado descargado (SDF · SMILES · CSV)",
                              type=TIPOS_ARCHIVO, key=f"{key}_up")
    return entry, up


def buscar_derivados_fuente(entry, smiles_core, max_records, log, subido=None):
    """Ejecuta la búsqueda de derivados en la fuente elegida (API o archivo local)."""
    if entry["modo"] == "api":
        return entry["src"].buscar(smiles_core, max_records=max_records, log=log)
    if subido is None:
        raise ValueError("Sube primero un archivo (SDF/SMILES/CSV) para esta fuente.")
    return _sources.LOCAL.buscar_local(smiles_core, subido.getvalue(), subido.name,
                                       max_records=max_records, log=log)


def zip_de(archivos: dict) -> bytes:
    import zipfile
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for nombre, contenido in archivos.items():
            z.writestr(nombre, contenido)
    return buf.getvalue()


# ── Moléculas de ejemplo ──────────────────────────────────────────────────────
CABEZAS = {
    "2,2'-Bipiridina (bpy)": 1474,
    "1,10-Fenantrolina (phen)": 1318,
    "2-(2-Piridil)imidazol (pyim)": 589313,
    "Di(2-piridil)amina (dipa)": 14547,
    "Etilendiamina (en)": 3301,
    "2,2':6',2''-Terpiridina (terpy)": 70848,
    "Clioquinol (CQ)": 2788,
    "Glicina (Gly)": 750,
    "Ácido iminodiacético (IDA)": 8897,
    "Tris(2-piridilmetil)amina (TPMA)": 379259,
}


# ── Referencias (ACS / APA) ───────────────────────────────────────────────────
REFERENCIAS = [
    {"tool": "RDKit", "doc": "https://www.rdkit.org/docs/",
     "acs": "Landrum, G. RDKit: Open-Source Cheminformatics Software. "
            "https://www.rdkit.org (accedido 2026). DOI: 10.5281/zenodo.591637.",
     "apa": "Landrum, G. (s. f.). RDKit: Open-source cheminformatics. https://www.rdkit.org"},
    {"tool": "PubChem", "doc": "https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest",
     "acs": "Kim, S.; et al. PubChem 2025 Update. Nucleic Acids Res. 2025, 53 (D1), "
            "D1516–D1525. DOI: 10.1093/nar/gkae1059.",
     "apa": "Kim, S., et al. (2025). PubChem 2025 update. Nucleic Acids Research, "
            "53(D1), D1516–D1525. https://doi.org/10.1093/nar/gkae1059"},
    {"tool": "Streamlit", "doc": "https://docs.streamlit.io/",
     "acs": "Streamlit Inc. Streamlit: A Faster Way to Build and Share Data Apps. "
            "https://streamlit.io (accedido 2026).",
     "apa": "Streamlit Inc. (2019). Streamlit [Software]. https://streamlit.io"},
    {"tool": "pandas", "doc": "https://pandas.pydata.org/docs/",
     "acs": "McKinney, W. Data Structures for Statistical Computing in Python. "
            "Proc. 9th Python in Science Conf. 2010, 51–56. DOI: 10.25080/Majora-92bf1922-00a.",
     "apa": "McKinney, W. (2010). Data structures for statistical computing in Python. "
            "Proceedings of the 9th Python in Science Conference, 51–56."},
    {"tool": "py3Dmol / 3Dmol.js", "doc": "https://3dmol.csb.pitt.edu/",
     "acs": "Rego, N.; Koes, D. 3Dmol.js: Molecular Visualization with WebGL. "
            "Bioinformatics 2015, 31 (8), 1322–1324. DOI: 10.1093/bioinformatics/btu829.",
     "apa": "Rego, N., & Koes, D. (2015). 3Dmol.js: molecular visualization with WebGL. "
            "Bioinformatics, 31(8), 1322–1324."},
    {"tool": "openpyxl / Pillow", "doc": "https://openpyxl.readthedocs.io/",
     "acs": "Gazoni, E.; Clark, C. openpyxl. https://openpyxl.readthedocs.io. "
            "Clark, A. Pillow (PIL Fork). https://pillow.readthedocs.io.",
     "apa": "Gazoni, E., & Clark, C. openpyxl [Software]. Clark, A. (2015). Pillow (PIL Fork)."},
]

REFERENCIAS_QM = [
    {"tool": "Gaussian 16", "doc": "https://gaussian.com/citation/",
     "acs": "Frisch, M. J.; et al. Gaussian 16, Revision C.01; Gaussian, Inc.: "
            "Wallingford, CT, 2016.",
     "apa": "Frisch, M. J., et al. (2016). Gaussian 16, Revision C.01 [Software]."},
    {"tool": "ORCA", "doc": "https://www.faccts.de/orca/",
     "acs": "Neese, F. Software Update: The ORCA Program System—Version 5.0. "
            "WIREs Comput. Mol. Sci. 2022, 12, e1606. DOI: 10.1002/wcms.1606.",
     "apa": "Neese, F. (2022). ORCA program system—Version 5.0. WIREs Comput. Mol. Sci., 12, e1606."},
    {"tool": "Psi4 1.4", "doc": "https://psicode.org/",
     "acs": "Smith, D. G. A.; et al. Psi4 1.4. J. Chem. Phys. 2020, 152, 184108. "
            "DOI: 10.1063/5.0006002.",
     "apa": "Smith, D. G. A., et al. (2020). Psi4 1.4. J. Chem. Phys., 152, 184108."},
]

# Reglas y descriptores del panel ADMET. Se citan aparte porque son criterios
# de la literatura, no software: quien publique usando el tamiz debe citarlos.
REFERENCIAS_ADMET = [
    {"tool": "Lipinski — regla de 5", "doc": "https://doi.org/10.1016/S0169-409X(00)00129-0",
     "acs": "Lipinski, C. A.; Lombardo, F.; Dominy, B. W.; Feeney, P. J. Experimental and "
            "Computational Approaches to Estimate Solubility and Permeability in Drug "
            "Discovery and Development Settings. Adv. Drug Deliv. Rev. 2001, 46, 3–26.",
     "apa": "Lipinski, C. A., Lombardo, F., Dominy, B. W., & Feeney, P. J. (2001). "
            "Adv. Drug Deliv. Rev., 46, 3–26."},
    {"tool": "Veber — biodisponibilidad oral", "doc": "https://doi.org/10.1021/jm020017n",
     "acs": "Veber, D. F.; et al. Molecular Properties That Influence the Oral "
            "Bioavailability of Drug Candidates. J. Med. Chem. 2002, 45, 2615–2623.",
     "apa": "Veber, D. F., et al. (2002). J. Med. Chem., 45, 2615–2623."},
    {"tool": "Egan — «huevo» de absorción", "doc": "https://doi.org/10.1021/jm000292e",
     "acs": "Egan, W. J.; Merz, K. M.; Baldwin, J. J. Prediction of Drug Absorption Using "
            "Multivariate Statistics. J. Med. Chem. 2000, 43, 3867–3877.",
     "apa": "Egan, W. J., Merz, K. M., & Baldwin, J. J. (2000). J. Med. Chem., 43, 3867–3877."},
    {"tool": "Clark — logBB (barrera hematoencefálica)",
     "doc": "https://doi.org/10.1021/js980402t",
     "acs": "Clark, D. E. In Silico Prediction of Blood–Brain Barrier Permeation. "
            "J. Pharm. Sci. 1999, 88, 815–821.",
     "apa": "Clark, D. E. (1999). J. Pharm. Sci., 88, 815–821."},
    {"tool": "QED — fármaco-similitud", "doc": "https://doi.org/10.1038/nchem.1243",
     "acs": "Bickerton, G. R.; Paolini, G. V.; Besnard, J.; Muresan, S.; Hopkins, A. L. "
            "Quantifying the Chemical Beauty of Drugs. Nat. Chem. 2012, 4, 90–98.",
     "apa": "Bickerton, G. R., et al. (2012). Nat. Chem., 4, 90–98."},
    {"tool": "SAscore — accesibilidad sintética",
     "doc": "https://doi.org/10.1186/1758-2946-1-8",
     "acs": "Ertl, P.; Schuffenhauer, A. Estimation of Synthetic Accessibility Score of "
            "Drug-like Molecules. J. Cheminform. 2009, 1, 8.",
     "apa": "Ertl, P., & Schuffenhauer, A. (2009). J. Cheminform., 1, 8."},
    {"tool": "PAINS — interferencia en ensayos",
     "doc": "https://doi.org/10.1021/jm901137j",
     "acs": "Baell, J. B.; Holloway, G. A. New Substructure Filters for Removal of Pan "
            "Assay Interference Compounds (PAINS). J. Med. Chem. 2010, 53, 2719–2740.",
     "apa": "Baell, J. B., & Holloway, G. A. (2010). J. Med. Chem., 53, 2719–2740."},
    {"tool": "TPSA — superficie polar topológica",
     "doc": "https://doi.org/10.1021/jm000942e",
     "acs": "Ertl, P.; Rohde, B.; Selzer, P. Fast Calculation of Molecular Polar Surface "
            "Area. J. Med. Chem. 2000, 43, 3714–3717.",
     "apa": "Ertl, P., Rohde, B., & Selzer, P. (2000). J. Med. Chem., 43, 3714–3717."},
]


# ── Panel ADMET reutilizable ──────────────────────────────────────────────────

def firma_conjunto(df, smiles_col=None):
    """
    Huella barata de un conjunto, para saber si sigue siendo el mismo sin
    compararlo entero. La usan el panel ADMET (para no recalcular de balde) y
    las exportaciones (para no entregar un archivo generado con otra cosa).

    Mira el número de filas y de columnas, y los 50 SMILES de cada extremo: dos
    filtros distintos pueden dejar el mismo número de moléculas, así que el
    recuento por sí solo no basta.
    """
    if df is None:
        return None
    if smiles_col and smiles_col in df.columns:
        col = df[smiles_col].astype(str)
    elif len(df.columns):
        col = df.iloc[:, 0].astype(str)
    else:
        return (len(df), 0)
    return (len(df), len(df.columns),
            hash(tuple(col.head(50))), hash(tuple(col.tail(50))))


def panel_admet(df, smiles_col, key="admet"):
    """
    Panel de perfil ADMET, pensado para engancharse en CUALQUIER punto de un
    cribado: recibe el conjunto que haya en pantalla, calcula los descriptores y
    las reglas de fármaco-similitud con `ligandnexus.admet`, y deja filtrar por
    ellas.

    Devuelve el DataFrame resultante: con las columnas ADMET añadidas si se han
    calculado (y filtrado, si se han activado filtros), o el original si no. Así
    la exportación posterior arrastra lo que el usuario esté viendo.
    """
    import pandas as pd

    if df is None or df.empty or not smiles_col or smiles_col not in df.columns:
        return df

    st.markdown("#### 💊 Perfil ADMET (opcional)")
    st.caption(
        "Descriptores y reglas de biodisponibilidad calculados **en local con "
        "RDKit** sobre las moléculas que hay ahora en pantalla. Es un tamiz de "
        "**priorización**, no una predicción farmacocinética: incumplir una regla "
        "señala una molécula, no la descarta.")

    k_df, k_firma = f"_{key}_df", f"_{key}_firma"
    k_pedido = f"_{key}_pedido"          # el usuario ya pidió el perfil aquí
    firma = firma_conjunto(df, smiles_col)
    if st.session_state.get(k_firma) != firma:
        # El conjunto cambió (otro filtro, otra búsqueda): lo ya calculado no vale.
        st.session_state.pop(k_df, None)

    # Si el perfil ya se pidió para este panel y TODAS las moléculas que quedan
    # están memorizadas, rehacerlo no cuesta nada: se rehace solo. Así acotar un
    # filtro no tira por la borda un cálculo de minutos ni obliga a volver a
    # pulsar el botón. Si aparecen moléculas nuevas (otra búsqueda), no: el
    # cálculo sería caro y debe pedirlo el usuario.
    gratis = admet.faltan_por_calcular(df[smiles_col].astype(str)) == 0
    rehacer = (st.session_state.get(k_pedido) and gratis
               and st.session_state.get(k_df) is None)

    if st.button(f"🧪 Calcular ADMET de {len(df):,} molécula(s)", key=f"{key}_calc") or rehacer:
        barra = st.progress(0.0, text="Calculando descriptores…")
        res = admet.evaluar_dataframe(
            df, smiles_col,
            progreso=lambda i, t: barra.progress(i / max(t, 1),
                                                 text=f"Calculando… {i:,}/{t:,}"))
        barra.empty()
        st.session_state[k_df] = res
        st.session_state[k_firma] = firma
        st.session_state[k_pedido] = True

    dfa = st.session_state.get(k_df)
    if dfa is None:
        return df

    r = admet.resumen(dfa)
    c = st.columns(6)
    c[0].metric("Lipinski", f"{r['Lipinski']}/{r['n']}",
                help="MW≤500, clogP≤5, HBD≤5, HBA≤10 — se admite una violación "
                     "(Lipinski 2001).")
    c[1].metric("Veber", f"{r['Veber']}/{r['n']}",
                help="Enlaces rotables ≤10 y TPSA ≤140 Å² (Veber 2002).")
    c[2].metric("Egan", f"{r['Egan']}/{r['n']}",
                help="TPSA ≤131,6 y clogP ≤5,88 (Egan 2000).")
    c[3].metric("Atraviesa BHE", r.get("BHE_atraviesa", 0),
                help="logBB > 0,3 por la ecuación de Clark (1999). Es una "
                     "correlación de dos descriptores: ordena bien una biblioteca, "
                     "pero no dictamina sobre una molécula concreta.")
    c[4].metric("Sintetizable", f"{r['SA_accesible']}/{r['n']}",
                help="SAscore ≤ 3,5 (Ertl y Schuffenhauer 2009).")
    c[5].metric("Alertas PAINS", r["PAINS"],
                help="Subestructuras que suelen dar falsos positivos en los "
                     "ensayos (Baell y Holloway 2010). Los polifenoles las "
                     "disparan a menudo: es un aviso, no una condena.")

    with st.expander("🔬 Filtrar por criterios ADMET"):
        f1, f2, f3, f4 = st.columns(4)
        solo_lip = f1.checkbox("Cumple Lipinski", key=f"{key}_lip")
        solo_veb = f2.checkbox("Cumple Veber", key=f"{key}_veb")
        solo_sa = f3.checkbox("SAscore ≤ 3,5", key=f"{key}_sa")
        sin_pains = f4.checkbox("Sin alertas PAINS", key=f"{key}_pains")
        bhe = st.multiselect(
            "Barrera hematoencefálica (logBB de Clark)",
            ["atraviesa", "intermedio", "no atraviesa"], key=f"{key}_bhe",
            help="Vacío = no filtrar por barrera.")

    out = dfa
    if solo_lip and "Lipinski" in out.columns:
        out = out[out["Lipinski"].fillna(False).astype(bool)]
    if solo_veb and "Veber" in out.columns:
        out = out[out["Veber"].fillna(False).astype(bool)]
    if solo_sa and "SA_accesible" in out.columns:
        out = out[out["SA_accesible"].fillna(False).astype(bool)]
    if sin_pains and "PAINS" in out.columns:
        out = out[~out["PAINS"].fillna(False).astype(bool)]
    if bhe and "BHE" in out.columns:
        out = out[out["BHE"].isin(bhe)]
    out = out.reset_index(drop=True)

    if len(out) != len(dfa):
        st.success(f"Filtro ADMET: **{len(out)}** de {len(dfa)} moléculas.")

    vista = [c_ for c_ in ("CID", smiles_col, "MW", "clogP", "TPSA", "TPSA_admet",
                           "HBD", "HBA", "RotB", "QED", "SAscore", "logBB", "BHE",
                           "Lipinski_viol", "PAINS_alerta") if c_ in out.columns]
    st.dataframe(out[vista] if vista else out, width="stretch", height=280)
    # `data` como invocable: Streamlit lo ejecuta al pulsar, no en cada rerun.
    # Serializar el CSV en cada pasada costaba 178 ms con 20.000 filas y las 28
    # columnas que deja ADMET (56 ms sin él), y `st.tabs` pinta todas.
    st.download_button("⬇️ Descargar tabla ADMET (CSV)",
                       data=lambda: out.to_csv(index=False),
                       file_name="admet.csv", mime="text/csv", key=f"{key}_csv",
                       on_click="ignore")

    st.caption(
        "⚠️ La **TPSA** de esta tabla es la de Ertl (solo N y O), que es la que "
        "espera la ecuación de Clark; si la base de datos ya traía su propia "
        "columna `TPSA`, la calculada aquí aparece como `TPSA_admet` porque no "
        "son la misma magnitud. Y ojo con la tensión entre criterios: **Clark "
        "premia el clogP alto y Lipinski lo castiga**, así que optimizar hacia el "
        "cerebro empuja en contra de la regla de 5 — conviene mirar las dos "
        "columnas juntas.")
    return out
