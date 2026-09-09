"""
screening.py — Depuración y selección curada (modo coordinación)
================================================================
Filtros de la biblioteca (metales, sales, cargas, núcleo, peso) y selección de
un subconjunto representativo con balance donador/aceptor (ED/EW) para estudios
de coordinación. La clasificación de sustituyentes vive en `chem.py`.
"""

from __future__ import annotations

import re

import pandas as pd
from rdkit import Chem

from . import chem

# Metales/semimetales para el filtro. Se descartan porque los derivados que ya
# traen un metal no sirven como ligandos libres del cribado.
METALES = [
    "Li", "Be", "Na", "Mg", "Al", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe",
    "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc",
    "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Cs", "Ba", "La", "Ce", "Pr",
    "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf",
    "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "Si",
    "B", "As", "Se", "Te", "At",
]


_METALES_SET = set(METALES)
# Una formula se parte en simbolos de elemento: mayuscula + minuscula opcional.
_RE_ELEMENTO = re.compile(r"[A-Z][a-z]?")


def _tiene_metal(formula) -> bool:
    """
    True si la formula contiene algun metal/semimetal.

    Se TOKENIZA la formula en simbolos de elemento en vez de buscar cada simbolo
    con una expresion regular. El `(?<![A-Z])` de la version anterior impedia que
    casara un simbolo precedido por otra mayuscula, y la notacion de Hill de
    PubChem produce justo eso (C21H16N2OSi, C10H8N2OSn, C6H4NOSe, C6H5NSb), asi
    que los compuestos de Si, Sn, Se y Sb se colaban por el filtro.
    """
    if not isinstance(formula, str):
        return False
    return any(el in _METALES_SET for el in _RE_ELEMENTO.findall(formula))


def _tiene_carga(smiles) -> bool:
    """
    True si la especie tiene carga NETA distinta de cero.

    Mirar si aparece un [X+] o [X-] en el SMILES no basta: hay grupos NEUTROS que
    se escriben con cargas formales separadas —**nitro** [N+](=O)[O-], N-oxidos,
    azidas, zwitteriones— y esa version los descartaba como si fueran iones.
    Medido sobre derivados reales de bipiridina de PubChem, tiraba el 2,3 % de la
    biblioteca, y no al azar: se llevaba justo los **nitroderivados**, el grupo
    electroatrayente de referencia que despues `seleccionar_curada` intenta
    balancear (Nitro esta en EW_GRUPOS). Ahora se suma la carga formal real.

    El tamiz por texto se conserva como atajo: si no hay ninguna carga formal
    escrita, la molecula es neutra y no hace falta parsearla.
    """
    s = str(smiles)
    if not re.search(r"\[[^\]]*[+\-][^\]]*\]", s):
        return False
    m = chem.mol_from_smiles(s)
    if m is None:
        return True                     # ilegible: se descarta por precaucion
    return Chem.GetFormalCharge(m) != 0


def aplicar_filtros(df, smiles_col, smarts_nucleo,
                    filtrar_metales=True, filtrar_sales=True,
                    filtrar_cargas=True, filtrar_nucleo=True,
                    mw_relativo=None, mw_padre=None, log=print):
    """Aplica los filtros (cada uno activable). Devuelve (df_filtrado, conteos)."""
    conteos = {}
    n_ini = len(df)

    if filtrar_metales and "MolecularFormula" in df.columns:
        mask = df["MolecularFormula"].apply(_tiene_metal)
        conteos["metal"] = int(mask.sum())
        df = df[~mask].copy()

    if filtrar_sales:
        mask = df[smiles_col].astype(str).str.contains(".", regex=False)
        conteos["sal_mezcla"] = int(mask.sum())
        df = df[~mask].copy()

    if filtrar_cargas:
        mask = df[smiles_col].apply(_tiene_carga)
        conteos["carga_formal"] = int(mask.sum())
        df = df[~mask].copy()

    if filtrar_nucleo and smarts_nucleo:
        core = chem.query_nucleo(smarts_nucleo)
        if core is not None:
            mask = df[smiles_col].apply(lambda s: chem.contiene_nucleo(s, core))
            conteos["nucleo_invalido"] = int((~mask).sum())
            df = df[mask].copy()

    if mw_relativo is not None and mw_padre is not None and "MolecularWeight" in df.columns:
        techo = float(mw_padre) + float(mw_relativo)
        mw = pd.to_numeric(df["MolecularWeight"], errors="coerce")
        mask = mw > techo
        conteos["mw_excede"] = int(mask.sum())
        df = df[~mask].copy()

    df = df.reset_index(drop=True)
    conteos["inicial"] = n_ini
    conteos["final"] = len(df)
    log(f"  Filtros: {n_ini:,} -> {len(df):,}")
    return df, conteos


# ── Clasificación de todo el DataFrame ────────────────────────────────────────

_COLS_BASE = ["CID", "IUPACName", "MolecularFormula", "MolecularWeight",
              "HBondDonorCount", "HBondAcceptorCount", "XLogP", "TPSA"]


def clasificar_dataframe(df, smiles_col, smarts_nucleo, min_atomos_sub=1,
                         gf_patrones=None, log=print):
    """Clasifica los sustituyentes de cada fila. Devuelve (df_clasificado, stats)."""
    gf_comp = chem.compilar_gf(gf_patrones) if gf_patrones else chem._GF_SMARTS_MOL

    registros, n_sin_match, max_subs = [], 0, 0
    total = len(df)
    for ri, (_, row) in enumerate(df.iterrows(), 1):
        smi = str(row.get(smiles_col, ""))
        subs = chem.extraer_sustituyentes(smi, smarts_nucleo, gf_comp)
        if not subs:
            n_sin_match += 1
        subs = [s for s in subs if s["n_atomos_sub"] >= min_atomos_sub and s["gf"] != "H"]
        subs.sort(key=lambda s: (s["dist_coord"], s["n_atomos_sub"]))
        max_subs = max(max_subs, len(subs))

        rec = {c: row.get(c, "") for c in _COLS_BASE}
        rec["smiles"] = smi
        rec["n_sust"] = len(subs)
        for i, s in enumerate(subs, 1):
            rec[f"sub_{i}_gf"] = s["gf"]
            rec[f"sub_{i}_pos"] = s["pos_label"]
            rec[f"sub_{i}_dist"] = s["dist_coord"]
            rec[f"sub_{i}_smiles"] = s["sub_smiles"]
        registros.append(rec)
        if ri % 200 == 0:
            log(f"  Clasificando: {ri:,}/{total:,}")

    for rec in registros:
        for i in range(1, max_subs + 1):
            for key in (f"sub_{i}_gf", f"sub_{i}_pos", f"sub_{i}_dist", f"sub_{i}_smiles"):
                rec.setdefault(key, "")

    df_result = pd.DataFrame(registros)

    gf_counts = {}
    for i in range(1, max_subs + 1):
        col = f"sub_{i}_gf"
        if col in df_result.columns:
            for gf in df_result[col]:
                if gf:
                    gf_counts[gf] = gf_counts.get(gf, 0) + 1

    stats = {
        "n_total": len(df_result), "max_subs": max_subs, "sin_match": n_sin_match,
        "gf_counts": dict(sorted(gf_counts.items(), key=lambda x: -x[1])),
    }
    log(f"  Clasificados: {len(df_result):,} | sin núcleo: {n_sin_match} | "
        f"máx sust.: {max_subs}")
    return df_result, stats


# ── Selección curada con balance ED/EW ────────────────────────────────────────

ED_GRUPOS = {
    "Metilo", "Etilo", "Propilo", "Alquilo", "Cicloalquilo", "Metileno",
    "Metoxi", "Alcoxi", "Hidroxi", "Amino-NH2", "Amino-NH", "Amino-N", "Arilo",
}
EW_GRUPOS = {
    "CF3", "CHF2", "Halo-F", "Halo-Cl", "Halo-Br", "Halo-I", "Nitro", "Cyano",
    "Isociano", "Sulfonyl", "Sulfinyl", "Carboxilo", "Ester", "Amida", "Cetona",
    "Aldehido", "Carbonilo", "Trifluoromethyl",
}


def clase_ed_ew(gf):
    if gf in ED_GRUPOS: return "ED"
    if gf in EW_GRUPOS: return "EW"
    return "N"


def seleccionar_curada(df_clas, max_por_mol=150, min_por_gf=2, max_por_gf=8,
                       pref_n_sust=(1, 2, 3), balance_ed_ew=True,
                       pct_min_ed_ew=0.20, log=print):
    """Selección representativa (variedad GF×posición + balance ED/EW)."""
    vacio = {"n_sel": 0, "ed": 0, "ew": 0, "neutro": 0, "gf_top": [], "pos_top": []}
    if df_clas is None or df_clas.empty:
        return pd.DataFrame(), vacio

    sub_cols = [c for c in df_clas.columns if c.startswith("sub_") and c.endswith("_gf")]
    if not sub_cols:
        log("  Sin columnas de sustituyentes — nada que seleccionar.")
        return pd.DataFrame(), vacio

    df = df_clas.copy()
    df["n_sust"] = pd.to_numeric(df.get("n_sust", 0), errors="coerce").fillna(0).astype(int)
    df["_gf1"] = df.get("sub_1_gf", "").fillna("").astype(str)
    df["_pos1"] = df.get("sub_1_pos", "").fillna("").astype(str)
    df["_ed_ew"] = df["_gf1"].apply(clase_ed_ew)

    orden = {n: i for i, n in enumerate(pref_n_sust)}
    df["_ord"] = df["n_sust"].apply(lambda n: orden.get(n, 99))
    df = df.sort_values(["_ord", "n_sust"]).reset_index(drop=True)

    seleccionados = set()
    for (gf1, _pos1), grp in df.groupby(["_gf1", "_pos1"], sort=False):
        if not gf1:
            continue
        n_tomar = min(max_por_gf, max(min_por_gf, len(grp) // 10 + 1))
        for idx in grp.index[:n_tomar]:
            seleccionados.add(idx)
        if len(seleccionados) >= max_por_mol:
            break

    if balance_ed_ew and len(seleccionados) >= 20:
        sel_tmp = df.loc[list(seleccionados)]
        min_cada = max(4, int(len(sel_tmp) * pct_min_ed_ew))
        n_ed = int((sel_tmp["_ed_ew"] == "ED").sum())
        n_ew = int((sel_tmp["_ed_ew"] == "EW").sum())
        if n_ed < min_cada:
            pool = df[(df["_ed_ew"] == "ED") & (~df.index.isin(seleccionados))]
            seleccionados.update(pool.index[:min_cada - n_ed])
        if n_ew < min_cada:
            pool = df[(df["_ed_ew"] == "EW") & (~df.index.isin(seleccionados))]
            seleccionados.update(pool.index[:min_cada - n_ew])

    sel = df.loc[sorted(seleccionados)].copy()
    if len(sel) > max_por_mol:
        sel = sel.sort_values(["_gf1", "_pos1", "_ord"]).head(max_por_mol)

    sel["ED_EW"] = sel["_gf1"].apply(clase_ed_ew)
    sel = sel.drop(columns=[c for c in sel.columns if c.startswith("_")])
    sel = sel.reset_index(drop=True)

    ed_counts = {"ED": 0, "EW": 0, "N": 0}
    for v in sel["ED_EW"]:
        ed_counts[v] = ed_counts.get(v, 0) + 1
    gf_counts, pos_counts = {}, {}
    for col in [c for c in sel.columns if c.startswith("sub_") and c.endswith("_gf")]:
        for v in sel[col]:
            if v: gf_counts[v] = gf_counts.get(v, 0) + 1
    for col in [c for c in sel.columns if c.startswith("sub_") and c.endswith("_pos")]:
        for v in sel[col]:
            if v: pos_counts[v] = pos_counts.get(v, 0) + 1

    stats = {
        "n_sel": len(sel), "ed": ed_counts["ED"], "ew": ed_counts["EW"],
        "neutro": ed_counts["N"],
        "gf_top": sorted(gf_counts.items(), key=lambda x: -x[1])[:8],
        "pos_top": sorted(pos_counts.items(), key=lambda x: -x[1])[:5],
    }
    log(f"  Seleccionados: {len(sel)} | ED={ed_counts['ED']} EW={ed_counts['EW']} "
        f"N={ed_counts['N']}")
    return sel, stats
