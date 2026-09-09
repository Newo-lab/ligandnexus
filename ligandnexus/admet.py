"""
admet.py — Perfil ADMET / de fármaco-similitud
==============================================
Descriptores fisicoquímicos y reglas empíricas de biodisponibilidad para
cualquier conjunto de moléculas del cribado. Todo se calcula **localmente con
RDKit**: no hace falta red ni servicios externos.

Es un tamiz de PRIORIZACIÓN, no una predicción farmacocinética: las reglas son
correlaciones estadísticas sobre fármacos orales conocidos, así que una molécula
que las incumple no queda descartada — queda señalada.

Reglas implementadas, con su fuente:
  · **Lipinski** «regla de 5» — Lipinski et al., Adv Drug Deliv Rev 46 (2001) 3.
    MW <= 500, clogP <= 5, donadores H <= 5, aceptores H <= 10. Se admite hasta
    UNA violación, como en el enunciado original.
  · **Veber** — Veber et al., J Med Chem 45 (2002) 2615.
    Enlaces rotables <= 10 y TPSA <= 140 A^2 (biodisponibilidad oral en rata).
  · **Egan** («huevo» de absorción) — Egan et al., J Med Chem 43 (2000) 3867.
    TPSA <= 131.6 y clogP <= 5.88.
  · **logBB de Clark** — Clark, J Pharm Sci 88 (1999) 815.
    logBB = -0.0148 * TPSA + 0.152 * clogP + 0.139.
    logBB > 0.3 atraviesa bien la barrera hematoencefálica; < -1.0 la atraviesa
    mal. Es el criterio de barrera usado en la tesis.
  · **QED** (fármaco-similitud cuantitativa) — Bickerton et al., Nat Chem 4
    (2012) 90. De 0 (mal) a 1 (bien).
  · **SAscore** (accesibilidad sintética) — Ertl y Schuffenhauer, J Cheminform 1
    (2009) 8. De 1 (fácil) a 10 (difícil); <= 3.5 se considera accesible.
  · **PAINS** (compuestos que interfieren en los ensayos) — Baell y Holloway,
    J Med Chem 53 (2010) 2719. Catálogo incluido en RDKit.

OJO con la tensión Lipinski/logBB: la ecuación de Clark PREMIA el clogP alto
mientras que Lipinski lo castiga, de modo que optimizar hacia el cerebro empuja
en contra de la regla de 5. Las dos columnas se leen juntas, nunca una sola.

DOS AVISOS MEDIDOS, para no leer estos números como más exactos de lo que son:

1. La **TPSA** de aquí es la de Ertl clásica (solo N y O), que es la definición
   con la que Clark ajustó su ecuación. La columna `TPSA` que trae PubChem puede
   no coincidir: incluye S y P en algunos compuestos, y además la percepción de
   aromaticidad difiere (para la cafeína, RDKit da 61.8 A^2 y PubChem 58.4).
   Por eso las columnas calculadas aquí llevan el sufijo `_admet` cuando chocan
   con una columna que ya venía de la base de datos: son magnitudes distintas y
   no deben pisarse.

2. El logBB de Clark es una **correlación de dos descriptores**, no un modelo
   fino. Comprobado contra fármacos conocidos, sitúa bien los extremos
   (ciclosporina -3.5, quercetina -1.5, levodopa -1.4: ninguno entra al cerebro)
   pero se queda corto con moléculas pequeñas que sí entran: la cafeína sale
   -0.93 y el donepezilo 0.23, ambos en la banda «intermedio» pese a ser
   fármacos del sistema nervioso central. Sirve para ORDENAR una biblioteca, no
   para dictaminar sobre una molécula concreta.
"""

from __future__ import annotations

import os
import sys

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, QED
from rdkit import RDLogger

from .chem import mol_from_smiles

RDLogger.DisableLog("rdApp.*")


# ── Cargas perezosas de los extras de RDKit ───────────────────────────────────

_sascorer = None
_sa_intentado = False


def _sa():
    """Módulo SAscore de los Contrib de RDKit (se carga una sola vez)."""
    global _sascorer, _sa_intentado
    if not _sa_intentado:
        _sa_intentado = True
        try:
            from rdkit.Chem import RDConfig
            ruta = os.path.join(RDConfig.RDContribDir, "SA_Score")
            if ruta not in sys.path:
                sys.path.append(ruta)
            import sascorer
            _sascorer = sascorer
        except Exception:
            _sascorer = None
    return _sascorer


_pains = None
_pains_intentado = False


def _catalogo_pains():
    """Catálogo PAINS de RDKit (se construye una sola vez: es caro)."""
    global _pains, _pains_intentado
    if not _pains_intentado:
        _pains_intentado = True
        try:
            from rdkit.Chem import FilterCatalog
            from rdkit.Chem.FilterCatalog import FilterCatalogParams
            p = FilterCatalogParams()
            p.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
            _pains = FilterCatalog.FilterCatalog(p)
        except Exception:
            _pains = None
    return _pains


# ── Umbrales, en un solo sitio para poder citarlos en la interfaz ─────────────

LIPINSKI = {"MW": 500.0, "clogP": 5.0, "HBD": 5, "HBA": 10}
VEBER = {"RotB": 10, "TPSA": 140.0}
EGAN = {"TPSA": 131.6, "clogP": 5.88}
SA_ACCESIBLE = 3.5
LOGBB_ALTO = 0.3        # por encima: atraviesa bien la barrera
LOGBB_BAJO = -1.0       # por debajo: la atraviesa mal


def logbb_clark(tpsa, clogp) -> float:
    """logBB por la ecuación de Clark (1999). Ver el encabezado del módulo."""
    return -0.0148 * float(tpsa) + 0.152 * float(clogp) + 0.139


def clase_bhe(logbb) -> str:
    """Lectura cualitativa del logBB (barrera hematoencefálica)."""
    if logbb is None:
        return ""
    if logbb > LOGBB_ALTO:
        return "atraviesa"
    if logbb < LOGBB_BAJO:
        return "no atraviesa"
    return "intermedio"


def descriptores(entrada) -> dict:
    """
    Descriptores fisicoquímicos de una molécula (SMILES o Mol de RDKit).
    Devuelve {} si la entrada no es interpretable.
    """
    mol = entrada if isinstance(entrada, Chem.Mol) else mol_from_smiles(entrada)
    if mol is None:
        return {}
    d = {
        "MW": round(Descriptors.MolWt(mol), 2),
        "clogP": round(Crippen.MolLogP(mol), 2),
        "TPSA": round(Descriptors.TPSA(mol), 1),
        "HBD": int(Descriptors.NumHDonors(mol)),
        "HBA": int(Descriptors.NumHAcceptors(mol)),
        "RotB": int(Descriptors.NumRotatableBonds(mol)),
        "Anillos": int(Descriptors.RingCount(mol)),
        "AtomosPesados": int(mol.GetNumHeavyAtoms()),
        "FraccionCsp3": round(Descriptors.FractionCSP3(mol), 3),
        "RefracMolar": round(Crippen.MolMR(mol), 2),
    }
    try:
        d["QED"] = round(QED.qed(mol), 3)
    except Exception:
        d["QED"] = None
    sa = _sa()
    try:
        d["SAscore"] = round(sa.calculateScore(mol), 2) if sa else None
    except Exception:
        d["SAscore"] = None
    return d


def evaluar(entrada) -> dict:
    """
    Perfil ADMET completo de una molécula: descriptores + reglas + logBB.
    Devuelve {} si la entrada no es interpretable.
    """
    mol = entrada if isinstance(entrada, Chem.Mol) else mol_from_smiles(entrada)
    if mol is None:
        return {}
    d = descriptores(mol)

    # Lipinski: se cuentan las violaciones y se admite UNA, como en el original.
    viol = []
    if d["MW"] > LIPINSKI["MW"]:
        viol.append("MW {} > {:.0f}".format(d["MW"], LIPINSKI["MW"]))
    if d["clogP"] > LIPINSKI["clogP"]:
        viol.append("clogP {} > {:.0f}".format(d["clogP"], LIPINSKI["clogP"]))
    if d["HBD"] > LIPINSKI["HBD"]:
        viol.append("HBD {} > {}".format(d["HBD"], LIPINSKI["HBD"]))
    if d["HBA"] > LIPINSKI["HBA"]:
        viol.append("HBA {} > {}".format(d["HBA"], LIPINSKI["HBA"]))
    d["Lipinski_viol"] = len(viol)
    d["Lipinski_detalle"] = "; ".join(viol)
    d["Lipinski"] = len(viol) <= 1

    d["Veber"] = d["RotB"] <= VEBER["RotB"] and d["TPSA"] <= VEBER["TPSA"]
    d["Egan"] = d["TPSA"] <= EGAN["TPSA"] and d["clogP"] <= EGAN["clogP"]

    lbb = logbb_clark(d["TPSA"], d["clogP"])
    d["logBB"] = round(lbb, 3)
    d["BHE"] = clase_bhe(lbb)

    d["SA_accesible"] = (d["SAscore"] is not None and d["SAscore"] <= SA_ACCESIBLE)

    cat = _catalogo_pains()
    if cat is not None:
        try:
            m = cat.GetFirstMatch(mol)
            d["PAINS"] = bool(m)
            d["PAINS_alerta"] = m.GetDescription() if m else ""
        except Exception:
            d["PAINS"], d["PAINS_alerta"] = None, ""
    else:
        d["PAINS"], d["PAINS_alerta"] = None, ""

    return d


# Columnas del informe, en el orden en que conviene leerlas.
COLUMNAS = ["MW", "clogP", "TPSA", "HBD", "HBA", "RotB", "FraccionCsp3",
            "QED", "SAscore", "logBB", "BHE",
            "Lipinski", "Lipinski_viol", "Veber", "Egan", "SA_accesible",
            "PAINS", "PAINS_alerta", "Lipinski_detalle"]


# ── Memoria de perfiles ya calculados ─────────────────────────────────────────
# El perfil depende SOLO del SMILES, nunca del conjunto en el que venga. Sin
# esta memoria, filtrar el cribado (mover un deslizador, endurecer un umbral)
# obligaba a repetir el cálculo entero: medido a 7,5 ms/molécula, son 5,1 min
# por cada retoque en una corrida de 40.000. Con ella, recalcular es gratis.
_MEMO: dict = {}
_MEMO_TOPE = 200_000        # más allá se vacía: es una caché, no un almacén.


def faltan_por_calcular(smiles) -> int:
    """Cuántos de estos SMILES NO están memorizados (0 = recalcular es gratis)."""
    return sum(1 for s in smiles if str(s) not in _MEMO)


def olvidar_memoria() -> None:
    """Vacía la memoria de perfiles."""
    _MEMO.clear()


def _nombres_sin_chocar(nuevas, ocupados):
    """
    Renombra las columnas de `nuevas` que choquen con `ocupados`, con el sufijo
    `_admet`. Si ESE nombre también está ocupado (pasa al perfilar dos veces el
    mismo conjunto: `TPSA` -> `TPSA_admet` -> `TPSA_admet` otra vez), numera.
    Devolver dos columnas con el mismo nombre rompía la exportación a Excel:
    `row.get()` entregaba una Serie y openpyxl la rechazaba.
    """
    ocupados = set(ocupados)
    mapa = {}
    for c in nuevas:
        if c not in ocupados:
            ocupados.add(c)
            continue
        n, i = c + "_admet", 2
        while n in ocupados:
            n, i = f"{c}_admet{i}", i + 1
        mapa[c] = n
        ocupados.add(n)
    return mapa


def evaluar_dataframe(df, smiles_col, prefijo="", progreso=None, memoria=True):
    """
    Devuelve una copia de `df` con las columnas ADMET añadidas, calculadas a
    partir de `smiles_col`. `progreso(i, total)` se llama cada 100 filas para
    poder pintar una barra. Las filas con SMILES ilegible quedan vacías.

    Si alguna columna calculada coincide de nombre con una que ya traía `df`
    (el caso típico es `TPSA`, que PubChem también entrega pero con otra
    definición), la NUEVA se renombra con el sufijo `_admet` en vez de duplicar
    el nombre: son magnitudes distintas y conviene poder compararlas.

    Con `memoria=True` los perfiles ya calculados se reutilizan entre llamadas
    (ver `_MEMO`), de modo que volver a perfilar un subconjunto es inmediato.
    """
    import pandas as pd

    total = len(df)
    filas = []
    for i, smi in enumerate(df[smiles_col].astype(str), 1):
        d = _MEMO.get(smi) if memoria else None
        if d is None:
            d = evaluar(smi) or {}
            if memoria:
                if len(_MEMO) >= _MEMO_TOPE:
                    _MEMO.clear()
                _MEMO[smi] = d
        filas.append(d)
        if progreso is not None and (i % 100 == 0 or i == total):
            progreso(i, total)
    add = pd.DataFrame(filas, index=df.index)
    add = add.reindex(columns=[c for c in COLUMNAS if c in add.columns])
    if prefijo:
        add = add.add_prefix(prefijo)
    add = add.rename(columns=_nombres_sin_chocar(add.columns, df.columns))
    return pd.concat([df.copy(), add], axis=1)


def resumen(df_admet, prefijo="") -> dict:
    """Conteos agregados, para las tarjetas de la interfaz."""
    def col(c):
        nombre = prefijo + c
        return df_admet[nombre] if nombre in df_admet.columns else None

    out = {"n": len(df_admet)}
    for regla in ("Lipinski", "Veber", "Egan", "SA_accesible"):
        c = col(regla)
        out[regla] = int(c.fillna(False).astype(bool).sum()) if c is not None else 0
    p = col("PAINS")
    out["PAINS"] = int(p.fillna(False).astype(bool).sum()) if p is not None else 0
    b = col("BHE")
    if b is not None:
        out["BHE_atraviesa"] = int((b == "atraviesa").sum())
        out["BHE_intermedio"] = int((b == "intermedio").sum())
        out["BHE_no"] = int((b == "no atraviesa").sum())
    return out
