"""
metales.py — Complejos de metales de transición con polifenoles
===============================================================
Coordinación de **Fe(II)/Fe(III) y Cu(I)/Cu(II)** por polifenoles. Reconoce los
cinco motivos quelantes O,O que cubren tanto los taninos hidrolizables como los
flavonoides: catecoles, galoílos (pirogalol, ácido gálico), salicilatos y las dos
**hidroxi-cetonas** de los flavonoides (3-OH/4-C=O y 5-OH/4-C=O).

Diseño (coherencia química):
  * El motor de donadores de `chem.py` es metal-agnóstico y trata cada -OH como un
    donador independiente. Para estos metales eso pierde lo esencial: **dos -OH orto
    (catecol) o el triplete del galoílo se desprotonan y forman UN quelato bidentado
    O,O** — el motivo catecolato es de los quelantes de Fe(III) más fuertes que se
    conocen (enterobactina), y el Cu(II) además lo oxida a o-quinona (química redox
    tipo Fenton, relevante en estrés oxidativo).
  * La carga del sitio la fija su química, no una constante: catecolato, galato y
    salicilato quelan como **dianiones** (dos -OH que se desprotonan), mientras que
    la hidroxi-cetona lo hace como **monoanión** (un fenolato + un C=O neutro que
    dona su par solitario sin perder H).
  * Dos modos de ensamblaje 3D como SEMILLA para DFT (no una optimización DFT):
      - **Mononuclear:** un metal con 1–3 ligandos (mono/bis/tris) en geometría
        octaédrica o plano-cuadrada según el metal.
      - **Polinuclear:** una sola molécula polidentada (p. ej. un tanino) con varios
        metales, uno por sitio de quelación — el caso real de los taninos.

Trabaja SIEMPRE sobre la molécula con la aromaticidad ya percibida por RDKit.

Validación externa (bloque [13] de `tests/verificar_motores.py`):
  Zaccaron, S.; Ganzerla, R.; Bortoluzzi, M. «Iron complexes with gallic acid: a
  computational study on coordination compounds of interest for the preservation
  of cultural heritage». J. Coord. Chem. 2013, 66, 1709–1719.
  DOI 10.1080/00958972.2013.790019.
  Contrastado contra ese trabajo y contra las distancias cristalográficas que
  cita (Wunderlich, Weber, Bergerhoff, Z. Anorg. Allg. Chem. 598/599 (1991) 371,
  el único complejo de tinta ferrogálica resuelto por rayos X):
    · el ácido gálico quela por los FENOLATOS (1 galoílo) y no por el carboxilato,
      que en el cristal PUENTEA dos hierros; y no hay salicilato porque sus -OH
      están en 3,4,5, es decir meta y para respecto al -COOH, no orto;
    · número de coordinación 6 («six Fe–O bonds»), completando con agua;
    · Fe–O 2,01 Å frente a 2,000 / 2,006 / 2,028 Å del cristal (media 2,011);
    · Fe(III) sexteto y Fe(II) quinteto, como sus referencias [Fe(EDTA)]⁻/²⁻.
  LIMITACIONES declaradas: este motor NO modela un oxígeno puenteando dos metales
  (su modo polinuclear pone un metal por sitio), y no predice densidades de espín
  — el hallazgo central del artículo, que esos Fe(III) formales se comportan como
  Fe(II) d⁶ de alto espín con ligandos radicalarios, solo sale de un cálculo DFT.
"""

from __future__ import annotations

from itertools import combinations

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger

from .chem import mol_from_smiles

RDLogger.DisableLog("rdApp.*")


# ── Motivos quelantes (tras percepción de aromaticidad) ───────────────────────
_CATECOL = Chem.MolFromSmarts(
    "[$([cX3]([OX2H,OX1-]))]:[$([cX3]([OX2H,OX1-]))]")
_GALOIL = Chem.MolFromSmarts(
    "[$(c[OX2H,OX1-])]:[$(c[OX2H,OX1-])]:[$(c[OX2H,OX1-])]")
_SALICIL = Chem.MolFromSmarts(
    "[OX2H,OX1-]-c:c-C(=O)[OX2H,OX1-]")

# Carbono de CETONA: lleva =O y su otro vecino pesado es carbono. Excluye
# explícitamente ácidos, ésteres y amidas (que ya cubre `_SALICIL` o que no
# quelan por esta vía).
_KETONA = ("[$([CX3,cX3](=[OX1])[#6]);"
           "!$([CX3,cX3](=[OX1])[OX2,OX1-]);"
           "!$([CX3,cX3](=[OX1])[NX3])]")
# o-hidroxi-cetona: el motivo quelante de los FLAVONOIDES, que el catecol y el
# salicilato no cubren. Dos variantes según el tamaño del anillo quelato:
#   · 5 miembros → 3-OH/4-C=O de los flavonoles (quercetina, kaempferol)
#   · 6 miembros → 5-OH/4-C=O de las flavonas y las 2'-hidroxichalconas
_HIDROXICETONA_5 = Chem.MolFromSmarts(
    f"[OX2H,OX1-]-[cX3,CX3]=,:{_KETONA}=[OX1]")
_HIDROXICETONA_6 = Chem.MolFromSmarts(
    f"[OX2H,OX1-]-[cX3,CX3]=,:[cX3,CX3]-,=,:{_KETONA}=[OX1]")

FUERTE = "fuerte"
MODERADO = "moderado"


# ── Parámetros de los metales soportados ──────────────────────────────────────
# Distancias M–O tomadas de promedios cristalográficos (CSD) para el motivo
# catecolato/O,O de metales de transición 3d:
#   q: carga del ion · d: electrones d · mult: multiplicidad de espín (alto espín) ·
#   mo: distancia M–O ecuatorial/principal (Å) · mo_ax: M–O axial (Jahn-Teller, Cu²⁺) ·
#   coord: politopo de coordinación (oct/tetragonal/td) · max_bi: nº máx. de
#   ligandos bidentados · geom: descripción legible.
# El ángulo de mordida O–M–O NO se tabula por metal: se DERIVA en `construir_complejo`
# del O···O cristalográfico que impone el tamaño del anillo quelato (2.61 Å para los
# de cinco miembros, 2.72 Å para los de seis) junto con la M–O del metal. Así un mismo
# metal abre la mordida al pasar de un catecolato a un salicilato, que es lo que se
# observa. El O···O no se lee del campo de fuerzas, que lo sobrestima por la repulsión
# O⁻/O⁻ del dianión libre; los O se fijan a esa geometría exacta y el UFF restringido
# reacomoda la tensión del anillo orgánico.
METALES = {
    "Fe3+": {"simbolo": "Fe", "q": 3, "d": 5, "mult": 6, "mo": 2.01,
             "coord": "oct", "max_bi": 3, "geom": "octaédrica (alto espín)",
             "etiqueta": "Fe(III)"},
    "Fe2+": {"simbolo": "Fe", "q": 2, "d": 6, "mult": 5, "mo": 2.14,
             "coord": "oct", "max_bi": 3, "geom": "octaédrica (alto espín)",
             "etiqueta": "Fe(II)"},
    "Cu2+": {"simbolo": "Cu", "q": 2, "d": 9, "mult": 2, "mo": 1.95,
             "mo_ax": 2.35, "axial_agua": True, "coord": "tetragonal", "max_bi": 2,
             "geom": "octaédrica tetragonal (Jahn-Teller)", "etiqueta": "Cu(II)"},
    "Cu+":  {"simbolo": "Cu", "q": 1, "d": 10, "mult": 1, "mo": 2.10,
             "coord": "td", "max_bi": 2, "geom": "tetraédrica", "etiqueta": "Cu(I)"},
}


def metal_info(metal):
    return METALES.get(metal, METALES["Fe3+"])


# ── Detección de sitios de quelación ──────────────────────────────────────────

def _es_carbonilo(o_atom):
    return any(b.GetBondType() == Chem.BondType.DOUBLE for b in o_atom.GetBonds())


def _o_de_carbono(mol, c_idx):
    """Índice del O (hidroxilo/fenolato) unido al carbono aromático `c_idx`."""
    for nb in mol.GetAtomWithIdx(c_idx).GetNeighbors():
        if nb.GetAtomicNum() == 8 and nb.GetDegree() <= 2 and not _es_carbonilo(nb):
            return nb.GetIdx()
    return None


def _desprotonar_sitio(mw, o_indices):
    """
    Desprotona los O de un sitio quelante **solo donde corresponde**: quita el H de
    los hidroxilos (fenol, carboxilo) y deja INTACTO el O carbonílico, que dona su
    par solitario sin perder H.

    Ponerle carga −1 a un O con doble enlace da una valencia imposible y hace
    fracasar `SanitizeMol`, así que la carga del sitio la fija la química, no una
    constante: catecolato/galato/salicilato son **dianiónicos** (dos −OH) y la
    hidroxi-cetona es **monoaniónica** (un −OH + un C=O neutro).
    Devuelve False si algún índice no es oxígeno.
    """
    for o in o_indices:
        a = mw.GetAtomWithIdx(int(o))
        if a.GetAtomicNum() != 8:
            return False
        if a.GetFormalCharge() < 0 or a.GetTotalNumHs() == 0:
            continue                      # ya desprotonado, o carbonílico: intacto
        a.SetFormalCharge(-1)
        a.SetNoImplicit(True)
        a.SetNumExplicitHs(0)
    return True


def _add_hidroxicetona(mol, patron, fuerza, etiqueta, sitios, o_usados):
    """Registra los sitios o-hidroxi-cetona que casen con `patron` sin reutilizar
    oxígenos ya asignados a otro sitio. Devuelve cuántos añadió."""
    n = 0
    for match in mol.GetSubstructMatches(patron):
        o_fenol, o_ceto = match[0], match[-1]
        if o_fenol in o_usados or o_ceto in o_usados:
            continue
        sitios.append({"tipo": etiqueta, "fuerza": fuerza,
                       "o_idx": (o_fenol, o_ceto), "c_idx": tuple(match[1:-1])})
        o_usados |= {o_fenol, o_ceto}
        n += 1
    return n


def analizar_polifenol(smiles) -> dict:
    """
    Detecta los sitios de quelación de un polifenol.

    Devuelve:
      {"mol", "sitios":[{tipo, fuerza, o_idx:(o1,o2), c_idx:(...)}, ...],
       "n_galoil", "n_catecol", "n_salicil", "n_hidroxicetona"}

    Motivos reconocidos, en orden de afinidad decreciente (galoílo · catecol ·
    hidroxi-cetona de 5 · salicilato · hidroxi-cetona de 6). El orden importa
    porque los sitios NO comparten oxígenos: el primero que reclama un O se lo
    queda. Así, en un flavonol el sitio 3-OH/4-C=O (el preferido) toma el
    carbonilo antes de que lo intente el 5-OH/4-C=O, que es lo que ocurre en la
    práctica. Los catecoles internos de un galoílo tampoco se cuentan aparte (el
    galoílo es un único sitio bidentado).
    """
    mol = mol_from_smiles(smiles)
    if mol is None:
        return {"mol": None, "sitios": [], "n_galoil": 0, "n_catecol": 0, "n_salicil": 0}

    sitios, o_usados, carbonos_galoil = [], set(), set()

    # 1) Galoílos (1,2,3-trihidroxi) → un sitio bidentado fuerte cada uno.
    vistos_gal, n_galoil = set(), 0
    for c1, c2, c3 in mol.GetSubstructMatches(_GALOIL):
        clave = frozenset((c1, c2, c3))
        if clave in vistos_gal:
            continue
        vistos_gal.add(clave)
        carbonos_galoil |= {c1, c2, c3}
        o1, o2 = _o_de_carbono(mol, c1), _o_de_carbono(mol, c2)
        if o1 is None or o2 is None or o1 in o_usados or o2 in o_usados:
            continue
        sitios.append({"tipo": "galoíl", "fuerza": FUERTE,
                       "o_idx": (o1, o2), "c_idx": (c1, c2, c3)})
        o_usados |= {o1, o2}
        n_galoil += 1

    # 2) Catecoles (o-difenol) que NO forman parte de un galoílo.
    vistos_cat, n_catecol = set(), 0
    for c1, c2 in mol.GetSubstructMatches(_CATECOL):
        clave = frozenset((c1, c2))
        if clave in vistos_cat:
            continue
        vistos_cat.add(clave)
        if c1 in carbonos_galoil and c2 in carbonos_galoil:
            continue
        o1, o2 = _o_de_carbono(mol, c1), _o_de_carbono(mol, c2)
        if o1 is None or o2 is None or o1 in o_usados or o2 in o_usados:
            continue
        sitios.append({"tipo": "catecol", "fuerza": FUERTE,
                       "o_idx": (o1, o2), "c_idx": (c1, c2)})
        o_usados |= {o1, o2}
        n_catecol += 1

    # 3) Hidroxi-cetona de anillo 5 (3-OH/4-C=O de los flavonoles): el fenolato y
    #    el carbonilo cierran un quelato de cinco miembros, tan fuerte como el
    #    catecolato. Va ANTES que la variante de seis porque ambas comparten el
    #    mismo O carbonílico y solo una puede ocuparlo.
    n_hc5 = _add_hidroxicetona(mol, _HIDROXICETONA_5, FUERTE,
                               "hidroxi-cetona-5", sitios, o_usados)

    # 4) Salicilatos (o-hidroxi-carboxilo) → sitio bidentado moderado.
    n_salicil = 0
    for match in mol.GetSubstructMatches(_SALICIL):
        o_fenol, o_carbox = match[0], match[-1]
        if o_fenol in o_usados or o_carbox in o_usados:
            continue
        c_ar = tuple(a for a in match if mol.GetAtomWithIdx(a).GetIsAromatic())
        sitios.append({"tipo": "salicilato", "fuerza": MODERADO,
                       "o_idx": (o_fenol, o_carbox), "c_idx": c_ar})
        o_usados |= {o_fenol, o_carbox}
        n_salicil += 1

    # 5) Hidroxi-cetona de anillo 6 (5-OH/4-C=O de flavonas y 2'-hidroxichalconas).
    #    Moderado: el 5-OH está trabado en un puente de hidrógeno intramolecular
    #    fuerte, así que compite con el metal por ese protón.
    n_hc6 = _add_hidroxicetona(mol, _HIDROXICETONA_6, MODERADO,
                               "hidroxi-cetona-6", sitios, o_usados)

    return {"mol": mol, "sitios": sitios, "n_galoil": n_galoil,
            "n_catecol": n_catecol, "n_salicil": n_salicil,
            "n_hidroxicetona": n_hc5 + n_hc6}


# ── Especiación metal–polifenol vs pH ─────────────────────────────────────────

def especiacion(metal, n_sitios_por_mol, pH=7.0) -> dict:
    """
    Estima la especie dominante (mono/bis/tris) y el color en función del pH para
    un ligando con `n_sitios_por_mol` sitios bidentados.
    """
    inf = metal_info(metal)
    tope = inf["max_bi"]

    if metal == "Fe3+":
        n_lig, especie, color = ((1, "mono", "verde-azulado") if pH < 3 else
                                 (2, "bis", "azul-violeta") if pH < 6.5 else
                                 (3, "tris", "rojo"))
        nota = ("Fe(III) es un ácido duro: afinidad enorme por los O aniónicos del "
                "catecolato/galato. El tris-catecolato de Fe(III) es de los quelatos "
                "más estables conocidos (log β ≳ 40). Es el viraje del hierro-galato "
                "(tinta ferrogálica): verde→azul→rojo al subir el pH.")
    elif metal == "Fe2+":
        n_lig, especie, color = ((1, "mono", "débil / incolora") if pH < 5 else
                                 (2, "bis", "pardo") if pH < 8 else
                                 (3, "tris", "pardo-rojizo"))
        nota = ("Fe(II) coordina más débil que Fe(III) y, con polifenoles y O₂, "
                "tiende a OXIDARSE a Fe(III) (el polifenol pasa a quinona). Complejos "
                "de alto espín y sensibles al aire.")
    elif metal == "Cu2+":
        n_lig, especie, color = ((1, "mono", "verde") if pH < 4 else
                                 (2, "bis", "pardo-verdoso"))
        nota = ("Cu(II) (d⁹, plano-cuadrado, Jahn-Teller) forma catecolatos estables, "
                "pero además OXIDA el catecol a o-quinona con reducción a Cu(I) "
                "(química redox tipo Fenton, generación de ROS) — relevante en el "
                "estrés oxidativo y la química del cobre en Alzheimer.")
    else:  # Cu+
        n_lig, especie, color = ((1, "mono", "incolora / inestable") if pH < 6 else
                                 (2, "bis", "parda"))
        nota = ("Cu(I) es un ácido BLANDO: prefiere donadores N/S; con los O duros del "
                "catecolato su afinidad es baja y al aire se oxida a Cu(II). Se incluye "
                "por completitud, pero no es el modo dominante con polifenoles.")

    n_lig = min(n_lig, tope)
    n_lig_efectivo = min(n_lig, 3) if n_sitios_por_mol >= 1 else 0
    return {"n_ligandos": n_lig_efectivo, "especie": especie, "color": color,
            "nota": nota, "metal": metal, "pH": pH, "geom": inf["geom"]}



def carga_multiplicidad(metal, n_ligandos, carga_ligando=-2, mult_manual=None):
    """
    Carga total y multiplicidad del complejo mononuclear.
      * catecolato/galato desprotonado en el sitio = −2.
      * mult por defecto: Fe(III) d⁵ AE → 6; Fe(II) d⁶ AE → 5; Cu(II) d⁹ → 2; Cu(I) d¹⁰ → 1.
    """
    inf = metal_info(metal)
    carga = inf["q"] + n_ligandos * carga_ligando
    mult = int(mult_manual) if mult_manual is not None else inf["mult"]
    return int(carga), int(mult)


# ── Imagen 2D con los sitios de quelación resaltados ──────────────────────────

_COLOR_SITIO = [
    (0.20, 0.60, 0.86), (0.16, 0.65, 0.42), (0.90, 0.49, 0.13),
    (0.75, 0.35, 0.75), (0.85, 0.20, 0.30), (0.20, 0.70, 0.70),
    (0.95, 0.75, 0.15), (0.45, 0.45, 0.85),
]


def imagen_sitios(smiles, size=(520, 400), dark=False):
    """Imagen 2D con cada sitio de quelación resaltado y numerado."""
    from io import BytesIO
    from PIL import Image
    from rdkit.Chem.Draw import rdMolDraw2D

    info = analizar_polifenol(smiles)
    mol = info["mol"]
    if mol is None:
        return None, []
    sitios = info["sitios"]

    highlight, colores = [], {}
    for n, s in enumerate(sitios):
        col = _COLOR_SITIO[n % len(_COLOR_SITIO)]
        for a in list(s["o_idx"]) + list(s["c_idx"]):
            highlight.append(a)
            colores[a] = col
        mol.GetAtomWithIdx(s["o_idx"][0]).SetProp("atomNote", str(n + 1))

    d2d = rdMolDraw2D.MolDraw2DCairo(*size)
    if dark:
        rdMolDraw2D.SetDarkMode(d2d)
    opts = d2d.drawOptions()
    opts.setBackgroundColour((0, 0, 0, 0))
    opts.highlightBondWidthMultiplier = 12
    opts.annotationFontScale = 0.9
    rdMolDraw2D.PrepareAndDrawMolecule(
        d2d, mol, highlightAtoms=highlight, highlightAtomColors=colores)
    d2d.FinishDrawing()
    return Image.open(BytesIO(d2d.GetDrawingText())), sitios


# ── Geometría: utilidades comunes ─────────────────────────────────────────────

def _kabsch(P, Q):
    H = P.T @ Q
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    return Vt.T @ np.diag([1, 1, d]) @ U.T


def _rot_eje(axis, theta):
    axis = axis / np.linalg.norm(axis)
    x, y, z = axis
    c, s = np.cos(theta), np.sin(theta)
    C = 1 - c
    return np.array([[c + x*x*C, x*y*C - z*s, x*z*C + y*s],
                     [y*x*C + z*s, c + y*y*C, y*z*C - x*s],
                     [z*x*C - y*s, z*y*C + x*s, c + z*z*C]])


def _compress(u1, u2, half):
    """
    Dadas dos direcciones de vértice `u1,u2` de un politopo (bisectriz común), las
    reemplaza por dos direcciones simétricas separadas un semiángulo `half` (rad)
    respecto de la bisectriz. Sirve para fijar el ángulo de mordida O–M–O al valor
    que impone el O···O rígido del ligando, sin mover la bisectriz del slot.
    """
    u1 = np.asarray(u1, float); u2 = np.asarray(u2, float)
    b = u1 + u2
    nb = np.linalg.norm(b)
    b = b / nb if nb > 1e-6 else u1 / (np.linalg.norm(u1) or 1.0)
    p = u1 - (u1 @ b) * b
    npp = np.linalg.norm(p)
    if npp < 1e-9:                       # u1 ∥ b: elige un perpendicular cualquiera
        p = np.cross(b, np.array([0., 0., 1.]))
        npp = np.linalg.norm(p)
        if npp < 1e-9:
            p = np.cross(b, np.array([0., 1., 0.])); npp = np.linalg.norm(p)
    p = p / npp
    c, s = np.cos(half), np.sin(half)
    return c * b + s * p, c * b - s * p


# ── Politopos de coordinación ideales (slots = aristas bidentadas) ────────────
# Cada "slot" es un par de direcciones de vértice donde puede ir un ligando
# bidentado; los slots sobrantes se completan con dos aguas monodentadas.

def _slots_octaedro():
    """3 aristas cis de un octaedro en disposición D3 (tris-quelato)."""
    thm = np.radians(54.735610)

    def v(az_deg, polar):
        a = np.radians(az_deg)
        return np.array([np.sin(polar) * np.cos(a),
                         np.sin(polar) * np.sin(a), np.cos(polar)])
    return [(v(120 * k, thm), v(60 + 120 * k, np.pi - thm)) for k in range(3)]


def _slots_tetragonal():
    """2 slots ecuatoriales (plano xy) para Cu(II) tetragonal (bisectrices ±x)."""
    def v(a):
        return np.array([np.cos(np.radians(a)), np.sin(np.radians(a)), 0.])
    return [(v(45), v(-45)), (v(135), v(225))]


def _slots_tetraedro():
    """2 slots ortogonales (D2d) de un tetraedro para Cu(I)."""
    t = 1.0 / np.sqrt(3.0)
    a = np.array([t, t, t]); b = np.array([t, -t, -t])
    c = np.array([-t, t, -t]); d = np.array([-t, -t, t])
    return [(a, b), (c, d)]


def _geometria_coord(metal, n_lig):
    """
    Reparte los slots del politopo del metal entre `n_lig` ligandos bidentados;
    los slots restantes se completan con aguas para saturar la coordinación.
    Devuelve (slots_ligando, slots_agua, dirs_axiales, inf).
    """
    inf = metal_info(metal)
    coord = inf.get("coord", "oct")
    slots = (_slots_octaedro() if coord == "oct"
             else _slots_tetragonal() if coord == "tetragonal"
             else _slots_tetraedro())
    n_lig = max(1, min(int(n_lig), len(slots)))
    axial = [np.array([0., 0., 1.]), np.array([0., 0., -1.])] if inf.get("axial_agua") else []
    return slots[:n_lig], slots[n_lig:], axial, inf


def _clashes(mol, ignora=None):
    conf = mol.GetConformer()
    pos = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
    ignora = ignora or set()
    n = 0
    for i, j in combinations(range(mol.GetNumAtoms()), 2):
        if mol.GetBondBetweenAtoms(i, j) or i in ignora or j in ignora:
            continue
        d = np.linalg.norm(pos[i] - pos[j])
        hi = mol.GetAtomWithIdx(i).GetAtomicNum() == 1
        hj = mol.GetAtomWithIdx(j).GetAtomicNum() == 1
        thr = 1.0 if (hi and hj) else (1.2 if (hi or hj) else 1.6)
        if d < thr:
            n += 1
    return n


def _posicion_metal(pos, o1, o2, c_ar, mo):
    """
    Coloca el metal equidistante (= mo) de los dos O donores, en el círculo de la
    bisectriz perpendicular a O1–O2, hacia afuera del anillo. Garantiza M–O = mo.
    """
    O1, O2 = pos[o1], pos[o2]
    mid = (O1 + O2) / 2.0
    eu = O2 - O1
    nu = np.linalg.norm(eu)
    eu = eu / (nu or 1.0)
    h = nu / 2.0
    r = np.sqrt(max(mo * mo - h * h, 0.09))
    rc = pos[list(c_ar)].mean(0) if len(c_ar) else mid
    outward = mid - rc
    outward = outward - (outward @ eu) * eu          # quitar la componente a lo largo de O1–O2
    n_out = np.linalg.norm(outward)
    ew = outward / n_out if n_out > 1e-6 else np.cross(eu, [0, 0, 1.0])
    return mid + r * ew


# ── Aguas de saturación de la esfera de coordinación ──────────────────────────

def _add_agua(combo, xyz, o_pos, m_pos, m_idx, frozen):
    """Añade una molécula de agua con el O en `o_pos`, dativa al metal `m_idx`.
    Los H se orientan alejándose del metal (par solitario apuntando al centro)."""
    o_pos = np.asarray(o_pos, float); m_pos = np.asarray(m_pos, float)
    o = combo.AddAtom(Chem.Atom(8)); xyz.append(o_pos.copy())
    d = o_pos - m_pos
    d = d / (np.linalg.norm(d) or 1.0)
    perp = np.cross(d, np.array([0., 0., 1.]))
    if np.linalg.norm(perp) < 1e-3:
        perp = np.cross(d, np.array([0., 1., 0.]))
    perp = perp / (np.linalg.norm(perp) or 1.0)
    ang = np.radians(52.25)                          # ½ del ángulo H–O–H (104.5°)
    for s in (1.0, -1.0):
        hd = np.cos(ang) * d + s * np.sin(ang) * perp
        h = combo.AddAtom(Chem.Atom(1)); xyz.append(o_pos + 0.96 * hd)
        combo.AddBond(o, h, Chem.BondType.SINGLE)
    combo.AddBond(o, m_idx, Chem.BondType.DATIVE)
    frozen.add(o)
    return o


# ── Relajación con la esfera de coordinación congelada (UFF restringido) ──────

def _relajar_restringido(mol, frozen):
    """
    Minimización UFF que relaja SOLO el ligando orgánico manteniendo **congelada**
    la esfera de coordinación (metal + O donores/agua + C ipso del quelato).

    UFF NO tiene tipo de átomo para los metales aislados (p. ej. «Fe2+2»); si el
    metal entra en el cálculo, sus parámetros basura contaminan la minimización y
    llegan a arrastrar un H sobre el O. Por eso el cálculo se hace sobre una copia
    **sin los metales** (átomos eliminados, no solo sus enlaces dativos) y sin los
    enlaces dativos, sanitizada: así los ligandos (fenolatos) y las aguas quedan
    perfectamente parametrizados y solo se relaja lo orgánico, con los átomos de la
    esfera de coordinación fijos. El metal no se mueve (se copia tal cual).
    Devuelve (mol, ok) con las coordenadas relajadas en el `mol` original.
    """
    try:
        rw = Chem.RWMol(mol)
        for b in list(rw.GetBonds()):
            if b.GetBondType() == Chem.BondType.DATIVE:
                rw.RemoveBond(b.GetBeginAtomIdx(), b.GetEndAtomIdx())
        metal_idx = sorted((a.GetIdx() for a in rw.GetAtoms()
                            if a.GetSymbol() in ("Fe", "Cu")), reverse=True)
        for mi in metal_idx:
            rw.RemoveAtom(mi)
        removed = set(metal_idx)
        old2new = {}
        for old in range(mol.GetNumAtoms()):
            if old not in removed:
                old2new[old] = old - sum(1 for r in removed if r < old)
        om = rw.GetMol()
        Chem.SanitizeMol(om)
        if not AllChem.UFFHasAllMoleculeParams(om):
            return mol, False
        ff = AllChem.UFFGetMoleculeForceField(om, ignoreInterfragInteractions=False)
        if ff is None:
            return mol, False
        for old in frozen:
            if old in old2new:
                ff.AddFixedPoint(old2new[old])
        ff.Initialize()
        ff.Minimize(maxIts=2000, energyTol=1e-6, forceTol=1e-4)
        src, dst = om.GetConformer(), mol.GetConformer()
        for old, new in old2new.items():
            dst.SetAtomPosition(old, src.GetAtomPosition(new))
        return mol, True
    except Exception:
        return mol, False


# ── Ensamblaje MONONUCLEAR (un metal, 1–3 ligandos) ───────────────────────────

def _prep_ligando(smiles, o_site):
    """Ligando desprotonado en los dos O del sitio, embebido en 3D."""
    m = mol_from_smiles(smiles)
    if m is None:
        return None, None
    o1, o2 = o_site
    mw = Chem.RWMol(m)
    if not _desprotonar_sitio(mw, (o1, o2)):
        return None, None
    m2 = mw.GetMol()
    try:
        Chem.SanitizeMol(m2)
    except Exception:
        return None, None
    m2 = Chem.AddHs(m2)
    p = AllChem.ETKDGv3()
    p.randomSeed = 42
    if AllChem.EmbedMolecule(m2, p) != 0:
        p.useRandomCoords = True
        if AllChem.EmbedMolecule(m2, p) != 0:
            return None, None
    AllChem.MMFFOptimizeMolecule(m2, maxIters=1000)
    return m2, (o1, o2)


def _twist(mol, ligand_atoms, pivotes):
    """Gira cada ligando sobre su eje O–O para minimizar choques entre ligandos."""
    conf = mol.GetConformer()
    pos = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
    for mine, (oa, ob) in zip(ligand_atoms, pivotes):
        axis = pos[ob] - pos[oa]
        origin = pos[oa]
        movibles = [a for a in mine if a not in (oa, ob)]
        otros = [a for a in range(mol.GetNumAtoms()) if a not in mine]
        best_t, best_pen = 0.0, 1e9
        for deg in range(-45, 46, 3):
            R = _rot_eje(axis, np.radians(deg))
            newp = {a: origin + R @ (pos[a] - origin) for a in movibles}
            pen = sum((1.9 - np.linalg.norm(newp[a] - pos[o])) ** 2
                      for a in movibles for o in otros
                      if np.linalg.norm(newp[a] - pos[o]) < 1.9)
            if pen < best_pen:
                best_pen, best_t = pen, np.radians(deg)
        R = _rot_eje(axis, best_t)
        for a in movibles:
            pos[a] = origin + R @ (pos[a] - origin)
    for i in range(mol.GetNumAtoms()):
        conf.SetAtomPosition(i, pos[i].tolist())
    return mol


def construir_complejo(smiles, o_site, n_ligandos, metal="Fe3+"):
    """
    Geometría 3D de un complejo MONONUCLEAR metal–polifenol (mono/bis/tris) sobre el
    politopo de coordinación rígido del metal, con la esfera de coordinación exacta
    (M–O cristalográfico, mordida derivada del ligando) y el orgánico relajado por
    UFF restringido. La coordinación se satura con aguas.
    Devuelve (mol, info) o (None, motivo).
    """
    inf = metal_info(metal)
    n_ligandos = max(1, min(inf["max_bi"], 3, int(n_ligandos)))
    lig_slots, water_slots, axial, _ = _geometria_coord(metal, n_ligandos)
    mo = inf["mo"]
    mo_ax = inf.get("mo_ax", mo)

    combo = Chem.RWMol()
    xyz = []
    fe = combo.AddAtom(Chem.Atom(inf["simbolo"]))
    xyz.append(np.zeros(3))
    odon, ligand_atoms, pivotes, frozen, bites = [], [], [], {fe}, []

    for (s1, s2) in lig_slots:
        lig, ox = _prep_ligando(smiles, o_site)
        if lig is None:
            return None, "No se pudo preparar el ligando (embedding fallido)."
        conf = lig.GetConformer()
        pos = np.array([list(conf.GetAtomPosition(i)) for i in range(lig.GetNumAtoms())])
        oa, ob = ox
        ring = [a.GetIdx() for a in lig.GetAtoms() if a.GetIsAromatic()] or [oa, ob]
        rc = pos[ring].mean(0)
        m = (pos[oa] + pos[ob]) / 2
        d_ring = np.linalg.norm(rc - m) or 1.0
        # Mordida O–M–O IMPUESTA: se fija el O···O cristalográfico según el TAMAÑO del
        # anillo quelato (5 miembros catecolato/galato ≈2.61 Å → ~81° con Fe³⁺;
        # 6 miembros salicilato ≈2.72 Å → mordida más abierta), y M–O = mo. La mordida
        # sale de la geometría, no del campo de fuerzas (que infla el O···O del dianión).
        anillo_quelato = len(Chem.GetShortestPath(lig, oa, ob)) + 1   # + el metal
        ooo = 2.72 if anillo_quelato >= 6 else 2.61
        half = float(np.arcsin(min(ooo / (2.0 * mo), 1.0)))
        bites.append(round(float(np.degrees(2 * half)), 1))
        n1, n2 = _compress(s1, s2, half)
        ta, tb = mo * n1, mo * n2
        tm = (ta + tb) / 2
        bis = tm / (np.linalg.norm(tm) or 1.0)
        P = np.array([pos[oa], pos[ob], rc])
        Q = np.array([ta, tb, tm + d_ring * bis])
        R = _kabsch(P - P.mean(0), Q - Q.mean(0))
        pos_new = (R @ (pos - P.mean(0)).T).T + Q.mean(0)
        pos_new[oa], pos_new[ob] = ta, tb        # O exactamente en el politopo rígido

        idxmap, mine = {}, []
        for i, at in enumerate(lig.GetAtoms()):
            na = Chem.Atom(at.GetAtomicNum())
            na.SetFormalCharge(at.GetFormalCharge())
            j = combo.AddAtom(na)
            idxmap[i] = j
            xyz.append(pos_new[i])
            mine.append(j)
        for b in lig.GetBonds():
            combo.AddBond(idxmap[b.GetBeginAtomIdx()], idxmap[b.GetEndAtomIdx()],
                          b.GetBondType())
        for o in ox:
            combo.AddBond(idxmap[o], fe, Chem.BondType.DATIVE)
            odon.append(idxmap[o])
            frozen.add(idxmap[o])          # solo el O donor: el UFF reacomoda el C–O
        ligand_atoms.append(mine)
        pivotes.append((idxmap[oa], idxmap[ob]))

    # Saturar la coordinación con aguas: cada slot bidentado vacío aporta DOS aguas
    # monodentadas (una por vértice) + las axiales del Jahn-Teller (Cu²⁺).
    # Se llevan por separado: las ecuatoriales están a la M-O normal y NO tienen
    # nada que ver con el Jahn-Teller; solo las axiales de Cu(II) están alargadas.
    # Mezclarlas hacía que la interfaz atribuyera una distorsión Jahn-Teller a
    # Fe(III) d5 y a Cu(I) d10, que por capa llena/semillena no la tienen.
    aguas_eq, aguas_ax = [], []
    for (u1, u2) in water_slots:
        for u in (u1, u2):
            aguas_eq.append(_add_agua(combo, xyz, mo * np.asarray(u),
                                      np.zeros(3), fe, frozen))
    for u in axial:
        aguas_ax.append(_add_agua(combo, xyz, mo_ax * np.asarray(u),
                                  np.zeros(3), fe, frozen))
    aguas = aguas_eq + aguas_ax

    mol = combo.GetMol()
    mol.GetAtomWithIdx(fe).SetFormalCharge(inf["q"])
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, c in enumerate(xyz):
        conf.SetAtomPosition(i, [float(c[0]), float(c[1]), float(c[2])])
    mol.AddConformer(conf, assignId=True)
    mol.UpdatePropertyCache(strict=False)
    Chem.GetSSSR(mol)
    _twist(mol, ligand_atoms, pivotes)
    mol, relajado = _relajar_restringido(mol, frozen)

    conf = mol.GetConformer()
    pos = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
    d = lambda o: round(float(np.linalg.norm(pos[o] - pos[fe])), 2)
    fe_o = [d(o) for o in odon]
    info = {"metal_idx": [fe], "odon": odon, "carga": int(Chem.GetFormalCharge(mol)),
            "multiplicidad": inf["mult"], "fe_o": fe_o, "choques": _clashes(mol),
            "n_ligandos": n_ligandos, "metal": metal, "n_metales": 1,
            "modo": "mononuclear", "bite": bites,
            "mo_agua_eq": [d(o) for o in aguas_eq],
            "mo_agua_ax": [d(o) for o in aguas_ax],
            "jahn_teller": bool(inf.get("axial_agua")),
            "n_aguas": len(aguas), "relajado": relajado}
    return mol, info


# ── Ensamblaje POLINUCLEAR (una molécula, varios metales) ─────────────────────

def construir_polinuclear(smiles, sitios_idx, metal="Fe3+"):
    """
    Geometría 3D de un complejo POLINUCLEAR: UNA molécula polidentada (p. ej. un
    tanino) con un metal en cada uno de los sitios de `sitios_idx`. Es el modelo
    realista para taninos, que enlazan varios metales a la vez. Tras colocar los
    metales (M–O = mo exacto) se relaja con UFF restringido dejando congelados los
    metales y los O donores.

    LIMITACIÓN deliberada: a diferencia del modo mononuclear, aquí los metales NO se
    saturan con agua, porque los sitios de un tanino quedan en orientaciones
    arbitrarias y unas aguas colocadas a ciegas chocarían con la propia molécula.
    Cada centro sale por tanto sólo bidentado (2-coordinado), muy insaturado: antes
    de la optimización DFT conviene completar su esfera de coordinación con agua o
    disolvente explícito.
    Devuelve (mol, info) o (None, motivo).
    """
    inf = metal_info(metal)
    mo = inf["mo"]
    analisis = analizar_polifenol(smiles)
    sitios = analisis["sitios"]
    # Un mismo sitio no puede alojar dos metales: se descartan repetidos (y se
    # conserva el orden de selección) además de los índices fuera de rango.
    vistos = set()
    sitios_idx = [i for i in (sitios_idx or [])
                  if isinstance(i, int) and 0 <= i < len(sitios)
                  and not (i in vistos or vistos.add(i))]
    if not sitios_idx:
        return None, "No hay sitios válidos seleccionados."

    base = mol_from_smiles(smiles)
    if base is None:
        return None, "SMILES no válido."
    mw = Chem.RWMol(base)
    triples = []
    for si in sitios_idx:
        o1, o2 = sitios[si]["o_idx"]
        _desprotonar_sitio(mw, (o1, o2))
        triples.append((o1, o2, tuple(sitios[si]["c_idx"])))
    m2 = mw.GetMol()
    try:
        Chem.SanitizeMol(m2)
    except Exception as e:
        return None, f"No se pudo desprotonar la molécula: {e}"
    m2 = Chem.AddHs(m2)
    p = AllChem.ETKDGv3()
    p.randomSeed = 42
    if AllChem.EmbedMolecule(m2, p) != 0:
        p.useRandomCoords = True
        if AllChem.EmbedMolecule(m2, p) != 0:
            return None, "No se pudo generar la geometría 3D de la molécula."
    AllChem.MMFFOptimizeMolecule(m2, maxIters=800)

    combo = Chem.RWMol(m2)
    conf = combo.GetConformer()
    pos = np.array([list(conf.GetAtomPosition(i)) for i in range(m2.GetNumAtoms())])
    metal_idxs, odon, frozen = [], [], set()
    for (o1, o2, c_ar) in triples:
        Mpos = _posicion_metal(pos, o1, o2, c_ar, mo)
        mi = combo.AddAtom(Chem.Atom(inf["simbolo"]))
        combo.GetAtomWithIdx(mi).SetFormalCharge(inf["q"])
        conf.SetAtomPosition(mi, [float(Mpos[0]), float(Mpos[1]), float(Mpos[2])])
        combo.AddBond(o1, mi, Chem.BondType.DATIVE)
        combo.AddBond(o2, mi, Chem.BondType.DATIVE)
        metal_idxs.append(mi)
        odon += [o1, o2]
        frozen |= {mi, o1, o2}

    mol = combo.GetMol()
    mol.UpdatePropertyCache(strict=False)
    Chem.GetSSSR(mol)
    mol, relajado = _relajar_restringido(mol, frozen)

    conf = mol.GetConformer()
    pos = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
    fe_o = [round(float(np.linalg.norm(pos[o] - pos[m])), 2)
            for m, (o1, o2, _) in zip(metal_idxs, triples) for o in (o1, o2)]
    mm = [round(float(np.linalg.norm(pos[a] - pos[b])), 2)
          for a, b in combinations(metal_idxs, 2)]
    n_m = len(metal_idxs)
    carga = int(Chem.GetFormalCharge(mol))
    # Multiplicidad: acoplamiento no interactuante de n centros de alto espín.
    n_no_apareados = {"Fe3+": 5, "Fe2+": 4, "Cu2+": 1, "Cu+": 0}.get(metal, 5)
    mult = n_m * n_no_apareados + 1
    info = {"metal_idx": metal_idxs, "odon": odon, "carga": carga,
            "multiplicidad": mult, "fe_o": fe_o, "mm": mm,
            "choques": _clashes(mol), "n_metales": n_m, "metal": metal,
            "modo": "polinuclear", "n_ligandos": n_m, "relajado": relajado,
            "n_aguas": 0, "saturado": False, "jahn_teller": False}
    return mol, info


# ── Exportación de la geometría del complejo ──────────────────────────────────

def complejo_xyz(mol, comentario=""):
    conf = mol.GetConformer()
    lns = [str(mol.GetNumAtoms()), comentario]
    for a in mol.GetAtoms():
        p = conf.GetAtomPosition(a.GetIdx())
        lns.append(f"{a.GetSymbol():<2s} {p.x:>14.8f} {p.y:>14.8f} {p.z:>14.8f}")
    return "\n".join(lns) + "\n"


def complejo_xyz_lineas(mol):
    conf = mol.GetConformer()
    out = []
    for a in mol.GetAtoms():
        p = conf.GetAtomPosition(a.GetIdx())
        out.append(f"{a.GetSymbol():<2s} {p.x:>14.8f} {p.y:>14.8f} {p.z:>14.8f}")
    return out


def complejo_pdb(mol):
    """Bloque PDB del confórmero para el visor 3D (los enlaces dativos se omiten)."""
    try:
        return Chem.MolToPDBBlock(Chem.Mol(mol))
    except Exception:
        conf = mol.GetConformer()
        out = []
        for i, a in enumerate(mol.GetAtoms(), 1):
            p = conf.GetAtomPosition(a.GetIdx())
            out.append(f"HETATM{i:>5d} {a.GetSymbol():<3s} LIG A   1    "
                       f"{p.x:8.3f}{p.y:8.3f}{p.z:8.3f}  1.00  0.00          "
                       f"{a.GetSymbol():>2s}")
        out.append("END")
        return "\n".join(out) + "\n"
