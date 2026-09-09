"""
chem.py — Utilidades de RDKit para LigandNexus
==============================================
Percepción de núcleos coordinantes, identificación de átomos donadores
(coordinantes), extracción y clasificación de sustituyentes.

Diseño clave: se trabaja SIEMPRE sobre la molécula con la aromaticidad ya
percibida por RDKit (MolFromSmiles), porque PubChem (2025) entrega los SMILES en
forma Kekulé y un patrón alifático no casaría con los anillos aromáticos.
"""

from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")


# ── Núcleo / subestructura ────────────────────────────────────────────────────

def mol_from_smiles(smiles):
    """Mol de RDKit o None. Acepta cualquier cosa convertible a str."""
    if not smiles:
        return None
    return Chem.MolFromSmiles(str(smiles))


def query_nucleo(patron):
    """
    Mol-query del núcleo coordinante, robusto. Prefiere MolFromSmiles (percibe
    aromaticidad); cae a MolFromSmarts si el patrón es un SMARTS con comodines.
    """
    if not patron:
        return None
    m = Chem.MolFromSmiles(str(patron))
    if m is not None:
        return m
    return Chem.MolFromSmarts(str(patron))


def contiene_nucleo(smiles, core_query) -> bool:
    m = mol_from_smiles(smiles)
    return bool(m and core_query is not None and m.HasSubstructMatch(core_query))


# ── Identificación de átomos coordinantes (donadores) ─────────────────────────
#
# Heurística química para señalar qué heteroátomos pueden coordinar un metal.
# La distinción crítica es N tipo piridina (par libre disponible, buen donador)
# vs N tipo pirrol (par libre dentro del sistema aromático, mal donador salvo
# desprotonación) y el N de amida (par deslocalizado hacia el C=O).

# Niveles de capacidad donadora (para ordenar y colorear).
DONOR_FUERTE = "fuerte"
DONOR_MODERADO = "moderado"
DONOR_DEBIL = "débil"
NO_DONOR = "no donador"

_ORDEN_NIVEL = {DONOR_FUERTE: 0, DONOR_MODERADO: 1, DONOR_DEBIL: 2, NO_DONOR: 3}


def _vecino_carbonilo(atom) -> bool:
    """True si `atom` está unido a un carbono que tiene un doble enlace a O (C=O)."""
    for nb in atom.GetNeighbors():
        if nb.GetAtomicNum() == 6:
            for b in nb.GetBonds():
                if b.GetBondType() == Chem.BondType.DOUBLE:
                    otro = b.GetOtherAtom(nb)
                    if otro.GetAtomicNum() == 8:
                        return True
    return False


def _n_oxigenos(atom) -> int:
    """Cuántos oxígenos cuelgan de `atom` (para reconocer nitro, sulfonilo…)."""
    return sum(1 for nb in atom.GetNeighbors() if nb.GetAtomicNum() == 8)


def _unido_a_sulfonilo(atom) -> bool:
    """True si `atom` cuelga de un S muy oxidado (sulfonamida, sulfonato)."""
    return any(nb.GetAtomicNum() == 16 and _n_oxigenos(nb) >= 2
               for nb in atom.GetNeighbors())


def clasificar_donador(atom) -> tuple:
    """
    Clasifica un átomo como donador. Devuelve (es_candidato, nivel, etiqueta).
    Solo N, O y S se consideran donadores típicos en química de coordinación.
    """
    Z = atom.GetAtomicNum()
    nH = atom.GetTotalNumHs()
    carga = atom.GetFormalCharge()
    aromatico = atom.GetIsAromatic()
    grado = atom.GetDegree()
    en_anillo = atom.IsInRing()

    # ---- NITRÓGENO ----
    if Z == 7:
        # Nitrógenos SIN par libre disponible: nitro (-NO₂), N-óxido y amonio
        # cuaternario. El N queda formalmente positivo con la valencia saturada,
        # así que NO coordina — antes caían en el cajón de «amina sp³ fuerte» y un
        # simple nitrobenceno salía como donador fuerte.
        if carga > 0 and nH == 0:
            n_oxi = _n_oxigenos(atom)
            if n_oxi >= 2:
                return True, NO_DONOR, "N de nitro (sin par libre)"
            if n_oxi == 1:
                return True, NO_DONOR, "N-óxido (sin par libre)"
            if grado >= 4:
                return True, NO_DONOR, "N amonio cuaternario (sin par libre)"
        if aromatico:
            # Piridina: N aromático, sin H, dos enlaces en el anillo → par libre libre.
            if nH == 0 and grado == 2:
                return True, DONOR_FUERTE, "N piridínico (aromático)"
            # Pirrol: N aromático con H o con tres enlaces → par libre en el anillo.
            return True, DONOR_DEBIL, "N pirrólico (par libre en el anillo)"
        if _vecino_carbonilo(atom):
            # Amida: par libre deslocalizado hacia el C=O.
            return True, DONOR_DEBIL, "N amídico (deslocalizado)"
        if _unido_a_sulfonilo(atom):
            # Sulfonamida: el par se deslocaliza hacia el S(=O)₂.
            return True, DONOR_DEBIL, "N sulfonamídico (deslocalizado)"
        # Verificar imina (N=C) mirando enlaces dobles.
        for b in atom.GetBonds():
            if b.GetBondType() == Chem.BondType.DOUBLE and \
               b.GetOtherAtom(atom).GetAtomicNum() == 6:
                return True, DONOR_MODERADO, "N imínico (C=N)"
        # Azo (-N=N-): donador flojo, no una amina sp³.
        for b in atom.GetBonds():
            if b.GetBondType() == Chem.BondType.DOUBLE and \
               b.GetOtherAtom(atom).GetAtomicNum() == 7:
                return True, DONOR_DEBIL, "N azo (-N=N-)"
        # Nitrilo (N#C).
        for b in atom.GetBonds():
            if b.GetBondType() == Chem.BondType.TRIPLE:
                return True, DONOR_DEBIL, "N nitrilo (C≡N)"
        # Amina sp3 (primaria/secundaria/terciaria).
        return True, DONOR_FUERTE, "N amínico (sp³)"

    # ---- OXÍGENO ----
    if Z == 8:
        # Los O de un nitro o de un sulfonilo/sulfonato tienen el par muy
        # repartido por el grupo: su capacidad donadora es marginal y no debe
        # confundirse con la de un carboxilato o un fenolato.
        for nb in atom.GetNeighbors():
            if nb.GetAtomicNum() == 7 and _n_oxigenos(nb) >= 2:
                return True, DONOR_DEBIL, "O de nitro (deslocalizado)"
            if nb.GetAtomicNum() == 16:
                n_ox = _n_oxigenos(nb)
                if n_ox >= 2:
                    return True, DONOR_DEBIL, "O sulfonilo/sulfonato"
                if n_ox == 1:
                    return True, DONOR_MODERADO, "O sulfinilo (S=O)"
        if carga < 0:
            return True, DONOR_FUERTE, "O carboxilato/alcóxido (aniónico)"
        # Carbonilo (=O) frente a O de tipo éter/hidroxilo.
        doble_o = any(b.GetBondType() == Chem.BondType.DOUBLE for b in atom.GetBonds())
        if doble_o:
            return True, DONOR_MODERADO, "O carbonílico (C=O)"
        if nH >= 1:
            return True, DONOR_MODERADO, "O hidroxilo/fenol"
        return True, DONOR_DEBIL, "O éter"

    # ---- AZUFRE ----
    if Z == 16:
        if carga < 0 or nH >= 1:
            return True, DONOR_FUERTE, "S tiol/tiolato"
        # Sulfona/sulfóxido: S muy oxidado, mal donador por S.
        n_oxi = sum(1 for nb in atom.GetNeighbors() if nb.GetAtomicNum() == 8)
        if n_oxi >= 1:
            return True, NO_DONOR, "S oxidado (sulfonilo/sulfinilo)"
        return True, DONOR_MODERADO, "S tioéter"

    # Fósforo (fosfinas) como extra.
    if Z == 15:
        return True, DONOR_MODERADO, "P fosfina"

    return False, NO_DONOR, ""


def identificar_donadores(smiles, solo_relevantes=True) -> dict:
    """
    Identifica los átomos coordinantes de una molécula.

    Devuelve un dict:
      {
        "mol": mol,
        "donadores": [ {idx, elemento, nivel, etiqueta, aromatico, en_anillo}, ... ],
        "resumen": {"fuerte": n, "moderado": n, "débil": n, "no donador": n},
      }
    Los donadores se devuelven ordenados por nivel (más fuerte primero) y luego
    por índice, que es el orden natural para una numeración posterior.
    """
    mol = mol_from_smiles(smiles)
    if mol is None:
        return {"mol": None, "donadores": [], "resumen": {}}

    donadores = []
    resumen = {DONOR_FUERTE: 0, DONOR_MODERADO: 0, DONOR_DEBIL: 0, NO_DONOR: 0}
    for atom in mol.GetAtoms():
        es_cand, nivel, etiqueta = clasificar_donador(atom)
        if not es_cand:
            continue
        resumen[nivel] += 1
        if solo_relevantes and nivel == NO_DONOR:
            continue
        donadores.append({
            "idx": atom.GetIdx(),
            "elemento": atom.GetSymbol(),
            "nivel": nivel,
            "etiqueta": etiqueta,
            "aromatico": atom.GetIsAromatic(),
            "en_anillo": atom.IsInRing(),
        })
    donadores.sort(key=lambda d: (_ORDEN_NIVEL[d["nivel"]], d["idx"]))
    return {"mol": mol, "donadores": donadores, "resumen": resumen}


def estimar_sitios(mol, donadores) -> int:
    """
    Estima el número de SITIOS de coordinación (denticidad probable) a partir de
    los donadores fuertes/moderados. Corrige el conteo de los grupos con varios
    oxígenos sobre un mismo átomo central (carboxilato, nitro, sulfonato,
    fosfato): esos O son un único sitio, porque el grupo coordina por uno de
    ellos o puentea, pero cuenta como una sola posición de anclaje.
    """
    utiles = [d for d in donadores if d["nivel"] in (DONOR_FUERTE, DONOR_MODERADO)]
    idx_util = {d["idx"] for d in utiles}
    contados = set()
    sitios = 0
    for d in utiles:
        if d["idx"] in contados:
            continue
        contados.add(d["idx"])
        sitios += 1
        atom = mol.GetAtomWithIdx(d["idx"])
        # Si es un O unido a un átomo central que lleva otro O donador, ambos son
        # el mismo sitio. Antes solo se miraba el CARBONO, así que los oxígenos de
        # un nitro o de un sulfonato se contaban por separado pese a que la
        # docstring prometía corregirlo.
        if atom.GetAtomicNum() == 8:
            for c in atom.GetNeighbors():
                if c.GetAtomicNum() in (6, 7, 15, 16):
                    for o in c.GetNeighbors():
                        if (o.GetAtomicNum() == 8 and o.GetIdx() in idx_util
                                and o.GetIdx() != d["idx"]):
                            contados.add(o.GetIdx())
    return sitios


# Colores por nivel para resaltar en la imagen 2D (RGB 0-1).
_COLOR_NIVEL = {
    DONOR_FUERTE:   (0.13, 0.55, 0.13),   # verde
    DONOR_MODERADO: (0.85, 0.55, 0.10),   # ámbar
    DONOR_DEBIL:    (0.80, 0.30, 0.30),   # rojo apagado
    NO_DONOR:       (0.60, 0.60, 0.60),   # gris
}


def imagen_donadores(smiles, size=(440, 340), numerar=True, dark=False):
    """
    Imagen 2D con los átomos coordinantes resaltados por color según su fuerza
    (verde=fuerte, ámbar=moderado, rojo=débil) y numerados en el orden de
    `identificar_donadores` (útil para la numeración posterior de posiciones).
    Con dark=True usa paleta clara sobre fondo transparente (modo oscuro).
    Devuelve (PIL.Image, lista_donadores) o (None, []) si el SMILES es inválido.
    """
    from io import BytesIO
    from PIL import Image
    from rdkit.Chem.Draw import rdMolDraw2D

    info = identificar_donadores(smiles)
    mol = info["mol"]
    if mol is None:
        return None, []
    donadores = info["donadores"]

    highlight = [d["idx"] for d in donadores]
    colores = {d["idx"]: _COLOR_NIVEL[d["nivel"]] for d in donadores}

    if numerar:
        for n, d in enumerate(donadores, 1):
            mol.GetAtomWithIdx(d["idx"]).SetProp("atomNote", str(n))

    d2d = rdMolDraw2D.MolDraw2DCairo(*size)
    if dark:
        rdMolDraw2D.SetDarkMode(d2d)
    opts = d2d.drawOptions()
    opts.setBackgroundColour((0, 0, 0, 0))      # transparente: se integra al tema
    opts.highlightBondWidthMultiplier = 12
    opts.annotationFontScale = 0.9
    rdMolDraw2D.PrepareAndDrawMolecule(
        d2d, mol, highlightAtoms=highlight, highlightAtomColors=colores)
    d2d.FinishDrawing()
    img = Image.open(BytesIO(d2d.GetDrawingText()))
    return img, donadores


# ── Clasificación de sustituyentes (para el modo de coordinación) ─────────────

GF_PATRONES = [
    ("CF3", "[CX4](F)(F)F"), ("CHF2", "[CX4](F)(F)"),
    ("Halo-F", "[F]"), ("Halo-Cl", "[Cl]"), ("Halo-Br", "[Br]"), ("Halo-I", "[I]"),
    ("N-oxido", "[N+][O-]"), ("Nitro", "[N+](=O)[O-]"), ("Nitroso", "[NX2]=O"),
    ("Cyano", "C#N"), ("Isociano", "[N+]#[C-]"),
    ("Sulfonyl", "[SX4](=O)(=O)"), ("Sulfinyl", "[SX3](=O)"),
    ("Tiol", "[SX2H]"), ("Sulfide", "[SX2]"),
    ("Amida", "C(=O)[NH]"), ("Ester", "C(=O)OC"), ("Carboxilo", "C(=O)O"),
    ("Cetona", "C(=O)[C]"), ("Aldehido", "[CH]=O"), ("Carbonilo", "C=O"),
    ("Metoxi", "[OX2][CH3]"), ("Alcoxi", "[OX2][CX4]"), ("Hidroxi", "[OX2H]"),
    ("Oxo", "[OX1]"),
    ("Amino-NH2", "[NH2]"), ("Amino-NH", "[NH]"), ("Amino-N", "[NX3;!H]"),
    ("Imino", "[NX2]=C"),
    ("Alquino", "C#C"), ("Alqueno", "C=C"),
    ("Cicloalquilo", "[CX4;R]"), ("Metilo", "[CH3]"), ("Alquilo", "[CX4]"),
]


def compilar_gf(gf_patrones):
    out = []
    for n, p in gf_patrones:
        m = Chem.MolFromSmarts(p)
        if m is not None:
            out.append((n, m))
    return out


_GF_SMARTS_MOL = compilar_gf(GF_PATRONES)


def clasificar_gf_desde_mol(mol, frag_atoms, gf_compilados=None):
    """Clasifica el grupo funcional de un fragmento usando la aromaticidad ya
    percibida en la molécula original (evita SMILES de fragmento inválidos)."""
    if gf_compilados is None:
        gf_compilados = _GF_SMARTS_MOL
    if not frag_atoms:
        return "H"

    atoms = [mol.GetAtomWithIdx(a) for a in frag_atoms]
    ans = [a.GetAtomicNum() for a in atoms]

    if len(frag_atoms) == 1:
        atom = atoms[0]
        an = atom.GetAtomicNum()
        nH = atom.GetTotalNumHs()
        if an == 9:  return "Halo-F"
        if an == 17: return "Halo-Cl"
        if an == 35: return "Halo-Br"
        if an == 53: return "Halo-I"
        if an == 8:  return "Hidroxi" if nH >= 1 else "Oxo"
        if an == 7:
            if nH >= 2: return "Amino-NH2"
            if nH == 1: return "Amino-NH"
            return "Amino-N"
        if an == 16: return "Tiol" if nH >= 1 else "Sulfide"
        if an == 6:
            if nH >= 3: return "Metilo"
            if nH == 2: return "Metileno"
            return "Alquilo"
        return f"Atomo-{atom.GetSymbol()}"

    arom = [a for a in atoms if a.GetIsAromatic()]
    if arom:
        has_het = any(a.GetAtomicNum() != 6 for a in arom)
        return "Heteroarilo" if has_het else "Arilo"

    frag_set = set(frag_atoms)
    for nombre, patt in gf_compilados:
        for match in mol.GetSubstructMatches(patt):
            if frag_set.intersection(match):
                return nombre

    an_set = set(ans)
    if an_set <= {9}:  return "Halo-F"
    if an_set <= {17}: return "Halo-Cl"
    if an_set <= {35}: return "Halo-Br"
    if an_set <= {53}: return "Halo-I"
    if an_set <= {8}:  return "Oxo"
    if an_set <= {7}:  return "N-imino"
    if an_set <= {6}:  return "Alquilo"
    if 6 in an_set and 8 in an_set and 7 not in an_set: return "C-O"
    if 6 in an_set and 7 in an_set:  return "C-N"
    if 6 in an_set and 16 in an_set: return "C-S"
    return "Mixto"


def _pos_label(dist):
    return {0: "alpha (en coord.)", 1: "alpha (1 enlace)", 2: "beta (2 enlaces)",
            3: "gamma (3 enlaces)", 4: "delta (4 enlaces)"}.get(dist, f"{dist} enlaces")


def extraer_sustituyentes(smiles_mol, smarts_core, gf_compilados=None) -> list:
    """
    Lista de sustituyentes del núcleo. Para cada uno: sub_smiles, n_atomos_sub,
    gf, dist_coord (enlaces al donador más cercano), pos_label, core_atom_idx.
    Los donadores del núcleo se toman de la percepción química (N/O), coherente
    con `identificar_donadores`.
    """
    mol = mol_from_smiles(smiles_mol)
    if mol is None:
        return []
    core = query_nucleo(smarts_core)
    if core is None:
        return []
    match = mol.GetSubstructMatch(core)
    if not match:
        return []

    core_set = set(match)
    coord_set = {
        match[i] for i in range(len(match))
        if mol.GetAtomWithIdx(match[i]).GetAtomicNum() in (7, 8)
    }

    resultado, seen = [], set()
    for core_smarts_idx, core_mol_idx in enumerate(match):
        core_atom = mol.GetAtomWithIdx(core_mol_idx)
        for neighbor in core_atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            if n_idx in core_set or n_idx in seen:
                continue
            # BFS del fragmento sustituyente.
            visited, queue, frag_set = set(), [n_idx], set()
            while queue:
                curr = queue.pop(0)
                if curr in visited:
                    continue
                visited.add(curr)
                frag_set.add(curr)
                for nb in mol.GetAtomWithIdx(curr).GetNeighbors():
                    if nb.GetIdx() not in core_set and nb.GetIdx() not in visited:
                        queue.append(nb.GetIdx())
            seen.update(frag_set)

            sub_smiles = Chem.MolFragmentToSmiles(
                mol, list(frag_set), isomericSmiles=True, allHsExplicit=False)

            if coord_set:
                dists = []
                for c_idx in coord_set:
                    if core_mol_idx == c_idx:
                        dists.append(1)
                    else:
                        path = Chem.GetShortestPath(mol, core_mol_idx, c_idx)
                        dists.append(len(path) - 1)
                dist = min(dists)
            else:
                dist = -1

            gf = clasificar_gf_desde_mol(mol, frag_set, gf_compilados)
            resultado.append({
                "sub_smiles": sub_smiles, "n_atomos_sub": len(frag_set),
                "gf": gf, "dist_coord": dist, "pos_label": _pos_label(dist),
                "core_atom_idx": core_smarts_idx,
            })
    return resultado
