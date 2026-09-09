"""
geometry.py — Geometrías 3D e inputs de química cuántica
========================================================
Genera una geometría 3D pre-optimizada con campo de fuerzas (MMFF94, fallback
UFF) como punto de partida — NO una optimización DFT — y arma los archivos de
entrada para Gaussian, ORCA y Psi4.
"""

from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import AllChem

from .chem import mol_from_smiles


def _params_etkdg(semilla):
    p = AllChem.ETKDGv3()
    p.randomSeed = semilla
    p.numThreads = 0
    p.maxIterations = 1000
    p.pruneRmsThresh = 0.5
    return p


def generar_3d(smiles, n_confs=10, semilla=42, max_iter=2000):
    """
    Geometría 3D para un SMILES. Con n_confs>1 explora varios confórmeros y
    devuelve el de menor energía.

    Devuelve (mol_con_H, conf_id, metodo, energia_kcalmol, n_confs_efectivos)
    o (None, None, motivo, None, 0) si falla. `n_confs_efectivos` es el número
    de confórmeros realmente generados tras la poda por RMSD (puede ser < n_confs).
    """
    mol = mol_from_smiles(smiles)
    if mol is None:
        return None, None, "SMILES_invalido", None, 0
    mol = Chem.AddHs(mol)

    cids = list(AllChem.EmbedMultipleConfs(
        mol, numConfs=max(1, int(n_confs)), params=_params_etkdg(semilla)))
    if not cids:
        p2 = AllChem.ETKDGv3()
        p2.randomSeed = semilla + 999
        p2.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, p2) != 0:
            return None, None, "embedding_fallido", None, 0
        cids = [0]

    if AllChem.MMFFHasAllMoleculeParams(mol):
        metodo = "MMFF94"
        res = AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=max_iter, numThreads=0)
    else:
        metodo = "UFF"
        res = AllChem.UFFOptimizeMoleculeConfs(mol, maxIters=max_iter, numThreads=0)

    conf_ids = [c.GetId() for c in mol.GetConformers()]
    datos = [(ener, cid_, (flag == 0)) for (flag, ener), cid_ in zip(res, conf_ids)]
    convergidos = [d for d in datos if d[2]]
    pool = convergidos if convergidos else datos
    best_e, best_cid, best_conv = min(pool, key=lambda x: x[0])
    if not best_conv:
        metodo += "_no_convergido"

    return mol, best_cid, metodo, float(best_e), len(conf_ids)


def xyz_lineas(mol, conf_id):
    conf = mol.GetConformer(conf_id)
    out = []
    for atom in mol.GetAtoms():
        p = conf.GetAtomPosition(atom.GetIdx())
        out.append(f"{atom.GetSymbol():<2s} {p.x:>14.8f} {p.y:>14.8f} {p.z:>14.8f}")
    return out


def mol_a_xyz(mol, conf_id, comentario=""):
    lns = xyz_lineas(mol, conf_id)
    return f"{mol.GetNumAtoms()}\n{comentario}\n" + "\n".join(lns) + "\n"


def mol_a_molblock(mol, conf_id):
    return Chem.MolToMolBlock(mol, confId=conf_id)


def mol_a_sdf(mol, conf_id):
    return Chem.MolToMolBlock(mol, confId=conf_id) + "$$$$\n"


def carga_formal(mol):
    return Chem.GetFormalCharge(mol)


def pdb_block(mol, conf_id):
    """Bloque PDB del confórmero (para el visor 3D interactivo del navegador)."""
    return Chem.MolToPDBBlock(mol, confId=conf_id)


# ── Generadores de input ──────────────────────────────────────────────────────

def construir_input_gaussian(nombre, carga, mult, lineas_xyz, *, mem="8GB",
                             nprocshared=8, usar_chk=True, metodo="M06-2X",
                             base="6-311+G(d,p)", keywords="", titulo=None,
                             secciones_extra=""):
    route = f"#p {metodo}/{base}".rstrip()
    if keywords.strip():
        route += " " + keywords.strip()
    out = []
    if usar_chk:
        out.append(f"%chk={nombre}.chk")
    out += [f"%mem={mem}", f"%nprocshared={int(nprocshared)}", route, "",
            titulo or nombre, "", f"{int(carga)} {int(mult)}"]
    out += list(lineas_xyz)
    if secciones_extra.strip():
        out.append("")
        out.append(secciones_extra.rstrip("\n"))
    out.append("")
    return "\n".join(out) + "\n"


def construir_input_orca(nombre, carga, mult, lineas_xyz, *, metodo="M06-2X",
                         base="def2-TZVP", keywords="TightSCF", nprocs=8,
                         maxcore=3000, bloques_extra=""):
    linea = f"! {metodo} {base} {keywords}".rstrip()
    out = [f"# {nombre}", linea, f"%pal nprocs {int(nprocs)} end",
           f"%maxcore {int(maxcore)}"]
    if bloques_extra.strip():
        out.append(bloques_extra.rstrip("\n"))
    out.append(f"* xyz {int(carga)} {int(mult)}")
    out += list(lineas_xyz)
    out.append("*")
    out.append("")
    return "\n".join(out) + "\n"


def construir_input_psi4(nombre, carga, mult, lineas_xyz, *, memoria="8 GB",
                         base="6-311+G(d,p)", metodo="m06-2x", extra=""):
    out = [f"# {nombre}", f"memory {memoria}", "", "molecule {",
           f"{int(carga)} {int(mult)}"]
    out += list(lineas_xyz)
    out += ["}", "", f"set basis {base}"]
    if extra.strip():
        out.append(extra.rstrip("\n"))
    out.append(f"energy('{metodo}')")
    out.append("")
    return "\n".join(out) + "\n"
