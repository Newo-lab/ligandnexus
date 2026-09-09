"""
LigandNexus — cribado virtual de derivados moleculares a partir de PubChem.

Paquete organizado por responsabilidad:
  pubchem    · cliente PUG-REST (con control de tasa y reintentos)
  sources    · registro de bases de datos intercambiables
  chem       · RDKit: núcleos, donadores coordinantes, sustituyentes
  metales    · complejos Fe/Cu–polifenol (catecol, galoílo, salicilato,
               hidroxi-cetona de flavonoides)
  screening  · filtros + selección curada (modo coordinación)
  geometry   · geometrías 3D + inputs Gaussian/ORCA/Psi4
  excel      · exportación a Excel
  admet      · descriptores y reglas ADMET (Lipinski, Veber, logBB, QED, SA)
"""

__version__ = "0.2.5"

from . import pubchem, sources, chem, metales, screening, geometry, excel, admet  # noqa: F401
