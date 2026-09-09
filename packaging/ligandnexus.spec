# -*- mode: python ; coding: utf-8 -*-
# Spec de PyInstaller para LigandNexus (app Streamlit + RDKit).
import os
from PyInstaller.utils.hooks import collect_all, copy_metadata, collect_submodules

PROJ = os.path.abspath(os.path.join(os.getcwd()))

datas, binaries, hiddenimports = [], [], []

# Recolectar paquetes completos (datos + binarios + submódulos).
for pkg in ("streamlit", "rdkit", "py3Dmol", "altair", "pandas", "openpyxl",
            "PIL", "pyarrow"):
    d, b, h = collect_all(pkg)
    datas += d; binaries += b; hiddenimports += h

# Metadatos que Streamlit consulta en tiempo de ejecución.
for pkg in ("streamlit", "rdkit", "pandas", "numpy", "altair", "pyarrow",
            "packaging", "requests", "click", "tornado"):
    try:
        datas += copy_metadata(pkg)
    except Exception:
        pass

hiddenimports += collect_submodules("streamlit")

# El código de la app y sus recursos, colocados en la raíz del bundle.
datas += [
    (os.path.join(PROJ, "app.py"), "."),
    (os.path.join(PROJ, "ligandnexus"), "ligandnexus"),
    (os.path.join(PROJ, "views"), "views"),
    (os.path.join(PROJ, ".streamlit"), ".streamlit"),
    (os.path.join(PROJ, "examples"), "examples"),
]

block_cipher = None

a = Analysis(
    [os.path.join(PROJ, "packaging", "run_ligandnexus.py")],
    pathex=[PROJ],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    runtime_hooks=[],
    excludes=["tkinter", "matplotlib", "PyQt5", "PySide2"],
    cipher=block_cipher,
    noarchive=False,
)
pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz, a.scripts, [], exclude_binaries=True,
    name="LigandNexus", debug=False, bootloader_ignore_signals=False,
    strip=False, upx=False, console=True,
)
coll = COLLECT(
    exe, a.binaries, a.zipfiles, a.datas,
    strip=False, upx=False, name="LigandNexus",
)
