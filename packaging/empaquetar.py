# -*- coding: utf-8 -*-
"""
empaquetar.py — arma la carpeta distribuible de LigandNexus v0.2.5 y su .zip.

Copia el proyecto SIN el entorno virtual, sin cachés, sin los backups y sin las
carpetas de trabajo de las auditorías. Incluye `tests/` a propósito: las 259
comprobaciones empíricas son parte de lo que se entrega.
"""
import pathlib
import shutil
import zipfile

# El origen es la raíz del propio proyecto (este archivo vive en packaging/),
# de modo que el script funciona en cualquier máquina y no filtra una ruta
# personal al repositorio público.
ORIGEN = pathlib.Path(__file__).resolve().parent.parent
ESCRITORIO = ORIGEN.parent
VERSION = "0.2.5"
DESTINO = ESCRITORIO / f"LigandNexus_v{VERSION}"
ZIP = ESCRITORIO / f"LigandNexus_v{VERSION}.zip"

# Lo que NO viaja.
EXCLUIR_DIR = {".venv", "__pycache__", ".git", ".pytest_cache", ".ipynb_checkpoints"}
EXCLUIR_PREFIJO = ("_trabajo_", "_backup_", "_DESCARTADO_")
# `.pre_*` es como se hacen aquí las copias antes de editar (main.tex.pre_x,
# README.md.pre_repo…). Sin esta línea viajaban al paquete y al .zip.
EXCLUIR_SUFIJO = (".pyc", ".pyo", ".log", ".bak", ".orig", ".rej")
EXCLUIR_INFIJO = (".pre_",)
# Capturas de trabajo del empaquetado anterior: no aportan al usuario final.
EXCLUIR_ARCHIVO = {"run_live.log"}


def se_excluye(ruta: pathlib.Path) -> bool:
    for parte in ruta.parts:
        if parte in EXCLUIR_DIR or parte.startswith(EXCLUIR_PREFIJO):
            return True
    if any(x in ruta.name for x in EXCLUIR_INFIJO):
        return True
    if ruta.name in EXCLUIR_ARCHIVO or ruta.name.endswith(EXCLUIR_SUFIJO):
        return True
    # Las capturas de pantalla de la verificación no van en la entrega.
    if ruta.parts[:1] == ("packaging",) and ruta.name.startswith(("_shot", "_test")):
        return True
    return False


if DESTINO.exists():
    shutil.rmtree(DESTINO)
DESTINO.mkdir(parents=True)

n = 0
for f in sorted(ORIGEN.rglob("*")):
    if not f.is_file():
        continue
    rel = f.relative_to(ORIGEN)
    if se_excluye(rel):
        continue
    salida = DESTINO / rel
    salida.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(f, salida)
    n += 1

print(f"copiados {n} archivos a {DESTINO.name}")

# El LEEME vive en packaging/ dentro del repo, pero en la entrega va a la raíz,
# que es donde lo va a buscar quien reciba el zip.
leeme = DESTINO / "packaging" / "LEEME_PRIMERO.txt"
assert leeme.exists(), "falta packaging/LEEME_PRIMERO.txt"
shutil.move(str(leeme), str(DESTINO / "LEEME_PRIMERO.txt"))
print("  LEEME_PRIMERO.txt movido a la raiz del paquete")

# ── Comprobaciones de que el paquete NO sale roto ─────────────────────────────
bat = DESTINO / "Iniciar_LigandNexus.bat"
crudo = bat.read_bytes()
escapes = crudo.count(b"^(")
puerto = crudo.lower().count(b"reinici") + crudo.lower().count(b"ocupad")
assert escapes >= 4, f"el .bat del paquete NO lleva los parentesis escapados ({escapes})"
assert puerto >= 2, "el .bat del paquete no trae la logica de puerto ocupado"
assert b"\r\n" in crudo, "el .bat debe ir con finales de linea CRLF"
try:
    crudo.decode("ascii")
except UnicodeDecodeError as e:
    raise AssertionError(f"el .bat debe ir en ASCII puro: {e}")
print(f"  .bat verificado: {escapes} escapes '^(', CRLF, ASCII puro")

for imprescindible in ("app.py", "requirements.txt", "LICENSE", "README.md",
                       "ligandnexus/__init__.py", "views/hierro.py",
                       "tests/verificar_motores.py", ".streamlit/config.toml"):
    assert (DESTINO / imprescindible).exists(), f"falta {imprescindible}"
assert not (DESTINO / ".venv").exists(), "se colo el entorno virtual"
assert not list(DESTINO.rglob("__pycache__")), "se colaron cachés"
respaldos = [str(f.relative_to(DESTINO)) for f in DESTINO.rglob("*")
             if f.is_file() and (".pre_" in f.name or f.suffix in (".bak", ".orig"))]
assert not respaldos, f"se colaron archivos de respaldo: {respaldos[:3]}"
version = (DESTINO / "ligandnexus" / "__init__.py").read_text(encoding="utf-8")
assert f'__version__ = "{VERSION}"' in version, "la version del paquete no es la esperada"
print("  estructura verificada")

# ── .zip ──────────────────────────────────────────────────────────────────────
if ZIP.exists():
    ZIP.unlink()
with zipfile.ZipFile(ZIP, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as z:
    for f in sorted(DESTINO.rglob("*")):
        if f.is_file():
            z.write(f, pathlib.Path(DESTINO.name) / f.relative_to(DESTINO))

with zipfile.ZipFile(ZIP) as z:
    malo = z.testzip()
    assert malo is None, f"zip corrupto en {malo}"
    nombres = z.namelist()
    restos = [x for x in nombres if "__pycache__" in x or x.endswith(".pyc")]
    assert not restos, f"se colaron cachés en el zip: {restos[:3]}"
print(f"  {ZIP.name}: {len(nombres)} archivos, {ZIP.stat().st_size/1024:.0f} KB, integridad OK")


def limpiar_cache():
    """
    Borra los `__pycache__` de la carpeta entregable.

    Hay que llamarlo DESPUÉS de correr la suite desde dentro del paquete: al
    importar `ligandnexus` y `views` desde ahí, Python deja los .pyc en la
    carpeta, que ya estaba copiada limpia y con el zip cerrado. No rompe nada
    (Python los regenera), pero la carpeta que se entrega no debe llevarlos.
    """
    n = 0
    for d in sorted(DESTINO.rglob("__pycache__"), reverse=True):
        shutil.rmtree(d, ignore_errors=True)
        n += 1
    for f in DESTINO.rglob("*.pyc"):
        f.unlink(missing_ok=True)
    return n


if __name__ == "__main__":
    print("\nSi vas a correr la suite desde el paquete, llama luego a limpiar_cache().")
    print("Y OJO: si la app está corriendo DESDE la carpeta del paquete, el rmtree")
    print("inicial falla con PermissionError [WinError 32] — para el servidor antes.")
