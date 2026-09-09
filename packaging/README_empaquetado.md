# 📦 Empaquetado y distribución de LigandNexus

Tres formas de entregar la app, de la más simple a la más autónoma.

## 1. Lanzador `.bat` (recomendado para uso propio)

`Iniciar_LigandNexus.bat` (en la raíz del proyecto). La primera vez crea un entorno
local (`.venv`) e instala las dependencias; luego solo abre la app.

- **Requiere:** Python 3.10+ instalado (con «Add Python to PATH») e internet la
  primera vez.
- **Ventaja:** sencillo, ligero, fácil de mantener.
- **Nota de red institucional:** el `.bat` instala con `--trusted-host pypi.org
  --trusted-host files.pythonhosted.org` porque algunos proxys corporativos
  rompen la verificación de certificados de pip.

## 2. Ejecutable `.exe` con PyInstaller

Recetas y scripts: `packaging/run_ligandnexus.py` (arranque) y
`packaging/ligandnexus.spec` (spec).

```bash
.venv/Scripts/python -m PyInstaller packaging/ligandnexus.spec --noconfirm
# resultado: dist/LigandNexus/LigandNexus.exe  (carpeta autocontenida)
```

- **Ventaja:** no requiere que el usuario tenga Python.
- **Estado en esta máquina (2026-07-06):** el bundle se ensambla, pero el `.exe`
  final **desaparece justo tras crearse** (error `FileNotFoundError` en el paso
  `set_exe_build_timestamp`), de forma **reproducible en dos intentos**. El
  responsable es **Norton 360** (que además tiene relevada la protección de Windows
  Defender): pone en cuarentena el *bootloader* de PyInstaller — un **falso positivo
  muy conocido**. Para compilar aquí habría que añadir una **exclusión en Norton**
  para la carpeta de compilación (requiere permiso del usuario).
- **Implicación para distribuir:** aunque se genere en otra máquina, el `.exe` de
  PyInstaller suele activar alertas antivirus en los equipos de destino. Si se opta
  por esta vía, conviene **firmar el ejecutable** (certificado de firma de código)
  para reducir los falsos positivos, y añadir una **exclusión** en el equipo donde
  se compila.
- Para diagnosticar/compilar sin interferencias: compilar en una carpeta excluida del
  antivirus, o en un equipo con el antivirus temporalmente en pausa (con permiso).

## 3. Portable con Python embebido (mejor «sin instalar nada»)

La opción más fiable para entregar a quien no tiene Python **sin** los falsos
positivos de PyInstaller: una carpeta con una distribución **embeddable de Python**
(python.org → *Windows embeddable package*) + las librerías + un lanzador.

Bosquejo (`build_portable` — pendiente de automatizar):
1. Descargar `python-3.11.x-embed-amd64.zip` de python.org y descomprimir en `python/`.
2. Habilitar `site` (descomentar `import site` en `python311._pth`).
3. Instalar las dependencias en esa distribución con `get-pip.py` + `pip install -r
   requirements.txt --target python/Lib/site-packages`.
4. Copiar `app.py`, `ligandnexus/`, `views/`, `.streamlit/`, `examples/`.
5. Lanzador `.bat`/`.exe` que ejecute `python\python.exe -m streamlit run app.py`.

El resultado es un `.zip` que se descomprime y se ejecuta con doble clic, sin
instalar nada y sin internet (salvo las consultas a PubChem, que sí lo requieren).

---

### Recomendación

- **Para Owen y colaboradores con Python:** el **lanzador `.bat`** (opción 1).
- **Para el jurado o quien no tenga Python:** el **portable embebido** (opción 3),
  más fiable que el `.exe` de PyInstaller frente a los antivirus.
- El **`.exe`** (opción 2) queda como posibilidad documentada; su distribución masiva
  exige firma de código para no disparar antivirus.
