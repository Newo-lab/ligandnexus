# 🔗 LigandNexus

> Cribado virtual de derivados moleculares a partir de PubChem — versión **0.2.5**

**LigandNexus** recupera los **derivados** de una molécula base desde **PubChem**, los
depura y los prepara para estudios computacionales. Está pensado para usarse **sin
programar**: la interfaz web expone todo como controles.

La versión 0.2 reorganiza la herramienta en **cinco modos independientes**, para que
cada quien use solo lo que necesita:

| Modo | Para qué sirve |
|------|----------------|
| 🩸 **Metales + polifenoles** | Coordinación de Fe(II)/Fe(III) y Cu(I)/Cu(II) por **taninos y flavonoides**: detecta los sitios quelantes O,O, estima la especiación mono/bis/tris frente al pH y construye el complejo 3D **mono o polinuclear** como semilla para DFT. |
| 🔎 **Cribado general** | Recuperar derivados de PubChem con sus propiedades, filtrarlos por rango y exportarlos a CSV/Excel. Sin suponer química de coordinación. |
| 🧭 **Análisis de coordinación** | Pipeline completo para ligandos: depurar, clasificar sustituyentes por grupo funcional y **posición respecto al átomo coordinante**, y seleccionar un conjunto representativo con balance donador/aceptor. |
| 🎯 **Átomos coordinantes** | Identificar y **numerar los átomos donadores** (N, O, S, P) de una molécula, coloreados por fuerza donadora. Base para numerar posiciones. |
| 🧬 **Geometrías 3D** | Pre-optimizar geometrías con campo de fuerzas, verlas en un **visor 3D interactivo** y generar inputs de Gaussian / ORCA / Psi4. |

## 🔧 Correcciones de la v0.2.5

Auditoría de la pestaña **Exportar** tras la implementación de ADMET. Los motores de
química estaban sanos (113/113 comprobaciones); los fallos vivían en el **cableado de la
interfaz**, que hasta ahora no tenía pruebas. Ahora sí: la suite pasa de 113 a **259**
comprobaciones — el bloque `[12]` conduce las vistas reales con `AppTest` y las fuentes
simuladas (sin red), el `[13]` contrasta el motor de metales contra un artículo externo,
y los `[14]`–`[15]` verifican los polifenoles de ejemplo y el caso Fe(III) + ácido elágico.

Después llegó una **auditoría de seguridad**, y el bloque `[16]` es su resultado: 28
comprobaciones que **reproducen cada ataque** y exigen que el parche lo neutralice. Se
validaron reintroduciendo los fallos a propósito, y la primera versión de esas pruebas
tenía a su vez un hueco — comprobaba que el saneador funcionara, pero no que el visor
llegara a llamarlo, que es exactamente el tipo de descuido que causó los bugs del `[12]`.

- **Lo que se exporta ya no puede ser de otra molécula.** El perfil ADMET se guardaba en
  sesión y no se borraba nunca, así que al cambiar de molécula la pestaña Exportar seguía
  ofreciéndolo —por defecto, y a veces como única opción— con los derivados anteriores.
- **La descarga coincide con lo que se ve.** Un filtro ADMET que no dejaba ninguna
  molécula se descartaba en silencio: la pantalla decía «0 de N» y la descarga entregaba
  las N. El rótulo dice además de qué etapa del cribado sale el perfil.
- **El botón «Descargar Excel» ya no desaparece.** Vivía dentro del `if st.button(...)`,
  de modo que el propio clic de descarga lo borraba. Ahora el archivo se guarda en sesión
  atado al conjunto con el que se generó, y si cambian las opciones se avisa en vez de
  esfumarse sin explicación.
- **Perfilar dos veces ya no rompe el Excel.** El sufijo `_admet` se aplicaba sin mirar si
  ese nombre ya estaba ocupado; la columna duplicada hacía que `openpyxl` rechazara una
  Serie en vez de un valor.
- **Acotar un filtro ya no tira el cálculo.** Los perfiles se memorizan por SMILES (solo
  dependen de la estructura), así que rehacerlos es gratis y se rehacen solos: antes,
  mover un deslizador obligaba a repetir 5 min de cálculo en una corrida de 40.000.
- Los CSV se arman **al pulsar** y no en cada reejecución (con ADMET son 28 columnas en
  vez de 6: 178 ms frente a 56 ms cada vez, con 20.000 filas, y `st.tabs` pinta todas
  las pestañas siempre).
- **Geometrías 3D ya ve los conjuntos filtrados por ADMET.** Quien acotaba por Lipinski o
  por barrera hematoencefálica no podía llevarse *esos* candidatos a 3D: la vista solo
  ofrecía las listas sin perfilar. Ahora ofrece los dos perfiles (coordinación y cribado
  general), con el rótulo diciendo de qué etapa salen, y descarta los conjuntos vacíos.
- **El motor de metales, validado contra un artículo externo.** 24 comprobaciones nuevas
  contra Zaccaron *et al.*, *J. Coord. Chem.* **2013**, *66*, 1709–1719 (tintas
  ferrogálicas): sitios de quelación del ácido gálico, número de coordinación 6,
  distancias Fe–O cristalográficas y estados de espín. El detalle de las 24 comprobaciones
  está en el bloque `[13]` de `tests/verificar_motores.py`.
- **El ácido elágico de los ejemplos era un regioisómero.** Misma fórmula (C₁₄H₆O₈), mismo
  peso y los mismos dos catecoles, pero con los hidroxilos de un anillo corridos una
  posición: no lo delataba ni la fórmula ni el recuento de sitios, sí la **simetría** (el
  ácido elágico real tiene sus dos mitades equivalentes) y, aguas abajo, un Fe···Fe
  desplazado 0,41 Å en el complejo polinuclear. Corregido contra PubChem CID 5281855, y de
  paso el EGCG, que venía sin sus dos centros (2*R*,3*R*). **Los ocho polifenoles de
  ejemplo quedan verificados contra el SMILES de su CID.**

## ✨ Novedades de la v0.2

- **Múltiples bases de datos** (~20 a elegir; PubChem por defecto). Con búsqueda de
  subestructura **en vivo** en PubChem, ChEMBL y COCONUT, y para el resto (LOTUS,
  NPASS, NP Atlas, ChEBI, HMDB, ZINC, DrugBank…) la vía **descarga el volcado +
  búsqueda local con RDKit**. Además puedes subir **tu propio archivo** (SDF/SMILES/CSV).
  La lógica vive en `ligandnexus/sources.py` (registro + fuentes intercambiables).
- **Arquitectura modular** (`ligandnexus/` como paquete + `views/` para la interfaz).
- **Cliente de PubChem robusto**: control de tasa (~5 peticiones/s), reintentos con
  espera creciente ante «servidor ocupado» y aviso si un lote no se descarga (la v0.1
  perdía filas en silencio).
- **Identificador de átomos coordinantes** con heurística química (distingue N
  piridínico de pirrólico, detecta amidas, clasifica O y S).
- **Visor 3D interactivo** en el navegador (py3Dmol).
- **Panel ADMET** (`ligandnexus/admet.py`), disponible en **cualquier punto del
  cribado**: descriptores fisicoquímicos y reglas de biodisponibilidad calculadas en
  local con RDKit (Lipinski, Veber, Egan, **logBB de Clark** para la barrera
  hematoencefálica, QED, SAscore y alertas PAINS), con filtros por cada criterio.
- **Complejos metal–polifenol** (`ligandnexus/metales.py`): reconoce los **cinco motivos
  quelantes O,O** —catecol, galoílo, salicilato y las dos **hidroxi-cetonas** de los
  flavonoides (3-OH/4-C=O y 5-OH/4-C=O)— y arma la geometría sobre el politopo
  cristalográfico del metal (M–O tabulada, mordida derivada del tamaño del anillo
  quelato, Jahn-Teller en Cu²⁺) con la periferia relajada por UFF restringido.
- **Caché** de las consultas a PubChem (no se repiten al reejecutar la app).
- **Filtro de cargas corregido**: ya no confunde el guion de enlace simple del SMILES
  con una carga formal.

## ▶️ Cómo ejecutar

### Windows (un clic)
Doble clic en **`Iniciar_LigandNexus.bat`**. La primera vez prepara un entorno local e
instala las dependencias (requiere **Python 3.10+** y conexión a internet); las
siguientes veces solo abre la app.

> ¿No tiene Python? Descárguelo de <https://www.python.org/downloads/> y marque
> *"Add Python to PATH"*.

#### ⚠️ Cómo detener la app (y por qué a veces «no abría»)

**Cerrar la pestaña del navegador NO detiene la aplicación.** El servidor sigue
corriendo en su ventana de consola, ocupando el puerto 8501. Para detenerlo hay que
**cerrar esa ventana negra** (o pulsar Ctrl+C dentro de ella).

Esto importa porque Streamlit, si encuentra el puerto ocupado, **no da error**: se
queda esperando en silencio a que se libere. Antes, volver a ejecutar el lanzador
dejaba la ventana aparentemente colgada y apilaba procesos en cada intento. Ahora el
`.bat` lo detecta al arrancar y ofrece tres opciones:

| Opción | Qué hace |
|--------|----------|
| **ENTER** | Abre en el navegador la instancia que ya está corriendo. |
| **R** | La reinicia — útil justo después de actualizar el código. |
| **S** | Sale sin tocar nada. |

Si algún día el puerto queda ocupado por un proceso huérfano, basta con cerrar la
ventana de consola de LigandNexus, o desde PowerShell:

```powershell
Get-NetTCPConnection -LocalPort 8501 -State Listen | ForEach-Object { Stop-Process -Id $_.OwningProcess -Force }
```

> **Nota para quien edite el `.bat`:** debe permanecer en **ASCII puro y con saltos
> CRLF**, y los paréntesis dentro de un `echo` han de ir escapados (`^(` `^)`). cmd.exe
> cuenta los paréntesis del texto dentro de un bloque, y un `(...)` seguido de un punto
> aborta el script entero con *«No se esperaba . en este momento»* — en todas las
> ejecuciones, aunque el entorno ya exista.

### Manual (cualquier sistema)
```bash
pip install -r requirements.txt
streamlit run app.py
```

### Ejecutable (.exe) y portable
Consulte **`packaging/README_empaquetado.md`** para generar un ejecutable de Windows o
un portable autocontenido.

## 📁 Estructura

```
LigandNexus/
├─ app.py                  Punto de entrada (navegación multipágina)
├─ ligandnexus/            Motor (lógica pura, reutilizable)
│  ├─ pubchem.py           Cliente PUG-REST (tasa + reintentos)
│  ├─ sources.py           Registro de bases de datos intercambiables
│  ├─ chem.py              RDKit: núcleos, donadores, sustituyentes
│  ├─ metales.py           Complejos metal–polifenol (Fe/Cu) + geometría
│  ├─ admet.py             Descriptores y reglas ADMET (Lipinski, logBB, QED)
│  ├─ screening.py         Filtros + selección curada
│  ├─ geometry.py          Geometrías 3D + inputs QM
│  └─ excel.py             Exportación a Excel
├─ views/                  Interfaz (una vista por modo)
├─ tests/                  Verificación empírica de motores e interfaz
├─ .streamlit/config.toml  Tema
├─ examples/               Salidas de muestra
├─ packaging/              Guías y scripts de empaquetado
├─ requirements.txt
├─ Iniciar_LigandNexus.bat Lanzador de un clic (Windows)
├─ LICENSE                 MIT
└─ README.md
```

## 📐 Notas técnicas

- Las **geometrías 3D** son una pre-optimización con campo de fuerzas (MMFF94/UFF),
  **no** una optimización DFT.
- La identificación de átomos coordinantes es una **heurística** química; para casos
  límite conviene el criterio del especialista.
- El **panel ADMET** es un tamiz de **priorización**, no una predicción
  farmacocinética: sus reglas son correlaciones estadísticas sobre fármacos orales
  conocidos. La `TPSA` que calcula es la de Ertl (solo N y O), que es la que espera
  la ecuación de Clark, y puede diferir de la que trae la base de datos.
- Los **complejos metal–polifenol** son **semillas geométricas** para optimizar después
  con DFT, no predicciones termodinámicas. En modo polinuclear los metales quedan solo
  bidentados: hay que completarles la esfera de coordinación antes del cálculo.
- Los datos provienen de **PubChem** (dominio público).

## 🔒 Seguridad y red

- **El servidor escucha solo en `127.0.0.1`.** La aplicación no tiene autenticación,
  así que no debe quedar expuesta a la red local. Está fijado en `.streamlit/config.toml`
  y en el lanzador; si necesita abrirla a otras máquinas, ponga usted el control de acceso.
- **La conexión con PubChem, ChEMBL y COCONUT verifica el certificado TLS.** No es un
  detalle formal: si alguien alterase las respuestas por el camino, los SMILES que
  devuelve la herramienta serían otros y **nada lo advertiría** — y esas estructuras
  acaban en cálculos DFT.
  Si su institución usa un proxy que rompe la cadena de certificados, puede
  desactivarla **a conciencia**, asumiendo que los datos podrían venir alterados:

  ```bat
  set LIGANDNEXUS_TLS_INSEGURO=1        :: Windows
  export LIGANDNEXUS_TLS_INSEGURO=1     # Linux / macOS
  ```

- **Los archivos que usted sube se tratan como no confiables.** El título de un `.sdf`
  y los identificadores y encabezados de un `.csv` son texto libre: se filtran antes de
  llegar al visor 3D y se neutralizan antes de exportarse a Excel, donde una celda que
  empiece por `=`, `+`, `-` o `@` se ejecutaría al abrir el libro en otra máquina.
- El tamaño máximo de subida es de **25 MB** (`maxUploadSize` en `.streamlit/config.toml`).

## 📄 Licencia

[MIT](LICENSE) © 2026 O. Molina

## ✍️ Cómo citar

> O. Molina. *LigandNexus: cribado virtual de derivados moleculares a partir de
> PubChem*, v0.2.5, 2026. <https://github.com/Newo-lab/ligandnexus>

El archivo [`CITATION.cff`](CITATION.cff) lleva estos mismos datos en el formato
que GitHub y Zenodo entienden.
