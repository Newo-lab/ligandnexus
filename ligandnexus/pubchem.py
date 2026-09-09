"""
pubchem.py — Cliente robusto de PubChem PUG-REST
================================================
Acceso a PubChem con las cortesías que el servicio pide y que la versión 0.1 no
tenía: control de tasa de peticiones (máx. ~5/s), reintentos con espera creciente
ante «servidor ocupado» (HTTP 503/429) y manejo tolerante del nombre de la clave
del SMILES (PubChem migró CanonicalSMILES -> SMILES/ConnectivitySMILES en 2025).

Todas las funciones aceptan `log=print` para reportar progreso a la interfaz.

Ref.: https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest
"""

from __future__ import annotations

import json
import os
import time
import ssl
import threading
import urllib.request
import urllib.parse
import urllib.error

import pandas as pd

PUGREST = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

# Propiedades que la interfaz ofrece descargar (subconjunto de las de PubChem).
PROPIEDADES_DISPONIBLES = [
    "IUPACName", "MolecularFormula", "MolecularWeight",
    "HBondDonorCount", "HBondAcceptorCount", "XLogP", "TPSA",
    "SMILES",
]

# La verificacion TLS va ACTIVA por defecto. Si alguien alterase la respuesta
# de PubChem por el camino, los SMILES que devuelve esta herramienta serian
# otros y nada lo advertiria: el riesgo aqui no es de confidencialidad (no se
# envia ninguna credencial) sino de INTEGRIDAD de los datos cientificos.
#
# Quien este detras de un proxy corporativo que rompe la cadena de
# certificados puede desactivarla, a conciencia y asumiendo ese riesgo, con:
#     set LIGANDNEXUS_TLS_INSEGURO=1        (Windows)
#     export LIGANDNEXUS_TLS_INSEGURO=1     (Linux / macOS)
TLS_INSEGURO = os.environ.get("LIGANDNEXUS_TLS_INSEGURO") == "1"

_SSL_CTX = ssl.create_default_context()
if TLS_INSEGURO:
    _SSL_CTX.check_hostname = False
    _SSL_CTX.verify_mode = ssl.CERT_NONE

# ── Control de tasa ───────────────────────────────────────────────────────────
# PubChem pide no superar ~5 peticiones/segundo por IP. Mantenemos un intervalo
# mínimo entre llamadas, protegido por un lock para que sirva también cuando la
# interfaz dispara varias descargas.
_MIN_INTERVALO = 0.22          # s entre peticiones (~4.5/s)
_lock = threading.Lock()
_ultima_peticion = [0.0]


def _esperar_turno():
    with _lock:
        ahora = time.monotonic()
        espera = _MIN_INTERVALO - (ahora - _ultima_peticion[0])
        if espera > 0:
            time.sleep(espera)
        _ultima_peticion[0] = time.monotonic()


class PubChemError(RuntimeError):
    """Error de comunicación con PubChem tras agotar los reintentos."""


def _request(url, data=None, timeout=60, reintentos=4, log=None):
    """
    Petición HTTP con control de tasa y reintentos con backoff exponencial ante
    503/429 (servidor ocupado) y errores de red transitorios. Devuelve bytes.
    """
    headers = {"User-Agent": "LigandNexus/0.2 (research tool)"}
    if data is not None:
        headers["Content-Type"] = "application/x-www-form-urlencoded"

    ultimo_error = None
    for intento in range(reintentos):
        _esperar_turno()
        try:
            req = urllib.request.Request(url, data=data, headers=headers)
            with urllib.request.urlopen(req, timeout=timeout, context=_SSL_CTX) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            # 404 = no existe: es una respuesta legítima, no se reintenta.
            if e.code == 404:
                raise
            # 503/429 = ocupado o límite de tasa: reintentar con espera creciente.
            if e.code in (503, 429, 500):
                ultimo_error = e
                espera = 1.5 * (2 ** intento)
                if log:
                    log(f"    PubChem ocupado (HTTP {e.code}); reintento en {espera:.0f}s")
                time.sleep(espera)
                continue
            raise
        except (urllib.error.URLError, TimeoutError, ssl.SSLError) as e:
            ultimo_error = e
            espera = 1.5 * (2 ** intento)
            if log:
                log(f"    Red inestable ({type(e).__name__}); reintento en {espera:.0f}s")
            time.sleep(espera)
    raise PubChemError(f"PubChem no respondió tras {reintentos} intentos: {ultimo_error}")


def _get_json(url, **kw):
    return json.loads(_request(url, **kw))


# ── SMILES tolerante al nombre de la clave ────────────────────────────────────

def extraer_smiles(prop: dict) -> str:
    """SMILES de un dict de propiedades, tolerante al renombre de PubChem 2025."""
    for clave in ("SMILES", "ConnectivitySMILES", "CanonicalSMILES", "IsomericSMILES"):
        if prop.get(clave):
            return prop[clave]
    for k, v in prop.items():
        if "smiles" in k.lower() and v:
            return v
    return ""


# ── Identificación de la molécula ─────────────────────────────────────────────

def cid_a_smiles(cid, timeout=30) -> str:
    url = f"{PUGREST}/compound/cid/{int(cid)}/property/SMILES/JSON"
    data = _get_json(url, timeout=timeout)
    try:
        smi = extraer_smiles(data["PropertyTable"]["Properties"][0])
    except (KeyError, IndexError):
        smi = ""
    if not smi:
        raise PubChemError(f"No se obtuvo SMILES para el CID {cid}")
    return smi


def info_molecula(cid, timeout=30) -> dict:
    """Datos de verificación de un CID (nombre amigable, fórmula, MW, InChIKey…)."""
    props = "Title,IUPACName,MolecularFormula,MolecularWeight,InChIKey,SMILES"
    url = f"{PUGREST}/compound/cid/{int(cid)}/property/{props}/JSON"
    data = _get_json(url, timeout=timeout)
    p = data["PropertyTable"]["Properties"][0]
    return {
        "CID": cid,
        "Title": p.get("Title", ""),
        "IUPACName": p.get("IUPACName", ""),
        "MolecularFormula": p.get("MolecularFormula", ""),
        "MolecularWeight": p.get("MolecularWeight", ""),
        "InChIKey": p.get("InChIKey", ""),
        "SMILES": extraer_smiles(p),
    }


def buscar_cid_por_nombre(nombre, timeout=30) -> list:
    nombre = (nombre or "").strip()
    if not nombre:
        return []
    url = f"{PUGREST}/compound/name/{urllib.parse.quote(nombre)}/cids/JSON"
    try:
        data = _get_json(url, timeout=timeout)
        return data.get("IdentifierList", {}).get("CID", [])
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return []
        raise


def sinonimos(cid, n=6, timeout=30) -> list:
    url = f"{PUGREST}/compound/cid/{int(cid)}/synonyms/JSON"
    try:
        data = _get_json(url, timeout=timeout)
        return data["InformationList"]["Information"][0].get("Synonym", [])[:n]
    except Exception:
        return []


# ── Búsqueda de subestructura + descarga de propiedades ───────────────────────

def buscar_derivados(smiles_core, max_records=10000, poll_wait=4,
                     max_polls=120, log=print) -> list:
    """
    CIDs que contienen `smiles_core` como subestructura. Maneja la respuesta
    inmediata y la asíncrona (polling con ListKey). Devuelve lista de CIDs (ints).
    """
    url = f"{PUGREST}/compound/substructure/smiles/JSON"
    payload = urllib.parse.urlencode(
        {"smiles": smiles_core, "MaxRecords": max_records}).encode()

    log(f"Buscando subestructura: {smiles_core}  (MaxRecords={max_records:,})")
    try:
        data = json.loads(_request(url, data=payload, timeout=90, log=log))
    except urllib.error.HTTPError as e:
        cuerpo = e.read().decode(errors="replace")[:200]
        log(f"  HTTP {e.code}: {cuerpo}")
        return []
    except PubChemError as e:
        log(f"  {e}")
        return []

    if "IdentifierList" in data:
        cids = data["IdentifierList"]["CID"]
        log(f"  {len(cids):,} CIDs (respuesta inmediata)")
        return cids

    if "Waiting" in data:
        listkey = data["Waiting"]["ListKey"]
        poll_url = (f"{PUGREST}/compound/listkey/{listkey}/cids/JSON"
                    f"?MaxRecords={max_records}")
        log("  Búsqueda asíncrona, esperando resultados…")
        for intento in range(max_polls):
            time.sleep(poll_wait)
            try:
                data2 = _get_json(poll_url, timeout=30, log=log)
                if "IdentifierList" in data2:
                    cids = data2["IdentifierList"]["CID"]
                    log(f"  {len(cids):,} CIDs tras {intento + 1} sondeos")
                    return cids
            except urllib.error.HTTPError as e:
                if e.code == 202:      # aún procesando
                    continue
                log(f"  Error en sondeo: HTTP {e.code}")
                return []
        log("  Se agotó el tiempo de la búsqueda")
        return []

    log(f"  Respuesta inesperada: {list(data.keys())}")
    return []


def descargar_propiedades(cids, propiedades=None, batch_size=200,
                          log=print) -> pd.DataFrame:
    """
    Descarga las `propiedades` de cada CID por lotes. A diferencia de la v0.1,
    si un lote falla se registra en el log Y se cuenta, para no perder filas en
    silencio. Devuelve (DataFrame, n_fallidos_estimado) via atributo .attrs.
    """
    if propiedades is None:
        propiedades = list(PROPIEDADES_DISPONIBLES)
    props_str = ",".join(propiedades)

    filas, n_fallidos = [], 0
    total = len(cids)
    for i in range(0, total, batch_size):
        batch = cids[i:i + batch_size]
        cid_str = ",".join(map(str, batch))
        url = f"{PUGREST}/compound/cid/{cid_str}/property/{props_str}/JSON"
        try:
            data = _get_json(url, timeout=60, log=log)
            filas.extend(data["PropertyTable"]["Properties"])
        except Exception as e:
            n_fallidos += len(batch)
            log(f"  ⚠️ Lote {i // batch_size + 1} falló y se omitió ({len(batch)} CIDs): {e}")
        log(f"  Propiedades: {min(i + batch_size, total):,}/{total:,}")
    df = pd.DataFrame(filas)
    df.attrs["n_fallidos"] = n_fallidos
    if n_fallidos:
        log(f"  ⚠️ {n_fallidos} CIDs no se pudieron descargar (lotes fallidos).")
    return df
