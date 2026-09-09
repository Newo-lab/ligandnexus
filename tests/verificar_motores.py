# -*- coding: utf-8 -*-
"""
verificar_motores.py — Validacion empirica de LigandNexus
========================================================
Contrasta cada motor contra DATOS EXTERNOS, no contra si mismo:

  · chem.py      denticidades documentadas de 18 ligandos clasicos
                 (bipiridina, terpiridina, EDTA, ciclam, salen...)
  · screening.py filtros de carga y de metales sobre casos limite reales
  · metales.py   distancias M-O y angulos de mordida cristalograficos
                 (tris-catecolato de Fe(III)) y orden de radios de Shannon
  · geometry.py  longitudes y angulos de enlace experimentales
                 (benceno, etano, agua, metano) y formato de los inputs QM
  · admet.py     descriptores publicados de farmacos conocidos
  · metales.py   complejo Fe-galato del articulo de Zaccaron (2013) y las
                 distancias cristalograficas que cita
  · excel.py     integridad del libro exportado
  · views/hierro los polifenoles de ejemplo, contra el SMILES de su CID de
                 PubChem (el del acido elagico era un regioisomero)
  · views/       cableado de la interfaz conducido con AppTest y las
                 fuentes simuladas: invalidacion del estado de sesion,
                 exportaciones que deben coincidir con lo que se ve en
                 pantalla, y humo de las siete vistas

Ejecutar desde la raiz del proyecto:
    .venv/Scripts/python.exe tests/verificar_motores.py
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import io, numpy as np
from itertools import combinations
from rdkit import Chem
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
import ligandnexus as lf
from ligandnexus import chem, screening, metales, geometry, admet, excel

TOT=[0,0]
def chk(c, m):
    TOT[1]+=1; TOT[0]+=bool(c)
    if not c: print(f"      FALLA: {m}")

print(f"LigandNexus v{lf.__version__} — verificacion empirica de TODOS los motores")
print("="*72)

print("\n[1] chem.py — denticidad vs valores EXPERIMENTALES")
LIG=[("bipiridina","c1ccc(-c2ccccn2)nc1",2),("fenantrolina","c1cnc2c(c1)ccc1cccnc21",2),
("terpiridina","c1ccnc(-c2cccc(-c3ccccn3)n2)c1",3),("etilendiamina","NCCN",2),
("dietilentriamina","NCCNCCN",3),("glicina","NCC(=O)O",2),("IDA","OC(=O)CNCC(=O)O",3),
("NTA","OC(=O)CN(CC(=O)O)CC(=O)O",4),("EDTA","OC(=O)CN(CC(=O)O)CCN(CC(=O)O)CC(=O)O",6),
("8-HQ","Oc1cccc2cccnc12",2),("acac","CC(=O)CC(C)=O",2),("oxalico","OC(=O)C(=O)O",2),
("catecol","Oc1ccccc1O",2),("piridina","c1ccncc1",1),("imidazol","c1c[nH]cn1",1),
("ciclam","C1CNCCNCCNCCNC1",4),("TPMA","c1ccc(CN(Cc2ccccn2)Cc2ccccn2)nc1",4),
("salen","Oc1ccccc1C=NCCN=Cc1ccccc1O",4)]
for n,s,e in LIG:
    d=chem.identificar_donadores(s); chk(chem.estimar_sitios(d["mol"],d["donadores"])==e,f"{n}")
print(f"    {len(LIG)} ligandos con denticidad documentada")

print("\n[2] chem.py — nitrogenos SIN par libre (bug corregido)")
for n,s,e in [("nitrobenceno","c1ccccc1[N+](=O)[O-]",0),("4-nitropiridina","[O-][N+](=O)c1ccncc1",1),
 ("amonio cuaternario","C[N+](C)(C)C",0),("sulfonamida","CS(=O)(=O)N",0),
 ("bencenosulfonato","c1ccccc1S(=O)(=O)[O-]",0),("anilina","Nc1ccccc1",1)]:
    d=chem.identificar_donadores(s); chk(chem.estimar_sitios(d["mol"],d["donadores"])==e,n)
print("    6 casos")

print("\n[3] screening.py — filtro de cargas por carga NETA (bug corregido)")
for n,s,e in [("nitrobenceno","C1=CC=C(C=C1)[N+](=O)[O-]",False),("N-oxido","[O-][N+]1=CC=CC=C1",False),
 ("azida","[N-]=[N+]=NC1=CC=CC=C1",False),("betaina","C[N+](C)(C)CC(=O)[O-]",False),
 ("piridinio","C1=CC=[NH+]C=C1",True),("benzoato","[O-]C(=O)c1ccccc1",True),("benceno","c1ccccc1",False)]:
    chk(screening._tiene_carga(s)==e,n)
print("    7 casos")

print("\n[4] screening.py — filtro de metales con notacion de Hill (bug corregido)")
for f,e in [("C10H8N2",0),("C6H5BrN2",0),("C6H5ClN",0),("C10H8FeN2",1),("C6H5NaO",1),
 ("C6H5BO2",1),("C21H16N2OSi",1),("C10H8N2OSn",1),("C6H4NOSe",1),("C6H5NSb",1)]:
    chk(screening._tiene_metal(f)==bool(e),f)
print("    10 formulas")

print("\n[5] metales.py — parametros vs CRISTALOGRAFIA y radios de Shannon")
chk(abs(metales.METALES["Fe3+"]["mo"]-2.01)<0.03,"Fe3+-O")
chk(abs(metales.METALES["Fe2+"]["mo"]-2.13)<0.06,"Fe2+-O")
chk(abs(metales.METALES["Cu2+"]["mo"]-1.95)<0.05,"Cu2+-O eq")
chk(abs(metales.METALES["Cu2+"]["mo_ax"]-2.35)<0.15,"Cu2+-O ax")
chk(metales.METALES["Fe3+"]["mo"]<metales.METALES["Fe2+"]["mo"],"orden Shannon")
mol,inf=metales.construir_complejo("Oc1ccccc1O",(0,7),3,"Fe3+")
P=np.array([list(mol.GetConformer().GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
fe=inf["metal_idx"][0]; feo=[float(np.linalg.norm(P[o]-P[fe])) for o in inf["odon"]]
def ang(i,j,k):
    a,b=P[i]-P[j],P[k]-P[j]
    return float(np.degrees(np.arccos(np.clip(a@b/(np.linalg.norm(a)*np.linalg.norm(b)),-1,1))))
angs=sorted(ang(i,fe,j) for i,j in combinations(inf["odon"],2))
chk(abs(np.mean(feo)-2.01)<0.02,"Fe-O armado"); chk(np.std(feo)<0.01,"dispersion Fe-O")
chk(len([a for a in angs if a<135])==12,"12 cis"); chk(len([a for a in angs if a>=135])==3,"3 trans")
chk(inf["carga"]==-3 and inf["multiplicidad"]==6,"carga/mult tris-catecolato")
chk(inf["choques"]==0,"sin choques")
print("    11 comprobaciones (tris-catecolato de Fe(III))")

print("\n[6] metales.py — motivos quelantes de flavonoides")
for n,s,e in [("quercetina","O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12",2),
 ("luteolina","O=c1cc(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12",2),("crisina","O=c1cc(-c2ccccc2)oc2cc(O)cc(O)c12",1),
 ("ac. galico","OC(=O)c1cc(O)c(O)c(O)c1",1),("aspirina","CC(=O)Oc1ccccc1C(=O)O",0),
 ("resorcinol","Oc1cccc(O)c1",0)]:
    chk(len(metales.analizar_polifenol(s)["sitios"])==e,n)
print("    6 casos")

print("\n[7] geometry.py — geometrias vs valores EXPERIMENTALES")
def geo(s):
    m,c,_,_,_=geometry.generar_3d(s,n_confs=5)
    return m,np.array([list(m.GetConformer(c).GetAtomPosition(i)) for i in range(m.GetNumAtoms())])
def D(P,i,j): return float(np.linalg.norm(P[i]-P[j]))
def A(P,i,j,k):
    a,b=P[i]-P[j],P[k]-P[j]
    return float(np.degrees(np.arccos(np.clip(a@b/(np.linalg.norm(a)*np.linalg.norm(b)),-1,1))))
m,P=geo("c1ccccc1")
cc=[D(P,b.GetBeginAtomIdx(),b.GetEndAtomIdx()) for b in m.GetBonds()
    if m.GetAtomWithIdx(b.GetBeginAtomIdx()).GetAtomicNum()==6 and m.GetAtomWithIdx(b.GetEndAtomIdx()).GetAtomicNum()==6]
chk(abs(np.mean(cc)-1.390)<0.02,"benceno C-C")
for s,e,t in (("CC",1.535,0.03),("C#C",1.203,0.03),("C=C",1.339,0.03)):
    m,P=geo(s); chk(abs(D(P,0,1)-e)<t,s)
m,P=geo("O"); h=[a.GetIdx() for a in m.GetAtoms() if a.GetAtomicNum()==1]
chk(abs(D(P,0,h[0])-0.958)<0.03,"agua O-H"); chk(abs(A(P,h[0],0,h[1])-104.5)<4,"agua H-O-H")
m,P=geo("C"); h=[a.GetIdx() for a in m.GetAtoms() if a.GetAtomicNum()==1]
chk(abs(A(P,h[0],0,h[1])-109.47)<2,"metano")
print("    7 medidas")

print("\n[8] geometry.py — formato de los inputs QM")
lin=["O 0.0 0.0 0.0","H 0.958 0.0 0.0"]
g=geometry.construir_input_gaussian("x",-3,6,lin); o=geometry.construir_input_orca("x",-3,6,lin)
chk(g.splitlines()[0].startswith("%chk="),"gaussian chk"); chk("-3 6" in g,"gaussian carga/mult")
chk(any(l.startswith("#p ") for l in g.splitlines()),"gaussian ruta")
chk("* xyz -3 6" in o,"orca xyz"); chk(o.rstrip().endswith("*"),"orca cierre")
print("    5 comprobaciones")

print("\n[9] admet.py — vs FARMACOS con valores publicados")
for n,s,mw,tp,lip in [("aspirina","CC(=O)Oc1ccccc1C(=O)O",180.16,63.6,True),
 ("ibuprofeno","CC(C)Cc1ccc(C(C)C(=O)O)cc1",206.28,37.3,True),
 ("paracetamol","CC(=O)Nc1ccc(O)cc1",151.16,49.3,True),
 ("quercetina","O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12",302.24,131.4,True),
 ("ciclosporina","CC[C@H]1NC(=O)[C@H]([C@H](O)[C@H](C)C/C=C/C)N(C)C(=O)[C@H](C(C)C)N(C)C(=O)[C@H](CC(C)C)N(C)C(=O)[C@H](CC(C)C)N(C)C(=O)[C@@H](C)NC(=O)[C@H](C)NC(=O)[C@H](CC(C)C)N(C)C(=O)[C@@H](NC(=O)[C@H](CC(C)C)N(C)C(=O)CN(C)C1=O)C(C)C",1202.61,278.8,False)]:
    d=admet.evaluar(s)
    chk(abs(d["MW"]-mw)<0.6,f"{n} MW"); chk(abs(d["TPSA"]-tp)<1.5,f"{n} TPSA")
    chk(d["Lipinski"]==lip,f"{n} Lipinski")
# extremos de logBB bien ordenados
lb=lambda s: admet.evaluar(s)["logBB"]
chk(lb("CC(C)Cc1ccc(C(C)C(=O)O)cc1") > lb("O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"),"orden logBB")
chk(admet.evaluar("O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12")["PAINS"] is True,"PAINS quercetina")
print("    17 comprobaciones")

print("\n[10] excel.py — exportacion integra")
import pandas as pd
df=pd.DataFrame({"CID":[2244,2519],"SMILES":["CC(=O)Oc1ccccc1C(=O)O","Cn1c(=O)c2c(ncn2C)n(C)c1=O"]})
buf=io.BytesIO(); excel.exportar_excel(df,"SMILES",buf,con_imagenes=False,log=lambda *a: None)
b=buf.getvalue(); chk(len(b)>3000,"excel generado")
chk(b[:2]==b"PK","excel es un zip valido (xlsx)")
import openpyxl
wb=openpyxl.load_workbook(io.BytesIO(b)); ws=wb.active
chk(ws.max_row==3,"3 filas (cabecera + 2)")
print("    3 comprobaciones")

print("\n[11] metales.py — arreglos del modo de metales (sacado de beta)")
# Jahn-Teller SOLO donde existe: Cu(II) d9. Fe(III) d5 y Cu(I) d10 no lo tienen.
for met, esp in (("Fe3+", False), ("Fe2+", False), ("Cu+", False), ("Cu2+", True)):
    _m, i = metales.construir_complejo("Oc1ccccc1O", (0, 7), 1, met)
    chk(i["jahn_teller"] == esp and bool(i["mo_agua_ax"]) == esp, f"jahn-teller {met}")
    chk(i["choques"] == 0, f"sin choques {met}")
    chk("mo_axial" not in i, "clave enganosa mo_axial retirada")
# Un mismo sitio no puede alojar dos metales, y la entrada mala se rechaza limpia.
for arg in (None, [], [99], ["x"]):
    _m, _i = metales.construir_polinuclear("Oc1ccccc1O", arg, "Fe3+")
    chk(_m is None, f"polinuclear rechaza {arg!r}")
_m, i = metales.construir_polinuclear("Oc1ccccc1O", [0, 0], "Fe3+")
chk(_m is not None and i["n_metales"] == 1 and i["carga"] == 1, "polinuclear deduplica sitios")
chk(i.get("saturado") is False, "polinuclear declara su insaturacion")
# La carga la fija la quimica del sitio: catecolato dianionico, hidroxi-cetona monoanionica.
Q = "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"
for sitio in metales.analizar_polifenol(Q)["sitios"]:
    _m, i = metales.construir_complejo(Q, sitio["o_idx"], 1, "Fe3+")
    esperada = 3 + (-2 if sitio["tipo"] in ("catecol", "galoil", "salicilato") else -1)
    chk(i["carga"] == esperada, f"carga del sitio {sitio['tipo']}")
# Codigo muerto retirado.
chk(not hasattr(metales, "especiacion_fe"), "especiacion_fe retirada")
chk(all("bite" not in v for v in metales.METALES.values()), "parametro bite muerto retirado")
# La mordida se DERIVA del tamano del anillo quelato: la de 6 abre mas que la de 5.
_m, i5 = metales.construir_complejo("Oc1ccccc1O", (0, 7), 1, "Fe3+")
CRI = "O=c1cc(-c2ccccc2)oc2cc(O)cc(O)c12"
sc = metales.analizar_polifenol(CRI)["sitios"][0]
_m, i6 = metales.construir_complejo(CRI, sc["o_idx"], 1, "Fe3+")
chk(max(i6["bite"]) > max(i5["bite"]), "el anillo de 6 abre mas la mordida que el de 5")
print("    22 comprobaciones")

print("\n[12] views — cableado de la interfaz (auditoria del 4-sep-2026)")
# Los bloques de arriba prueban los MOTORES. Estos fallos vivian en el CABLEADO
# de la interfaz —estado de sesion que no se invalidaba, botones que se
# borraban solos— y por eso pasaron desapercibidos entre 113 comprobaciones de
# quimica correctas. Se conducen las vistas reales con AppTest y las fuentes
# simuladas (sin red). Detalle en _trabajo_2026-09-04/INFORME_bugs_exportar_admet.md
from streamlit.testing.v1 import AppTest as _AT

_DIR = os.path.dirname(os.path.abspath(__file__))
_RAIZ = os.path.dirname(_DIR)
_ARN_COO = os.path.join(_DIR, "_arnes_coordinacion.py")
_ARN_GEN = os.path.join(_DIR, "_arnes_general.py")


def _ss(at, k):
    return at.session_state[k] if k in at.session_state else None


def _sel(at, key):
    return next((s for s in at.selectbox if s.key == key), None)


def _sel_lbl(at, frag):
    return next((s for s in at.selectbox if frag.lower() in (s.label or "").lower()), None)


def _btn(at, frag):
    return next((b for b in at.button if frag in b.label), None)


def _dl(at, frag, excluir=None):
    """Botones de descarga cuyo rotulo contiene `frag`.

    No se usa `at.download_button`: ese atajo solo existe en algunas versiones
    de Streamlit (1.58 lo retiro y la suite reventaba con AttributeError). El
    nodo sigue en el arbol con type == 'download_button' y con su rotulo, asi
    que se recorre a mano cuando el atajo no esta disponible.
    """
    def sirve(x):
        lab = getattr(x, "label", "") or ""
        return frag in lab and (excluir is None or excluir not in lab)

    atajo = getattr(at, "download_button", None)
    if atajo is not None:
        return [x for x in atajo if sirve(x)]

    hallados = []

    def recorrer(nodo):
        if getattr(nodo, "type", None) == "download_button" and sirve(nodo):
            hallados.append(nodo)
        hijos = getattr(nodo, "children", None)
        if isinstance(hijos, dict):
            hijos = list(hijos.values())
        for h in (hijos or []):
            recorrer(h)

    recorrer(at.main)
    return hallados


def _coo_con_admet():
    """Vista de coordinacion conducida hasta tener el perfil ADMET calculado."""
    at = _AT.from_file(_ARN_COO, default_timeout=180)
    at.run()
    at.text_input(key="coo_entrada").set_value("1318").run()
    at.button[0].click().run()                       # Verificar
    at.button(key="coo_buscar_der").click().run()    # Buscar derivados
    _se_et = _sel(at, "coo_admet_etapa")
    if _se_et is not None:
        _se_et.set_value("Todos los derivados (①)").run()
    at.button(key="coo_admet_coo_df_calc").click().run()
    return at


# (a) El perfil no debe sobrevivir a un cambio de molecula. Antes sobrevivia y
#     Exportar lo ofrecia POR DEFECTO —a veces como unica opcion— con los
#     derivados de la molecula anterior, sin ningun aviso.
at = _coo_con_admet()
_adm = _ss(at, "coo_admet")
chk(_adm is not None and len(_adm) == 4, "el perfil ADMET cubre los 4 derivados")
_se = _sel_lbl(at, "exportar")
chk(_se is not None and any("ADMET" in o for o in _se.options),
    "Exportar ofrece el perfil recien calculado")
at.text_input(key="coo_entrada").set_value("1474").run()
at.button[0].click().run()                           # cambia de molecula
chk(_ss(at, "coo_admet") is None, "cambiar de molecula borra el perfil ADMET")
at.button(key="coo_buscar_der").click().run()
_se = _sel_lbl(at, "exportar")
chk(_se is not None and not any("ADMET" in o for o in _se.options),
    "Exportar ya no ofrece el perfil de la molecula anterior")

# (b) Un filtro que no deja ninguna molecula es un resultado legitimo: lo que se
#     exporta tiene que coincidir con lo que dice la pantalla.
at = _coo_con_admet()
_ms = next((m for m in at.multiselect if m.key == "coo_admet_coo_df_bhe"), None)
if _ms is not None:
    _ms.set_value(["atraviesa"]).run()               # ninguna lo cumple
_adm = _ss(at, "coo_admet")
chk(_adm is not None and len(_adm) == 0, "un filtro que deja 0 se refleja en la exportacion")
_ms = next((m for m in at.multiselect if m.key == "coo_admet_coo_df_bhe"), None)
if _ms is not None:
    _ms.set_value(["no atraviesa"]).run()            # solo el dicarboxilico
_adm = _ss(at, "coo_admet")
chk(_adm is not None and len(_adm) == 1, "un filtro que deja 1 tambien se refleja")

# (c) Cambiar de etapa sin recalcular no puede dejar en pie el perfil de la otra.
at = _coo_con_admet()
_n = next((x for x in at.number_input if "peso" in (x.label or "").lower()), None)
if _n is not None:
    _n.set_value(50.0).run()
_ba = _btn(at, "Aplicar filtros")
if _ba is not None:
    _ba.click().run()
_cf = _ss(at, "coo_filt")
chk(_cf is not None and len(_cf) == 3, "el tope de peso descarta el derivado pesado")
_se_et = _sel(at, "coo_admet_etapa")
if _se_et is not None:
    _se_et.set_value("Lista depurada (②)").run()
chk(_ss(at, "coo_admet") is None, "cambiar de etapa no deja un perfil rancio")

# (d) El boton de descargar el Excel vivia DENTRO del `if st.button(...)`, asi que
#     desaparecia en el primer rerun, incluido el que provoca su propio clic.
at = _coo_con_admet()
_bg = _btn(at, "Generar Excel")
if _bg is not None:
    _bg.click().run()
chk(bool(_dl(at, "Descargar Excel")),
    "tras generar aparece el boton de descargar Excel")
# Cualquier re-ejecucion vale: el fallo era que el boton, al vivir dentro del
# `if st.button(...)`, desaparecia en el primer rerun (incluido el de su clic).
at.run()
chk(bool(_dl(at, "Descargar Excel")),
    "el boton sobrevive a una re-ejecucion posterior")
_cb = next((c for c in at.checkbox if "estructuras 2D" in c.label), None)
if _cb is not None:
    _cb.set_value(True).run()
chk(not _dl(at, "Descargar Excel"),
    "al cambiar las opciones el Excel viejo deja de ofrecerse")
chk(any("vuelva a generarlo" in c.value for c in at.caption),
    "y se explica por que, en vez de desaparecer en silencio")

# (e) El CSV se arma al pulsar (data invocable), no en cada rerun.
_dc = _dl(at, "CSV", excluir="ADMET")
chk(bool(_dc), "la pestana Exportar ofrece el CSV")
at.run()
chk(not list(at.exception), "el CSV diferido no lanza excepcion al re-ejecutar")

# (f) Acotar un filtro no puede tirar un calculo de minutos: el perfil se rehace
#     solo desde la memoria de admet.py (medido: 7,5 ms/molecula sin memoria).
atg = _AT.from_file(_ARN_GEN, default_timeout=180)
atg.run()
atg.text_input(key="gen_entrada").set_value("1318").run()
atg.button[0].click().run()
atg.button[1].click().run()
_bc = next((b for b in atg.button if b.key == "gen_admet_calc"), None)
if _bc is not None:
    _bc.click().run()
chk("_gen_admet_df" in atg.session_state, "perfil ADMET calculado en cribado general")
_sl = next((s for s in atg.slider if s.key == "gen_flt_MolecularWeight"), None)
if _sl is not None:
    _sl.set_range(_sl.value[0], 200.0).run()
chk("_gen_admet_df" in atg.session_state, "mover un filtro NO tira el perfil")
chk("_gen_admet_df" in atg.session_state
    and len(atg.session_state["_gen_admet_df"]) == 2,
    "y se rehace sobre las 2 que quedan")

# (g) Perfilar dos veces el mismo conjunto duplicaba la columna TPSA_admet, y
#     entonces row.get() devolvia una Serie que openpyxl rechazaba.
_d2 = pd.DataFrame({"SMILES": ["c1cnc2c(c1)ccc1cccnc21", "Cn1cnc2c1c(=O)n(C)c(=O)n2C"],
                    "TPSA": [25.8, 58.4]})
_a2 = admet.evaluar_dataframe(admet.evaluar_dataframe(_d2, "SMILES"), "SMILES")
chk(len(set(_a2.columns)) == len(_a2.columns), "perfilar dos veces no duplica columnas")
chk(all(c in _a2.columns for c in ("TPSA", "TPSA_admet", "TPSA_admet2")),
    "la TPSA de la base y las dos calculadas conviven con nombres distintos")
_bx = io.BytesIO()
try:
    excel.exportar_excel(_a2, "SMILES", _bx, con_imagenes=False, log=lambda *a: None)
    chk(_bx.getvalue()[:2] == b"PK", "y el Excel del doble perfil sale valido")
except Exception as _e:
    # Con columnas duplicadas, row.get() devuelve una Serie y openpyxl la rechaza.
    chk(False, f"y el Excel del doble perfil sale valido ({type(_e).__name__})")

# (h) La memoria de perfiles: lo ya calculado no se vuelve a calcular.
chk(admet.faltan_por_calcular(["c1cnc2c(c1)ccc1cccnc21"]) == 0, "el perfil queda memorizado")
chk(admet.faltan_por_calcular(["C1CCCCC1CCO"]) == 1, "una molecula nueva no lo esta")

# (i) Humo: las siete vistas tienen que renderizar sin excepcion.
_PLANTILLA = ('import sys\nsys.path.insert(0, r"{raiz}")\n'
              'import streamlit as st\nst.set_page_config(layout="wide")\n'
              'from views import common, {v}\ncommon.aplicar_tema()\n{v}.render()\n')
for _v in ("inicio", "general", "coordinacion", "donadores", "hierro", "geometria", "acerca"):
    _p = os.path.join(_DIR, f"_humo_{_v}.py")
    with open(_p, "w", encoding="utf-8") as _f:
        _f.write(_PLANTILLA.format(raiz=_RAIZ, v=_v))
    try:
        _a = _AT.from_file(_p, default_timeout=120)
        _a.run()
        chk(not list(_a.exception), f"la vista {_v} renderiza")
    except Exception as _e:
        chk(False, f"la vista {_v} renderiza ({type(_e).__name__})")
    finally:
        try:
            os.remove(_p)
        except OSError:
            pass

# (j) Geometrias 3D tiene que ofrecer los conjuntos YA FILTRADOS por ADMET: quien
#     acote por Lipinski o por barrera quiere llevarse ESOS candidatos a 3D.
#     Y nunca un conjunto vacio, que no sirve para generar geometrias.
_ARN_GEO = os.path.join(_DIR, "_arnes_geometria.py")
_dfa = admet.evaluar_dataframe(
    pd.DataFrame({"CID": [1318, 2519], "SMILES": ["c1cnc2c(c1)ccc1cccnc21",
                                                  "Cn1cnc2c1c(=O)n(C)c(=O)n2C"]}),
    "SMILES")

ag = _AT.from_file(_ARN_GEO, default_timeout=120)
ag.session_state["coo_admet"] = _dfa
ag.session_state["coo_admet_de"] = "Candidatos seleccionados"
ag.session_state["gen_perfil"] = _dfa
ag.run()
_r = next((r for r in ag.radio if r.key == "geo_origen"), None)
if _r is not None:
    _r.set_value("Un conjunto del cribado").run()
_sg = next((s for s in ag.selectbox if s.key == "geo_fuente"), None)
chk(_sg is not None and any("ADMET" in o for o in _sg.options),
    "Geometrias 3D ofrece el conjunto perfilado de coordinacion")
chk(_sg is not None and sum("ADMET" in o for o in _sg.options) == 2,
    "y tambien el del cribado general")
chk(_sg is not None and any("Candidatos seleccionados" in o for o in _sg.options),
    "el rotulo dice de que etapa sale el perfil")

ag = _AT.from_file(_ARN_GEO, default_timeout=120)
ag.session_state["coo_admet"] = _dfa.head(0)          # filtro que no dejo nada
ag.session_state["coo_sel"] = _dfa
ag.run()
_r = next((r for r in ag.radio if r.key == "geo_origen"), None)
if _r is not None:
    _r.set_value("Un conjunto del cribado").run()
_sg = next((s for s in ag.selectbox if s.key == "geo_fuente"), None)
chk(_sg is not None and not any("ADMET" in o for o in _sg.options),
    "un conjunto vacio no se ofrece para generar geometrias")

print("    33 comprobaciones")

print("\n[13] metales.py — vs ZACCARON 2013 (acido galico + hierro)")
# Fuente externa: S. Zaccaron, R. Ganzerla, M. Bortoluzzi, «Iron complexes with
# gallic acid: a computational study on coordination compounds of interest for
# the preservation of cultural heritage», J. Coord. Chem. 2013, 66, 1709-1719.
# DOI 10.1080/00958972.2013.790019. Las distancias son cristalograficas y vienen
# de su ref. [2] (Wunderlich, Weber, Bergerhoff, Z. Anorg. Allg. Chem. 598/599
# (1991) 371), el unico complejo de tinta ferrogalica caracterizado por rayos X.
GALICO = "OC(=O)c1cc(O)c(O)c(O)c1"        # H4L: acido 3,4,5-trihidroxibenzoico

_ag = metales.analizar_polifenol(GALICO)
# El articulo describe la quelacion POR LOS FENOLATOS ("three phenate-type
# oxygens of every ligand chelate two metal ions"); el carboxilato PUENTEA dos
# hierros, no quela uno. Y los OH estan en 3,4,5, o sea meta y para respecto al
# COOH: no hay motivo salicilato aunque convivan un OH y un COOH en la molecula.
chk(_ag["n_galoil"] == 1, "acido galico: 1 sitio galoilo")
chk(_ag["n_salicil"] == 0, "acido galico: NINGUN salicilato (los OH no son orto al COOH)")
chk(len(_ag["sitios"]) == 1 and _ag["sitios"][0]["tipo"] == "galoíl",
    "y el unico sitio es el galoilo")

_sitio = _ag["sitios"][0]["o_idx"]

# «slightly distorted octahedral geometry with SIX Fe-O bonds»: la esfera se
# completa con agua, asi que el numero de coordinacion es 6 con 1, 2 o 3 galatos.
for _n in (1, 2, 3):
    _m, _i = metales.construir_complejo(GALICO, _sitio, _n, "Fe3+")
    _cn = len(_i["fe_o"]) + len(_i["mo_agua_eq"]) + len(_i["mo_agua_ax"])
    chk(_cn == 6, f"Fe(III)-galato con {_n} ligando(s): numero de coordinacion 6")
    chk(_i["choques"] == 0, f"Fe(III)-galato con {_n} ligando(s): sin choques")

# Distancias Fe-O del cristal: 2.006 (carboxilato), 2.000 (fenolatos 3 y 5) y
# 2.028 A (fenolato 4). Media 2.011 A.
_m3, _i3 = metales.construir_complejo(GALICO, _sitio, 3, "Fe3+")
_CRIST = (2.000, 2.006, 2.028)
chk(all(min(_CRIST) - 0.03 <= d <= max(_CRIST) + 0.03 for d in _i3["fe_o"]),
    "todas las Fe-O caen en el rango cristalografico +-0.03 A")
chk(abs(sum(_i3["fe_o"]) / len(_i3["fe_o"]) - sum(_CRIST) / 3) < 0.02,
    "la Fe-O media coincide con la del cristal (<0.02 A)")

# Carga: Fe(3+) con tres galatos dianionicos -> -3.
chk(_i3["carga"] == -3, "tris-galato de Fe(III): carga -3")

# Multiplicidades. El articulo fija sus referencias: [Fe(EDTA)]- (Fe(III)) en
# estado SEXTETO y [Fe(EDTA)]2- (Fe(II)) en QUINTETO.
chk(metales.METALES["Fe3+"]["mult"] == 6, "Fe(III) alto espin: sexteto, como [Fe(EDTA)]-")
chk(metales.METALES["Fe2+"]["mult"] == 5, "Fe(II) alto espin: quinteto, como [Fe(EDTA)]2-")
chk(_i3["multiplicidad"] == 6, "el complejo de Fe(III) sale en sexteto")

# Invariante fisico: en una configuracion d^n caben como mucho n desapareados si
# n<=5, y 10-n si n>5. Cualquier multiplicidad por encima de eso es imposible.
# (El propio articulo se traiciona aqui: en dos sitios asigna «sextet» a
# [Fe(H2O)6]2+, que es d6 y admite 4 desapareados como maximo -> quinteto. Su
# pareja de EDTA, en la misma frase, si esta bien.)
for _k, _v in metales.METALES.items():
    _d = _v["d"]
    _max_desap = _d if _d <= 5 else 10 - _d
    chk(_v["mult"] <= _max_desap + 1, f"{_k}: multiplicidad posible para d{_d}")

# El articulo llama a la geometria «slightly distorted octahedral» por el empaque
# cristalino, no por Jahn-Teller: ni d5 ni d6 alto espin lo tienen.
for _met in ("Fe3+", "Fe2+"):
    _m, _i = metales.construir_complejo(GALICO, _sitio, 1, _met)
    chk(not _i["jahn_teller"], f"{_met}: sin Jahn-Teller (la distorsion del cristal es otra cosa)")

# Controles del motivo, para que «1 galoilo y 0 salicilatos» no sea casualidad.
chk(metales.analizar_polifenol("Oc1cccc(O)c1O")["n_galoil"] == 1,
    "control: pirogalol (el galoilo sin COOH) tambien da 1 galoilo")
chk(len(metales.analizar_polifenol("OC(=O)c1ccccc1")["sitios"]) == 0,
    "control negativo: el acido benzoico solo (COOH sin OH) no da ningun sitio")
chk(metales.analizar_polifenol("OC(=O)c1ccccc1O")["n_salicil"] == 1,
    "control positivo: el acido salicilico (OH ORTO al COOH) si da salicilato")

# Semilla de calculo con la carga y la multiplicidad del articulo.
_gjf = geometry.construir_input_gaussian(
    "fe_galato", _i3["carga"], _i3["multiplicidad"],
    metales.complejo_xyz_lineas(_m3), metodo="M06", base="SDD")
chk("\n-3 6\n" in _gjf, "el input de Gaussian lleva la carga y el sexteto correctos")
chk(_gjf.count("Fe") >= 1, "y contiene el hierro")

print("    24 comprobaciones")

print("\n[14] views/hierro.py — los polifenoles de ejemplo SON los de su CID")
# El respaldo del ACIDO ELAGICO era un REGIOISOMERO: misma formula C14H6O8 y los
# mismos dos catecoles, pero con los OH de un anillo corridos una posicion. No lo
# delataba ni la formula ni el recuento de sitios; si lo delata la SIMETRIA (el
# acido elagico real tiene sus dos mitades equivalentes) y, aguas abajo, un
# Fe...Fe desplazado 0,41 A en el complejo polinuclear. El del EGCG venia sin
# estereoquimica. Ambos verificados contra PubChem en vivo el 4-sep-2026; aqui se
# fijan sus formas canonicas para que la suite lo guarde SIN red.
from rdkit.Chem import rdMolDescriptors as _rdMD
from rdkit.Chem import inchi as _inchi
from views.hierro import POLIFENOLES as _POLI

# Ademas del SMILES canonico y la formula se fija el InChIKey, que es el
# guardian mas fuerte: un hash de la conectividad COMPLETA, estereoquimica
# incluida. Dos moleculas que solo se diferencian en donde esta un -OH dan
# claves distintas, que es justo lo que se le escapo al respaldo del elagico.
# Los tres valores se verificaron contra PubChem en vivo el 4-sep-2026.
_ESPERADO = {
    370:     ("O=C(O)c1cc(O)c(O)c(O)c1", "C7H6O5",
              "LNTHITQWFMADLM-UHFFFAOYSA-N"),
    289:     ("Oc1ccccc1O", "C6H6O2",
              "YCIMNLLNPGFGHC-UHFFFAOYSA-N"),
    1057:    ("Oc1cccc(O)c1O", "C6H6O3",
              "WQGWDDDVZFFDIG-UHFFFAOYSA-N"),
    72:      ("O=C(O)c1ccc(O)c(O)c1", "C7H6O4",
              "YQUVCSBJEUQKSH-UHFFFAOYSA-N"),
    5281855: ("O=c1oc2c(O)c(O)cc3c(=O)oc4c(O)c(O)cc1c4c23", "C14H6O8",
              "AFSDNFLWKVMVRB-UHFFFAOYSA-N"),
    65064:   ("O=C(O[C@@H]1Cc2c(O)cc(O)cc2O[C@@H]1c1cc(O)c(O)c(O)c1)"
              "c1cc(O)c(O)c(O)c1", "C22H18O11",
              # El segundo bloque, WIYYLYMNSA, codifica la estereoquimica: sin
              # los centros (2R,3R) definidos la clave saldria distinta.
              "WMBWREPUVVBILR-WIYYLYMNSA-N"),
    5280343: ("O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12", "C15H10O7",
              "REFJWTPEDVJJIY-UHFFFAOYSA-N"),
}
for _nom, (_cid, _smi) in _POLI.items():
    if _smi is None:
        continue                      # se resuelve por CID: no hay nada que fijar
    _m = chem.mol_from_smiles(_smi)
    chk(_m is not None, f"{_nom}: SMILES valido")
    if _m is None:
        continue
    _can, _frm, _key = _ESPERADO[_cid]
    chk(Chem.MolToSmiles(_m) == _can, f"{_nom}: es la molecula del CID {_cid}")
    chk(_rdMD.CalcMolFormula(_m) == _frm, f"{_nom}: formula {_frm}")
    chk(_inchi.MolToInchiKey(_m) == _key, f"{_nom}: InChIKey del CID {_cid}")

# El invariante que habria cazado el regioisomero sin necesidad de red: el acido
# elagico es SIMETRICO (22 atomos pesados en 11 clases de equivalencia). El
# impostor tenia 22 clases, o sea ninguna simetria.
_el = chem.mol_from_smiles(_POLI["Ácido elágico"][1])
_clases = len(set(Chem.CanonicalRankAtoms(_el, breakTies=False)))
chk(_el.GetNumAtoms() == 22 and _clases == 11,
    "acido elagico: simetrico, sus dos mitades son equivalentes")
# El EGCG es (2R,3R): sin los centros definidos la geometria 3D sale de un
# diastereomero cualquiera.
_eg = chem.mol_from_smiles(_POLI["EGCG (té verde)"][1])
chk(len(Chem.FindMolChiralCenters(_eg, includeUnassigned=True)) == 2,
    "EGCG: dos centros estereogenicos")
chk(all(a for _, a in Chem.FindMolChiralCenters(_eg, includeUnassigned=True)),
    "EGCG: y los dos estan ASIGNADOS")

print("    31 comprobaciones")

print("\n[15] metales.py — Fe(III) + acido elagico de punta a punta")
ELAGICO = _POLI["Ácido elágico"][1]
_ae = metales.analizar_polifenol(ELAGICO)
# Dilactona del acido hexahidroxidifenico: un catecol en cada anillo aromatico.
chk(_ae["n_catecol"] == 2, "acido elagico: 2 catecoles")
chk(_ae["n_galoil"] == 0 and _ae["n_salicil"] == 0 and _ae["n_hidroxicetona"] == 0,
    "y ningun otro motivo (sus lactonas no son salicilato ni hidroxi-cetona)")
chk(all(s["tipo"] == "catecol" for s in _ae["sitios"]), "los dos sitios son catecoles")

# MONONUCLEAR: un Fe(III) con 1, 2 o 3 elagatos; la esfera se completa con agua.
_s0 = _ae["sitios"][0]["o_idx"]
for _n, _q in ((1, 1), (2, -1), (3, -3)):
    _m, _i = metales.construir_complejo(ELAGICO, _s0, _n, "Fe3+")
    chk(_m is not None, f"Fe(III)+elagico mononuclear, {_n} ligando(s): se construye")
    _cn = len(_i["fe_o"]) + len(_i["mo_agua_eq"]) + len(_i["mo_agua_ax"])
    chk(_cn == 6, f"Fe(III)+elagico, {_n} ligando(s): numero de coordinacion 6")
    chk(_i["choques"] == 0, f"Fe(III)+elagico, {_n} ligando(s): sin choques")
    chk(_i["carga"] == _q, f"Fe(III)+elagico, {_n} ligando(s): carga {_q}")
    chk(_i["multiplicidad"] == 6, f"Fe(III)+elagico, {_n} ligando(s): sexteto")

# POLINUCLEAR: un Fe(III) en cada catecol, que es el caso real de este ligando.
_mp, _ip = metales.construir_polinuclear(ELAGICO, [0, 1], "Fe3+")
chk(_mp is not None, "Fe(III)+elagico polinuclear: se construye")
chk(_ip["n_metales"] == 2, "polinuclear: dos centros de Fe(III), uno por catecol")
chk(_ip["choques"] == 0, "polinuclear: sin choques")
chk(_ip["carga"] == 2, "polinuclear: carga +2 (2 Fe3+ y 2 catecolatos dianionicos)")
chk(_ip["multiplicidad"] == 11, "polinuclear: 2 x 5 desapareados + 1 = 11")
chk(all(abs(d - 2.01) < 0.05 for d in _ip["fe_o"]), "polinuclear: Fe-O impuesta a 2,01 A")
# La separacion Fe...Fe es la diagonal de la molecula: el regioisomero daba 10,95.
chk(len(_ip["mm"]) == 1 and 11.0 < _ip["mm"][0] < 11.8,
    "polinuclear: Fe...Fe en los 11,4 A que impone el elagico REAL")

# Semilla de calculo cuantico del complejo polinuclear.
_xyz = metales.complejo_xyz_lineas(_mp)
chk(len(_xyz) == _mp.GetNumAtoms(), "el .xyz del complejo tiene todos los atomos")
chk(sum(1 for l in _xyz if l.split()[0] == "Fe") == 2, "y los dos hierros")
_g = geometry.construir_input_gaussian("fe_elagico", _ip["carga"], _ip["multiplicidad"],
                                       _xyz, metodo="M06", base="SDD")
chk("\n2 11\n" in _g, "el input de Gaussian lleva carga +2 y multiplicidad 11")

print("    28 comprobaciones")

print("\n[16] seguridad — los vectores de la auditoria del 8-sep-2026")
# Estos NO son fallos de quimica: son de SUPERFICIE DE ATAQUE, y se colaron
# porque nadie miraba la entrada no confiable (el archivo que sube el usuario).
# Cada comprobacion REPRODUCE el ataque y exige que el parche lo neutralice; si
# alguien revierte un parche, aqui falla. Detalle en _trabajo_2026-09-08/.
import ssl as _ssl, tempfile as _tmp, zipfile as _zf, importlib as _il
import pandas as _pd
from rdkit.Chem import AllChem as _AC
from views import common as _cm

# ── H-01. XSS almacenado: el titulo de un SDF lo escribe quien hizo el archivo,
# RDKit lo guarda en _Name, MolToPDBBlock lo copia al registro COMPND y py3Dmol
# lo interpola SIN escapar en el HTML del visor.
_CARGA = "</script><script>alert(1)</script>"
_SDF = _CARGA + """
     RDKit          2D

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0
    1.5000    0.0000    0.0000 C   0  0
  1  2  1  0
M  END
$$$$
"""
_m = next(iter(Chem.ForwardSDMolSupplier(io.BytesIO(_SDF.encode()))))
chk(_m is not None and _m.GetProp("_Name") == _CARGA,
    "el SDF hostil se parsea y RDKit guarda el titulo en _Name")
_m3 = Chem.AddHs(_m); _AC.EmbedMolecule(_m3, randomSeed=42)
_pdb = Chem.MolToPDBBlock(_m3)
chk(_CARGA in _pdb, "el ataque es valido: la carga llega al PDB por COMPND")

_lim = _cm._limpiar_pdb(_pdb)
chk(_CARGA not in _lim, "_limpiar_pdb elimina la carga util")
chk("<" not in _lim and ">" not in _lim, "no sobrevive ningun < ni > en el bloque")
chk("COMPND" not in _lim, "el registro COMPND se descarta entero")
_nat = sum(1 for l in _lim.splitlines() if l.startswith(("ATOM", "HETATM")))
chk(_nat == _m3.GetNumAtoms(), "se conservan los atomos: el visor sigue sirviendo")
import py3Dmol as _p3
_v = _p3.view(width="100%", height=200); _v.addModel(_lim, "pdb")
chk(_CARGA not in _v._make_html(), "la carga no llega al HTML de py3Dmol")

# No basta con que _limpiar_pdb funcione: hay que probar que mostrar_3d LO LLAMA.
# Al reventar esa unica linea, las comprobaciones de arriba seguian en verde. Se
# intercepta py3Dmol y st.iframe para ver que llega de verdad al visor.
_cap = {}


class _VistaFalsa:
    def addModel(self, m, fmt=None): _cap["modelo"] = m
    def setStyle(self, *a, **k): pass
    def setBackgroundColor(self, *a, **k): pass
    def zoomTo(self): pass
    def _make_html(self): return "<div></div>"


_view_real, _iframe_real = _p3.view, _cm.st.iframe
try:
    _p3.view = lambda *a, **k: _VistaFalsa()
    _cm.st.iframe = lambda *a, **k: None
    _cm.mostrar_3d(_pdb)
finally:
    _p3.view, _cm.st.iframe = _view_real, _iframe_real
chk(_CARGA not in _cap.get("modelo", ""),
    "mostrar_3d SI llama al saneador (existir no basta: hay que estar cableado)")

# ── H-02. CWE-1236: Excel ejecuta como formula toda celda que empiece por
# = + - @. Los IDs y los NOMBRES DE COLUMNA salen del archivo que sube el
# usuario, y el .xlsx se le envia a un tercero, que es quien lo ejecutaria.
chk(_cm is not None and excel._celda_segura("=cmd|calc") == "'=cmd|calc",
    "_celda_segura neutraliza el = inicial")
for _v_ok in ("CCO", "etanol", 46.07, 1318, None):
    chk(excel._celda_segura(_v_ok) == _v_ok,
        f"_celda_segura NO toca un valor legitimo ({_v_ok!r})")

_df = _pd.DataFrame([{"ID": "=cmd|'/c calc.exe'!A1", "SMILES": "CCO", "CID": "1318"}])
_df = _df.rename(columns={"CID": "CID"})
_df["@SUM(1+1)"] = "-2+3+cmd"          # encabezado Y valor hostiles
_ruta = os.path.join(_tmp.mkdtemp(), "reg.xlsx")
excel.exportar_excel(_df, "SMILES", _ruta, con_imagenes=False, log=lambda *a: None)
with _zf.ZipFile(_ruta) as _z:
    _hoja = _z.read("xl/worksheets/sheet1.xml").decode("utf8")
    _rels = [n for n in _z.namelist() if "sheet1.xml.rels" in n]
    _link = bool(_rels) and b"pubchem" in _z.read(_rels[0])
chk("<f>" not in _hoja, "ninguna celda del libro exportado quedo como formula")
chk(_link, "el hipervinculo legitimo a PubChem sigue funcionando")

# ── H-03 y TLS. Configuracion de red: la app no tiene autenticacion, asi que
# el servidor no puede escuchar en 0.0.0.0 (el defecto de Streamlit es None).
_CFG = open(os.path.join(_RAIZ, ".streamlit", "config.toml"), encoding="utf-8").read()
chk('address = "127.0.0.1"' in _CFG, "config.toml ata el servidor a 127.0.0.1")
chk("maxUploadSize" in _CFG, "config.toml limita el tamano de subida")
_BAT = open(os.path.join(_RAIZ, "Iniciar_LigandNexus.bat"), "rb").read()
chk(b"--server.address 127.0.0.1" in _BAT, "el lanzador .bat tambien fija la direccion")
chk(all(b <= 127 for b in _BAT), "el .bat sigue siendo ASCII puro")
chk(b"\n" not in _BAT.replace(b"\r\n", b""), "el .bat sigue siendo CRLF, sin LF sueltos")

from ligandnexus import pubchem as _pc
chk(_pc.TLS_INSEGURO is False, "la verificacion TLS viene ACTIVA por defecto")
chk(_pc._SSL_CTX.verify_mode == _ssl.CERT_REQUIRED, "el contexto exige certificado")
chk(_pc._SSL_CTX.check_hostname is True, "y comprueba el nombre del servidor")
os.environ["LIGANDNEXUS_TLS_INSEGURO"] = "1"
try:
    _pc2 = _il.reload(_pc)
    chk(_pc2._SSL_CTX.verify_mode == _ssl.CERT_NONE,
        "la variable de entorno si permite desactivarla a conciencia")
finally:
    del os.environ["LIGANDNEXUS_TLS_INSEGURO"]
    _il.reload(_pc)

# ── H-04 y fugas. El piso de la dependencia transitiva y el canario de datos.
_REQ = open(os.path.join(_RAIZ, "requirements.txt"), encoding="utf-8").read()
chk("gitpython>=3.1.59" in _REQ, "requirements fija el piso de GitPython (29 avisos)")
chk("requests" not in _REQ, "no se declara requests, que no se usa")

import re as _re
_FUGAS = _re.compile(r"universidad|nacional de colombia|distrital|@gmail|owenmoli", _re.I)
_sucios = []
for _dp, _dn, _fn in os.walk(_RAIZ):
    if any(p in _dp for p in (".git", ".venv", "__pycache__", "_trabajo", "_backup")):
        continue
    for _f in _fn:
        # Este mismo archivo se excluye: lleva los terminos DENTRO del patron
        # de busqueda, asi que se denunciaria a si mismo.
        if _f == os.path.basename(__file__):
            continue
        if _f.endswith((".py", ".md", ".txt", ".toml", ".cff", ".bat")):
            try:
                if _FUGAS.search(open(os.path.join(_dp, _f), encoding="utf-8").read()):
                    _sucios.append(_f)
            except Exception:
                pass
chk(not _sucios, f"ningun archivo publicable nombra universidad ni correo {_sucios}")

print("    28 comprobaciones")

print("\n"+"="*72)
print(f"RESULTADO: {TOT[0]}/{TOT[1]} comprobaciones empiricas correctas")
print("TODOS LOS MOTORES VERIFICADOS" if TOT[0]==TOT[1] else ">>> HAY FALLOS <<<")
