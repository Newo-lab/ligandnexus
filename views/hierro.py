"""
Vista: Complejos de metales con polifenoles — Fe / Cu
=====================================================
Coordinación de Fe(II)/Fe(III) y Cu(I)/Cu(II) por polifenoles (catecoles,
galoílos, salicilatos e hidroxi-cetonas), los motivos quelantes de los taninos
hidrolizables y de los flavonoides.
Detecta los sitios de quelación reales (bidentados O,O), estima la especiación
mono/bis/tris frente al pH, y construye una geometría 3D del complejo —
**mononuclear** (un metal, varios ligandos) o **polinuclear** (un tanino con
varios metales)— como semilla para DFT.
"""

import pandas as pd
import streamlit as st

from ligandnexus import metales, geometry, chem
from . import common


# Polifenoles de ejemplo (CID de PubChem · SMILES de respaldo).
POLIFENOLES = {
    "Ácido gálico": (370, "OC(=O)c1cc(O)c(O)c(O)c1"),
    "Catecol": (289, "Oc1ccccc1O"),
    "Pirogalol": (1057, "Oc1cccc(O)c1O"),
    "Ácido protocatéquico": (72, "OC(=O)c1ccc(O)c(O)c1"),
    # SMILES canónico de PubChem CID 5281855. OJO: el respaldo anterior era un
    # regioisómero (misma fórmula y los mismos 2 catecoles, pero los OH de un
    # anillo corridos una posición, lo que rompía la simetría C2 de la molécula
    # real y desplazaba 0,41 Å el Fe···Fe del complejo polinuclear).
    "Ácido elágico": (5281855, "O=c1oc2c(O)c(O)cc3c(=O)oc4c(O)c(O)cc1c4c23"),
    "Ácido tánico (galotanino)": (16129778, None),
    # Con los dos centros (2R,3R) definidos, como en PubChem CID 65064: sin ellos
    # la geometría 3D salía de un diastereómero cualquiera.
    "EGCG (té verde)": (65064,
                        "O=C(O[C@@H]1Cc2c(O)cc(O)cc2O[C@@H]1c1cc(O)c(O)c(O)c1)"
                        "c1cc(O)c(O)c(O)c1"),
    "Quercetina": (5280343, "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"),
}

METAL_LABELS = {"Fe3+": "Fe(III)", "Fe2+": "Fe(II)", "Cu2+": "Cu(II)", "Cu+": "Cu(I)"}


def _resolver(entrada, modo):
    entrada = (entrada or "").strip()
    if not entrada:
        return None, None, "Escriba una entrada."
    if modo == "SMILES":
        if chem.mol_from_smiles(entrada) is None:
            return None, None, "SMILES no válido."
        return entrada, "(SMILES introducido)", None
    info, err = common.resolver_entrada(entrada, "CID" if modo == "CID" else "Nombre")
    if err:
        return None, None, err
    nombre = info.get("Title") or info.get("IUPACName") or f"CID {info['CID']}"
    return info["SMILES"], nombre, None


def _cargar_polifenol(nombre):
    """(smiles, nombre_a_mostrar, error) para un polifenol curado, por su CID.
    Usa el SMILES de respaldo (instantáneo) cuando existe; si no, resuelve el CID."""
    cid, smi = POLIFENOLES[nombre]
    if smi:
        return smi, f"{nombre} (CID {cid})", None
    info, err = common.resolver_entrada(str(cid), "CID")
    if err:
        return None, None, err
    return info["SMILES"], info.get("Title") or f"{nombre} (CID {cid})", None


def render():
    st.header("🩸 Complejos de metales con polifenoles")
    st.write(
        "Coordinación de **Fe(II)/Fe(III)** y **Cu(I)/Cu(II)** por polifenoles: "
        "**taninos hidrolizables** (ácido gálico, ácido tánico) y **flavonoides** "
        "(quercetina, luteolina). Un par de **-OH orto (catecol)**, el triplete de un "
        "**galoílo** o el **hidroxilo vecino a un carbonilo** cuentan como **un solo "
        "sitio bidentado O,O** — que es como estos metales realmente se unen.")
    st.caption("La esfera de coordinación se construye con **geometría cristalográfica** "
               "(M–O por metal y mordida O–M–O derivada del tamaño del anillo quelato, "
               "más el alargamiento Jahn-Teller en Cu²⁺) y la periferia orgánica se relaja "
               "con UFF restringido — son **semillas** para DFT, no predicciones "
               "termodinámicas.")

    # ── ⚡ Prueba rápida (lo primero: polifenol por CID + metal + pH) ──────────
    with st.container(border=True):
        st.markdown("#### ⚡ Prueba rápida")
        st.caption("Elige un **polifenol** (por su CID) y el **metal** en su estado de "
                   "oxidación; se carga y analiza al instante.")
        q1, q2, q3 = st.columns([2.3, 2.4, 1.3])
        with q1:
            pol = st.selectbox("Polifenol (por CID)", list(POLIFENOLES),
                               format_func=lambda n: f"{n} — CID {POLIFENOLES[n][0]}",
                               key="fe_quickpol")
        with q2:
            metal = st.radio("Metal (estado de oxidación)",
                             ["Fe3+", "Fe2+", "Cu2+", "Cu+"], horizontal=True,
                             key="fe_metal", format_func=lambda x: METAL_LABELS[x])
        with q3:
            pH = st.slider("pH", 1.0, 12.0, 7.0, 0.5, key="fe_pH")

        # Carga automática al cambiar el polifenol (instantánea con el SMILES de
        # respaldo; por CID si no lo hay). No pisa una búsqueda manual posterior.
        if st.session_state.get("_fe_lastpol") != pol:
            smi_new, nom_new, err = _cargar_polifenol(pol)
            if err:
                st.error(f"No se pudo cargar {pol}: {err}")
            else:
                st.session_state.update(fe_smiles=smi_new, fe_nombre=nom_new,
                                        _fe_lastpol=pol)
                for k in ("fe_mol", "fe_cinfo", "fe_lig3d"):
                    st.session_state.pop(k, None)

        with st.expander("🔎 …o buscar otra molécula (CID, nombre o SMILES)"):
            modo = st.radio("Tipo de entrada", ["CID de PubChem", "Nombre", "SMILES"],
                            horizontal=True, key="fe_modo")
            modo_key = {"CID de PubChem": "CID", "Nombre": "Nombre",
                        "SMILES": "SMILES"}[modo]
            mc1, mc2 = st.columns([3, 1])
            with mc1:
                ph_txt = {"CID": "Ej.: 370", "Nombre": "Ej.: tannic acid",
                          "SMILES": "Ej.: OC(=O)c1cc(O)c(O)c(O)c1"}[modo_key]
                entrada = st.text_input("Molécula", key="fe_entrada", placeholder=ph_txt,
                                        label_visibility="collapsed")
            with mc2:
                analizar = st.button("🔍 Analizar", type="primary", width="stretch")
            if analizar:
                smi_m, nom_m, err = _resolver(entrada, modo_key)
                if err:
                    st.error(err)
                else:
                    st.session_state.update(fe_smiles=smi_m, fe_nombre=nom_m,
                                            _fe_lastpol=pol)   # evita la recarga auto
                    for k in ("fe_mol", "fe_cinfo", "fe_lig3d"):
                        st.session_state.pop(k, None)

    # El complejo armado pertenece a UN metal: si se cambia el metal hay que
    # descartarlo, o la ficha de resultados sigue describiendo el complejo anterior.
    if st.session_state.get("_fe_lastmetal") != metal:
        st.session_state["_fe_lastmetal"] = metal
        for _k in ("fe_mol", "fe_cinfo"):
            st.session_state.pop(_k, None)

    smi = st.session_state.get("fe_smiles")
    if not smi:
        st.info("Elija un polifenol arriba para empezar.")
        return

    info = metales.analizar_polifenol(smi)
    sitios = info["sitios"]
    esp = metales.especiacion(metal, len(sitios), pH)

    # ── Lectura de especiación ────────────────────────────────────────────────
    st.markdown("**Especie dominante estimada**")
    if sitios:
        st.markdown(f"### {esp['especie']}-complejo · color {esp['color']}  ·  "
                    f"geometría {esp['geom']}")
    else:
        st.warning("No se detectaron sitios tipo catecol/galoílo/salicilato.")

    # ── 2D + tabla de sitios ──────────────────────────────────────────────────
    izq, der = st.columns([3, 2])
    with izq:
        st.markdown(f"**{st.session_state.get('fe_nombre','')}**")
        img, _ = metales.imagen_sitios(smi, size=(540, 400), dark=common.es_oscuro())
        if img:
            st.image(img, width="stretch")
        st.caption("Cada color numera un **sitio bidentado** de quelación.")
        # Ver la molécula libre en 3D (sin metal).
        if st.button("👁️ Ver la molécula en 3D (sin metal)", key="fe_ver3d"):
            with st.status("Generando geometría 3D de la molécula…"):
                mol3, conf, met, _en, _n = geometry.generar_3d(smi, n_confs=1)
                st.session_state["fe_lig3d"] = (geometry.pdb_block(mol3, conf)
                                                if mol3 is not None else None)
        if st.session_state.get("fe_lig3d"):
            common.mostrar_3d(st.session_state["fe_lig3d"], height=340, estilo="stick")
    with der:
        st.markdown("#### Sitios de quelación")
        if sitios:
            st.dataframe(pd.DataFrame([
                {"N.º": i + 1, "Tipo": s["tipo"], "Fuerza": s["fuerza"]}
                for i, s in enumerate(sitios)]), hide_index=True, width="stretch")
        # La cuarta columna va mas ancha: con cuatro iguales, la etiqueta
        # «Hidroxi-cetonas» sale cortada con puntos suspensivos.
        m1, m2, m3, m4 = st.columns([1, 1, 1, 1.6])
        m1.metric("Galoílos", info["n_galoil"])
        m2.metric("Catecoles", info["n_catecol"])
        m3.metric("Salicilatos", info["n_salicil"])
        m4.metric("Hidroxi-cetonas", info.get("n_hidroxicetona", 0),
                  help="Motivo de los flavonoides: un -OH vecino a un C=O. "
                       "«-5» es el 3-OH/4-C=O de los flavonoles (quelato de cinco "
                       "miembros, fuerte); «-6» es el 5-OH/4-C=O (seis miembros, "
                       "moderado porque ese -OH está trabado en un puente de "
                       "hidrógeno intramolecular).")
        st.caption(esp["nota"])

    if not sitios:
        st.stop()

    # ── Construcción 3D del complejo ──────────────────────────────────────────
    st.divider()
    st.subheader("🧬 Complejo 3D (semilla para DFT)")
    n_sitios = len(sitios)
    modo3d = st.radio(
        "Modo de complejo",
        ["Mononuclear (1 metal, varios ligandos)",
         "Polinuclear (1 molécula, varios metales — para taninos)"],
        key="fe_modo3d",
        help="Mononuclear: un centro metálico con mono/bis/tris ligandos idénticos. "
             "Polinuclear: una sola molécula polidentada con un metal en cada sitio "
             "(el caso real de los taninos, que enlazan varios metales).")
    es_poli = modo3d.startswith("Poli")
    inf_m = metales.metal_info(metal)

    if not es_poli:
        g1, g2, g3 = st.columns(3)
        idx_sitio = g1.selectbox("Sitio coordinante", range(n_sitios),
                                 format_func=lambda i: f"#{i+1} · {sitios[i]['tipo']}",
                                 key="fe_sitio")
        n_max = min(inf_m["max_bi"], 3)
        n_def = max(1, min(esp["n_ligandos"], n_max))
        n_lig = g2.select_slider("Nº de ligandos", options=list(range(1, n_max + 1)),
                                 value=n_def, key="fe_nlig",
                                 help=f"{METAL_LABELS[metal]} admite hasta {n_max} "
                                      f"ligandos bidentados ({inf_m['geom']}).")
        mult = g3.number_input("Multiplicidad de espín", 1, 30, inf_m["mult"],
                               key="fe_mult",
                               help="Fe(III) d⁵→6 · Fe(II) d⁶→5 · Cu(II) d⁹→2 · Cu(I) d¹⁰→1.")
    else:
        st.caption("Elija qué sitios ocupar con un metal. Por defecto, todos los "
                   "sitios fuertes (catecol/galoílo).")
        opciones = list(range(n_sitios))
        default = [i for i in opciones if sitios[i]["fuerza"] == metales.FUERTE] or opciones
        sel_sitios = st.multiselect(
            "Sitios a metalar", opciones, default=default, key="fe_polisitios",
            format_func=lambda i: f"#{i+1} · {sitios[i]['tipo']}")

    if st.button("🔩 Construir complejo 3D", type="primary", key="fe_build"):
        with st.status("Ensamblando el complejo…", expanded=False) as s:
            if es_poli:
                mol, cinfo = metales.construir_polinuclear(smi, sel_sitios, metal)
            else:
                mol, cinfo = metales.construir_complejo(
                    smi, sitios[idx_sitio]["o_idx"], int(n_lig), metal)
            if mol is None:
                s.update(label=f"No se pudo: {cinfo}", state="error")
                st.session_state.pop("fe_mol", None)
            else:
                st.session_state["fe_mol"] = mol
                st.session_state["fe_cinfo"] = cinfo
                st.session_state["fe_mult_user"] = (
                    cinfo["multiplicidad"] if es_poli else int(mult))
                s.update(label=f"Complejo listo · {mol.GetNumAtoms()} átomos",
                         state="complete")

    if st.session_state.get("fe_mol") is not None:
        mol = st.session_state["fe_mol"]
        cinfo = st.session_state["fe_cinfo"]
        # La multiplicidad se lee EN VIVO del control mientras esté a la vista: si se
        # tomara del momento en que se armó el complejo, cambiarla dejaría el .gjf con
        # el espín viejo y el cálculo DFT saldría equivocado sin ningún aviso.
        if cinfo["modo"] == "polinuclear":
            mult_u = int(cinfo["multiplicidad"])
        elif not es_poli:
            mult_u = int(mult)
        else:
            mult_u = int(st.session_state.get("fe_mult_user", inf_m["mult"]))

        a, b, c, d = st.columns(4)
        a.metric("Átomos", mol.GetNumAtoms())
        b.metric("Carga total", f"{cinfo['carga']:+d}")
        c.metric("Nº de metales", cinfo["n_metales"])
        feo = cinfo["fe_o"]
        d.metric("M–O ligando (Å)", f"{min(feo):.2f}–{max(feo):.2f}" if feo else "—")
        extra = [f"Multiplicidad (alto espín): **{mult_u}**"]
        bites = cinfo.get("bite")
        if bites:
            uniq = sorted(set(bites))
            extra.append("Mordida O–M–O: **" +
                         "/".join(f"{x:.0f}°" for x in uniq) + "**")
        eq_w, ax_w = cinfo.get("mo_agua_eq") or [], cinfo.get("mo_agua_ax") or []
        if eq_w or ax_w:
            partes = []
            if eq_w:
                partes.append(f"{len(eq_w)} ecuatorial(es) a "
                              f"{min(eq_w):.2f}–{max(eq_w):.2f} Å")
            if ax_w:
                partes.append(f"{len(ax_w)} axial(es) a {min(ax_w):.2f}–{max(ax_w):.2f} Å, "
                              "alargadas por Jahn-Teller")
            extra.append(f"{cinfo.get('n_aguas', 0)} H₂O de saturación: **"
                         + " · ".join(partes) + "**")
        if cinfo.get("mm"):
            extra.append(f"Distancias M–M: {', '.join(f'{x:.1f}' for x in cinfo['mm'])} Å")
        if cinfo.get("relajado"):
            extra.append("periferia orgánica relajada (UFF restringido)")
        st.caption(" · ".join(extra))
        if cinfo["choques"]:
            st.warning(f"El ensamblaje dejó {cinfo['choques']} contacto(s) corto(s); "
                       "la optimización DFT posterior los resolverá.")
        if cinfo.get("saturado") is False:
            st.info("**Esfera de coordinación incompleta.** En modo polinuclear cada "
                    "metal queda solo **bidentado** (2-coordinado): los sitios de un "
                    "tanino apuntan en direcciones arbitrarias y unas aguas puestas a "
                    "ciegas chocarían con la propia molécula. Antes de optimizar con "
                    "DFT conviene completar cada centro con agua o disolvente "
                    "explícito. El modo mononuclear sí satura la coordinación.")

        pdb = metales.complejo_pdb(mol)
        estilo = st.radio("Estilo", ["stick", "ball"], horizontal=True, key="fe_estilo")
        common.mostrar_3d(pdb, height=460, estilo=estilo)
        st.caption(f"{METAL_LABELS[cinfo['metal']]} · {cinfo['modo']} · "
                   f"{cinfo['n_metales']} centro(s) metálico(s).")

        # ── Descargas ─────────────────────────────────────────────────────────
        nombre_base = (st.session_state.get("fe_nombre", "complejo")
                       .replace(" ", "_").replace("(", "").replace(")", ""))
        etq_m = cinfo["metal"].replace("+", "").replace("3", "3").replace("2", "2")
        nombre_base = f"{etq_m}_{nombre_base}_{cinfo['n_metales']}M"
        comentario = (f"{cinfo['metal']} {cinfo['modo']} {cinfo['n_metales']}M "
                      f"carga={cinfo['carga']} mult={mult_u}")
        xyz = metales.complejo_xyz(mol, comentario)
        xyzl = metales.complejo_xyz_lineas(mol)
        gjf = geometry.construir_input_gaussian(
            nombre_base, cinfo["carga"], mult_u, xyzl,
            metodo="B3LYP", base="def2-SVP", keywords="opt",
            titulo=f"{comentario} — semilla LigandNexus")

        st.markdown("#### Descargas")
        st.caption("La esfera de coordinación se fija a datos cristalográficos (M–O del "
                   "metal y mordida O–M–O derivada del tamaño del anillo quelato) y el "
                   "orgánico se relaja con UFF congelando el "
                   "núcleo del quelato; sigue siendo una **semilla**: optimícela con DFT "
                   "(para 3d, B3LYP/TPSSh + def2 + dispersión es un punto de partida "
                   "razonable). En Cu²⁺ se añaden 2 aguas axiales (Jahn-Teller); en los "
                   "polinucleares la multiplicidad es el límite de alto espín no "
                   "interactuante, ajústela según el acoplamiento esperado.")
        dc1, dc2 = st.columns(2)
        dc1.download_button("⬇️ Geometría .xyz", data=xyz,
                            file_name=f"{nombre_base}.xyz", mime="text/plain",
                            width="stretch")
        dc2.download_button("⬇️ Input Gaussian .gjf", data=gjf,
                            file_name=f"{nombre_base}.gjf", mime="text/plain",
                            width="stretch")

    # ── Química de fondo ──────────────────────────────────────────────────────
    with st.expander("ℹ️ La química metal–polifenol (por qué esto es coherente)"):
        st.markdown(
            "- **Catecol** y **galoílo** (1,2,3-trihidroxi, como en el ácido gálico y "
            "los galotaninos): al desprotonarse forman **catecolato**, quelato "
            "bidentado O,O dianiónico. El *tris*-catecolato de Fe(III) tiene log β ≳ 40.\n"
            "- **Fe(III)** (ácido duro): afinidad enorme por estos O aniónicos → hasta "
            "*tris*, octaédrico. **Fe(II)** coordina más débil y se oxida a Fe(III) al aire.\n"
            "- **Cu(II)** (d⁹, plano-cuadrado, Jahn-Teller): forma catecolatos, pero "
            "además **oxida el catecol a o-quinona** (redox tipo Fenton, ROS) — clave en "
            "el estrés oxidativo y la química del cobre en Alzheimer. **Cu(I)** es blando "
            "(prefiere N/S); afinidad baja por O duros.\n"
            "- **Flavonoides** (quercetina, luteolina, crisina): además del catecol del "
            "anillo B, quelan por el **hidroxilo vecino al carbonilo**. El **3-OH/4-C=O** "
            "cierra un anillo de cinco y es un sitio fuerte; el **5-OH/4-C=O** cierra uno "
            "de seis y es moderado, porque ese -OH está trabado en un puente de hidrógeno "
            "intramolecular. Los dos comparten el mismo C=O, así que **solo uno** puede "
            "alojar el metal: se asigna al de cinco, que es el preferido. A diferencia del "
            "catecolato, este sitio quela como **monoanión** (fenolato + C=O neutro).\n"
            "- **Taninos hidrolizables** (p. ej. ácido tánico): muchos galoílos → enlazan "
            "**varios metales a la vez** (modo **polinuclear**). Aquí cada sitio recibe un "
            "metal; los centros quedan a la distancia que impone la geometría de la molécula.\n"
            "- **Especiación con el pH** (Fe): ácido→*mono*, neutro→*bis*, básico→*tris* "
            "(viraje verde→azul→rojo del hierro-galato / tinta ferrogálica).\n"
            "- **Construcción geométrica:** cada complejo se arma sobre el politopo del "
            "metal (octaedro / octaedro tetragonal / tetraedro) con **distancias M–O "
            "cristalográficas** (Fe³⁺ 2.01 Å, Fe²⁺ 2.14 Å, Cu²⁺ 1.95 Å ecuatorial + "
            "2.35 Å axial por **Jahn-Teller** —que solo tiene Cu²⁺ d⁹: Fe³⁺ d⁵ y Cu⁺ d¹⁰ "
            "no, por capa semillena y llena—, Cu⁺ 2.10 Å). La **mordida O–M–O** no se "
            "tabula por metal: se deriva del O···O que impone el anillo quelato (≈2.61 Å "
            "con cinco miembros, ≈2.72 Å con seis) junto con la M–O, así que un "
            "salicilato abre más la mordida que un catecolato sobre el mismo metal. La "
            "coordinación se satura con agua y la periferia orgánica se relaja con **UFF "
            "restringido** (núcleo del quelato congelado). El metal se excluye del campo "
            "de fuerzas porque UFF no lo parametriza bien.\n\n"
            "El modelo desprotona **solo** los O que coordinan y ajusta carga y "
            "multiplicidad; es una **semilla geométrica**, no una predicción termodinámica.")
