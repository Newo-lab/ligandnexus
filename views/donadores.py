"""
Vista: Identificador de átomos coordinantes
===========================================
Dada una molécula (CID, nombre o SMILES), señala qué heteroátomos pueden
coordinar un metal, los colorea por fuerza donadora y los numera — la base para
numerar posiciones en un análisis posterior.
"""

import pandas as pd
import streamlit as st

from ligandnexus import chem
from . import common


def _cargar_smiles(entrada, modo):
    """Devuelve (smiles, nombre, error)."""
    entrada = (entrada or "").strip()
    if not entrada:
        return None, None, "Escriba una entrada."
    if modo == "SMILES":
        m = chem.mol_from_smiles(entrada)
        if m is None:
            return None, None, "SMILES no válido."
        return entrada, "(SMILES introducido)", None
    info, err = common.resolver_entrada(entrada, modo)
    if err:
        return None, None, err
    nombre = info.get("Title") or info.get("IUPACName") or f"CID {info['CID']}"
    return info["SMILES"], nombre, None


def render():
    st.header("🎯 Identificador de átomos coordinantes")
    st.write(
        "Señala los heteroátomos (N, O, S, P) que pueden **coordinar un metal**, "
        "los colorea según su **fuerza donadora** y los **numera**. Útil para "
        "entender el modo de coordinación y para numerar posiciones antes de un "
        "análisis de sustituyentes.")

    modo = st.radio("Entrada", ["CID de PubChem", "Nombre", "SMILES"],
                    horizontal=True, key="don_modo")
    modo_key = {"CID de PubChem": "CID", "Nombre": "Nombre", "SMILES": "SMILES"}[modo]

    c1, c2 = st.columns([3, 1])
    with c1:
        ph = {"CID": "Ej.: 1474", "Nombre": "Ej.: 1,10-phenanthroline",
              "SMILES": "Ej.: c1ccnc(c1)-c1ccccn1"}[modo_key]
        entrada = st.text_input("Molécula", key="don_entrada", placeholder=ph,
                                label_visibility="collapsed")
    with c2:
        analizar = st.button("🔍 Analizar", type="primary", width="stretch")

    if analizar:
        smi, nombre, err = _cargar_smiles(entrada, modo_key)
        if err:
            st.error(err)
        else:
            st.session_state["don_smiles"] = smi
            st.session_state["don_nombre"] = nombre

    smi = st.session_state.get("don_smiles")
    if not smi:
        st.info("Introduzca una molécula y pulse **Analizar**.")
        return

    img, donadores = chem.imagen_donadores(smi, size=(520, 400), dark=common.es_oscuro())
    info = chem.identificar_donadores(smi, solo_relevantes=False)

    izq, der = st.columns([3, 2])
    with izq:
        st.markdown(f"**{st.session_state.get('don_nombre','')}**")
        if img:
            st.image(img, width="stretch")
        st.caption("Verde = donador fuerte · ámbar = moderado · rojo = débil. "
                   "Los números indican el orden de coordinación sugerido.")
    with der:
        st.markdown("#### Átomos coordinantes")
        if donadores:
            tabla = pd.DataFrame([
                {"N.º": i + 1, "Átomo": d["elemento"], "Fuerza": d["nivel"],
                 "Tipo": d["etiqueta"]}
                for i, d in enumerate(donadores)])
            st.dataframe(tabla, hide_index=True, width="stretch")
        else:
            st.warning("No se detectaron átomos coordinantes claros.")

        res = info["resumen"]
        m1, m2, m3 = st.columns(3)
        m1.metric("Fuertes", res.get("fuerte", 0))
        m2.metric("Moderados", res.get("moderado", 0))
        m3.metric("Débiles", res.get("débil", 0))

        # Sugerencia de denticidad (sitios de coordinación, no átomos sueltos:
        # un carboxilato cuenta como un solo sitio).
        sitios = chem.estimar_sitios(info["mol"], info["donadores"]) if info["mol"] else 0
        st.markdown("#### Lectura sugerida")
        etiqueta = {0: "", 1: "**Monodentado**", 2: "**Bidentado**",
                    3: "**Tridentado**", 4: "**Tetradentado**"}.get(sitios, "**Polidentado**")
        if sitios == 0:
            st.write("Sin donadores útiles: poco probable como ligando quelante.")
        elif sitios == 1:
            st.write(f"{etiqueta} probable (un solo punto de anclaje).")
        elif sitios == 2:
            st.write(f"{etiqueta} probable — candidato a quelato de 5–6 miembros si "
                     "los dos sitios están bien situados.")
        elif sitios == 3:
            st.write(f"{etiqueta} probable (tipo terpiridina / pinza).")
        else:
            st.write(f"{etiqueta} ({sitios} sitios): posible quelato envolvente.")
        st.caption(f"≈ {sitios} sitio(s) de coordinación "
                   "(los dos O de un carboxilato cuentan como uno).")

    with st.expander("ℹ️ ¿Cómo se decide la fuerza donadora?"):
        st.markdown(
            "- **N piridínico** (aromático, sin H): par libre disponible → *fuerte*.\n"
            "- **N amínico** (sp³): *fuerte*.\n"
            "- **N imínico** (C=N): *moderado*.\n"
            "- **N pirrólico** (aromático con H) y **N amídico** (junto a C=O): el par "
            "libre está deslocalizado → *débil* (solo coordina si se desprotona).\n"
            "- **O** carboxilato/alcóxido (aniónico): *fuerte*; carbonilo e hidroxilo: "
            "*moderado*; éter: *débil*.\n"
            "- **S** tiolato/tiol: *fuerte*; tioéter: *moderado*; S oxidado: descartado.\n\n"
            "Es una **heurística** basada en hibridación, aromaticidad y entorno; para "
            "casos límite conviene confirmarla con el criterio del químico.")
