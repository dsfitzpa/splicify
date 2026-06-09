"""
Priority-based host inference for plasmid tokenization.

Rules (highest priority first):
  1. Payload/module signals: lentiviral, AAV, baculovirus, T-DNA, yeast CEN/ARS,
     mammalian Pol II cassette, EBV episomal, SV40 replication, Sleeping Beauty,
     PiggyBac (mammalian), retroviral (MMLV/MSCV), Drosophila Ac5/P-element, etc.
     ANY such signal overrides bacterial backbone — bacterial selection markers
     (AmpR, KanR) reflect propagation, not host-of-use.

  2. Folder prior (Module_Library only): the file's parent directory under
     Module_Library_gb maps to a host category. Used only if rule 1 gave no
     strong signal; down-weighted so rule 1 always wins.

  3. Feature-level fallbacks: organism-hint features (plant promoters like
     CaMV 35S, insect-specific promoters like polh/p10, yeast selection like
     URA3/LEU2/HIS3/TRP1).

  4. If only bacterial replication + bacterial selection + bacterial MCS →
     bacterial.

  5. Otherwise → unknown.
"""
from __future__ import annotations
import re
from typing import Any

# Non-bacterial module signals — presence forces that host.
HOST_MODULE_SIGNALS = {
    # mammalian
    "lentiviral_payload": "mammalian",
    "lentiviral_cis_element": "mammalian",
    "mammalian_pol2_expression_cassette": "mammalian",
    "mammalian_lentiviral_expression_cassette": "mammalian",
    "mammalian_selection_cassette": "mammalian",
    "mammalian_replication": "mammalian",
    "aav_payload": "mammalian",
    "ebv_episomal_module": "mammalian",
    "sv40_replication_module": "mammalian",
    "floxed_cassette": "mammalian",
    "floxed_region": "mammalian",
    "lsl_cassette": "mammalian",
    "frt_flanked_cassette": "mammalian",
    "tet_regulator_cassette": "mammalian",
    "aid_degron_system": "mammalian",
    "fkbp_frb_dimerization": "mammalian",
    "integrase_landing_pad": "mammalian",

    # plant
    "tdna_module": "plant",

    # yeast
    "yeast_replication": "yeast",

    # bacterial / backbone-only (do NOT force bacterial — they're backbone)
    "bac_f_replicon": "bacterial_weak",  # weak hint
    "bacterial_replication": "bacterial_weak",
    "phage_replication": "bacterial_weak",

    # recombination cassettes are host-neutral
    "gateway_entry_cassette": None,
    "gateway_dest_cassette": None,
    "gateway_recombination": None,
    "lac_alpha_blue_white_module": "bacterial_weak",
    "ivt_cloning_cassette": None,
}

# Feature/payload-name keyword signals → host
# (keys are lowercased substrings; check feature name / payload_id / sseqid)
HOST_NAME_KEYWORDS = {
    "mammalian": {
        "lentivir", "lentiv", "lti ", "hiv-1", "hiv1", "hiv 1",
        "rre ", "rre-", "wpre", "cppt", "packaging signal",
        "5' ltr", "3' ltr", "sin ltr", "delta u3", "δu3", "deltau3",
        "5'ltr", "3'ltr",
        "mmlv", "mscv", "msv",
        "aav ", "aav-", "aavs", " itr", "5' itr", "3' itr",
        "ebna1", "ebna-1", "orip", "mini-orip",
        "sv40 poly", "bgh poly", "hgh poly",
        "cmv promoter", "cmv ", "cmv enhancer", "ef-1", "ef1a", "ef1α",
        "cag promoter", "cagg", "htlv", "sffv", "pgk promoter",
        "rsv ltr", "rous sarcoma",
        "cre recombinase", "cre/loxp", "loxp",
        "tet operator", "tet-on", "tet-off",
        "puror", "puromycin", "neor", "g418", "hygror", "hygromycin",
        "blasticidin", "zeocin", " bsd ",
        "piggybac", "piggy-bac", "pb itr",
        "sleeping beauty", "sb itr",
        "mclover", "mcherry", "egfp", "meyfp", "mneongreen",
        "dcas9", "cas9", "cas12", "cas13", "cas14", "abe", "cbe",
        "sgrna scaffold", "tracrrna", "pegrna",
    },
    "insect": {
        "baculovir", "baculo", "bacmid", "polh promoter", "polyhedrin",
        "p10 promoter", "op-ie", "ac5 promoter", "ac5.1",
        "ie-1 promoter", "copia ", "metallothionein promoter",
        "pfastbac", "pfast-bac", "bacrescue", "bac-to-bac",
        "tn7l", "tn7r",
    },
    "plant": {
        "t-dna", " tdna", "t-border", "right border", "left border",
        " lb ", " rb ", "rb repeat", "lb repeat",
        "pcambia", "cambia", "agrobacterium",
        "35s promoter", "camv 35s", "camv35s", "camv ",
        "nos promoter", "nos terminator", "nopaline synthase",
        "ocs promoter", "octopine synthase",
        "ubiquitin promoter zm", "zm ubi", "zmubi",
    },
    "yeast": {
        "ura3", "leu2", "his3", "trp1", "lys2", "ade2",
        "cen/ars", "cen6", "ars1", "ars4", "ars209",
        "2 micron", "2μ", "2u ori", "2-micron",
        "adh1 promoter", "gal1 promoter", "gal1-10",
        "gpd promoter", "tef promoter", "cyc1 terminator",
        "prs3", "prs4", "prs41", "prs42", "pyes",
    },
    "bacterial_weak": {
        "t7 promoter", "t3 promoter", "sp6 promoter",
        "lac promoter", "ptac", "trc promoter", "lacuv5",
        "lacz", "araC", "arabad", "arabinose",
    },
}

# Folder → host prior (Module_Library_gb subject directories)
FOLDER_PRIOR = {
    "Basic Cloning Vectors":           "bacterial",
    "CRISPR Plasmids":                 "mammalian",     # mostly mammalian
    "Fluorescent Protein Genes & Plasmids": "mammalian",
    "Gateway Cloning Vectors":         "bacterial",     # donor/dest are bacterial-propagated, target host varies
    "I.M.A.G.E. Consortium Plasmids":  "mammalian",     # cDNA clones, mostly mammalian
    "Insect Cell Vectors":             "insect",
    "Luciferase Vectors":              "mammalian",
    "Lucigen Vectors":                 "bacterial",
    "Mammalian Expression Vectors":    "mammalian",
}


def _tokens_contain_module(tokens: list[str], module_types: set[str]) -> bool:
    for t in tokens:
        if t.startswith("<MOD_OPEN:"):
            mt = t[len("<MOD_OPEN:"):].rstrip(">")
            if mt in module_types:
                return True
    return False


def _token_names(tokens: list[str]) -> list[str]:
    """Extract payload names from feature tokens."""
    names = []
    for t in tokens:
        if not (t.startswith("<") and ":" in t):
            continue
        # skip meta tokens
        prefix = t.split(":", 1)[0].lstrip("<")
        if prefix in ("BOS", "EOS", "TOPOLOGY", "LEN_BIN", "HOST", "SOURCE",
                      "ROTATION_IDX", "MOD_OPEN", "INT", "CLN", "VB_CAS",
                      "VB_FAMILY", "VB_SYS"):
            continue
        # Get the payload part
        parts = t.rstrip(">").split(":", 1)
        if len(parts) == 2:
            names.append(parts[1].lower())
    return names


def infer_host_priority(tokens: list[str],
                        folder_hint: str | None = None,
                        source: str = "") -> str:
    """
    Rule 1: Module-type signals override everything.
    Rule 2: Feature-name keywords (non-bacterial wins over bacterial_weak).
    Rule 3: Folder prior (Module_Library only; only when no non-bacterial signal).
    Rule 4: Bacterial if multiple bacterial_weak signals and no non-bacterial.
    Rule 5: Unknown.
    """
    # --- Rule 1: module signals ---
    mod_host_votes: dict[str, int] = {}
    for t in tokens:
        if not t.startswith("<MOD_OPEN:"):
            continue
        mt = t[len("<MOD_OPEN:"):].rstrip(">")
        host = HOST_MODULE_SIGNALS.get(mt)
        if host and host != "bacterial_weak":
            mod_host_votes[host] = mod_host_votes.get(host, 0) + 1

    # Strong non-bacterial module hit: prefer the most-frequent non-bacterial host
    if mod_host_votes:
        return max(mod_host_votes.items(), key=lambda kv: kv[1])[0]

    # --- Rule 2: feature-name keywords ---
    names = _token_names(tokens)
    joined = " ".join(names)
    name_votes: dict[str, int] = {}
    for host, kws in HOST_NAME_KEYWORDS.items():
        for kw in kws:
            if kw in joined:
                name_votes[host] = name_votes.get(host, 0) + 1

    # Non-bacterial keyword dominance
    non_bact = {h: c for h, c in name_votes.items()
                if h not in ("bacterial_weak",)}
    if non_bact:
        # Require ≥ 2 matches OR folder confirms to prevent a single stray word
        top = max(non_bact.items(), key=lambda kv: kv[1])
        if top[1] >= 2 or (folder_hint and FOLDER_PRIOR.get(folder_hint) == top[0]):
            return top[0]

    # --- Rule 3: folder prior ---
    if folder_hint:
        host = FOLDER_PRIOR.get(folder_hint)
        if host:
            return host

    # --- Rule 4: bacterial-weak signals ---
    bact_weak_module = _tokens_contain_module(tokens, {"bac_f_replicon",
                                                        "bacterial_replication",
                                                        "phage_replication",
                                                        "lac_alpha_blue_white_module"})
    bact_kws = name_votes.get("bacterial_weak", 0)
    if bact_weak_module or bact_kws >= 1:
        return "bacterial"

    return "unknown"
