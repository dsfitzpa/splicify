"""
Generate a plain-English summary of a plasmid's functional interactions.

Operates on the module-scoped interaction graph emitted by
`interaction_builder.build_interactions()`. Each interaction carries a
`rule_id` and `source_module`, letting us narrate relationships in the
user's mental model — cassette-by-cassette — rather than as a flat feature
list.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

SBO_LABELS = {
    "http://identifiers.org/biomodels.sbo/SBO:0000589": "genetic production",
    "http://identifiers.org/biomodels.sbo/SBO:0000183": "transcription",
    "http://identifiers.org/biomodels.sbo/SBO:0000184": "translation",
    "http://identifiers.org/biomodels.sbo/SBO:0000170": "stimulation",
    "http://identifiers.org/biomodels.sbo/SBO:0000169": "inhibition",
    "http://identifiers.org/biomodels.sbo/SBO:0000178": "cleavage",
    "http://identifiers.org/biomodels.sbo/SBO:0000177": "non-covalent binding",
    "http://identifiers.org/biomodels.sbo/SBO:0000168": "control",
    "http://identifiers.org/biomodels.sbo/SBO:0000182": "recombination",
    "http://identifiers.org/biomodels.sbo/SBO:0000179": "degradation",
}


def _participants_by_role(ix: Dict[str, Any]) -> Dict[str, List[Dict[str, Any]]]:
    out: Dict[str, List[Dict[str, Any]]] = {}
    for p in ix.get("participants", []) or []:
        role = (p.get("role") or "").lower()
        out.setdefault(role, []).append(p)
    return out


def _first(parts: List[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    return parts[0] if parts else None


def _span(p: Dict[str, Any]) -> str:
    s, e = p.get("start"), p.get("end")
    if s is None or e is None:
        return ""
    return f" ({s}–{e})"


# --------------------------------------------------------------------------- #
# Per-rule description templates
# --------------------------------------------------------------------------- #


def _describe_one(ix: Dict[str, Any]) -> str:
    roles = _participants_by_role(ix)
    rule = ix.get("rule_id", "")
    source = ix.get("source_module", "")

    stim = _first(roles.get("stimulator", []))
    tpl = _first(roles.get("template", []))
    mod = _first(roles.get("modifier", []))
    inh = _first(roles.get("inhibitor", []))
    react_list = roles.get("reactant", [])
    prod_list = roles.get("product", [])

    # Pol II cassette: promoter → CDS → polyA  (legacy rule, kept for compat)
    if rule == "INT-POL2-CAS-01" and stim and tpl:
        text = (f"**{stim['name']}**{_span(stim)} drives transcription of "
                f"**{tpl['name']}**{_span(tpl)}")
        if mod:
            text += f", terminated by **{mod['name']}**{_span(mod)}"
        if source == "mammalian_lentiviral_expression_cassette":
            text += " (lentiviral expression cassette)"
        return text + "."

    # Pol II upstream regulatory → CDS (peer interaction replacing INT-POL2-CAS-01)
    if rule == "INT-POL2-UR-01" and stim and tpl:
        return (f"**{stim['name']}**{_span(stim)} (upstream regulatory) drives "
                f"transcription of **{tpl['name']}**{_span(tpl)}.")

    # Pol II downstream regulatory modifies the CDS transcript
    if rule == "INT-POL2-DR-01" and tpl:
        mod_p = mod or (ix.get("participants", []) or [None])[0]
        if mod_p:
            return (f"**{mod_p['name']}**{_span(mod_p)} (downstream regulatory) "
                    f"terminates transcription and stabilizes the "
                    f"**{tpl['name']}**{_span(tpl)} mRNA.")

    # Lentiviral peer interactions: upstream drives payload, downstream modifies
    if rule == "INT-LENTI-UR-01" and stim and tpl:
        return (f"**{stim['name']}**{_span(stim)} drives transcription across "
                f"the lentiviral payload **{tpl['name']}**{_span(tpl)} "
                f"(integration-competent LTR-to-LTR region).")
    if rule == "INT-LENTI-DR-01" and tpl:
        mod_p = mod or (ix.get("participants", []) or [None])[0]
        if mod_p:
            return (f"**{mod_p['name']}**{_span(mod_p)} post-transcriptionally "
                    f"regulates the lentiviral payload "
                    f"**{tpl['name']}**{_span(tpl)} "
                    f"(WPRE stabilization and polyA termination).")

    # LAC-BW-01: lac promoter → lacZα
    if rule == "INT-LACBW-01" and stim and tpl:
        return (f"**{stim['name']}**{_span(stim)} transcribes the α-fragment "
                f"**{tpl['name']}**{_span(tpl)} that supports blue/white "
                f"α-complementation.")

    # LAC-BW-02: operator inhibits promoter
    if rule == "INT-LACBW-02" and inh and mod:
        return (f"**{inh['name']}**{_span(inh)} is bound by LacI to repress "
                f"**{mod['name']}**{_span(mod)}; IPTG relieves repression.")

    # LAC-BW-03: LacI binds operator
    if rule == "INT-LACBW-03" and stim and tpl:
        return (f"**{stim['name']}**{_span(stim)} encodes LacI, which binds "
                f"**{tpl['name']}**{_span(tpl)} to maintain the lac repressor loop.")

    # LAC-BW-04: MCS disrupts lacZα
    if rule == "INT-LACBW-04" and inh and mod:
        return (f"Inserting into **{inh['name']}**{_span(inh)} disrupts "
                f"**{mod['name']}**{_span(mod)} → white colony on X-gal/IPTG "
                f"(blue/white screening).")

    # Recombination — loxP / FRT / att / integrase
    if rule and rule.startswith("INT-REC-"):
        cargo = mod
        enzyme = None
        for p in ix.get("participants", []):
            if p.get("external") and p.get("role") == "stimulator":
                enzyme = p
                break
        flanks = react_list
        if cargo and enzyme and len(flanks) >= 2:
            return (f"**{cargo['name']}**{_span(cargo)} is flanked by "
                    f"{flanks[0]['name']} and {flanks[1]['name']} sites; "
                    f"{enzyme['name']} catalyzes excision, inversion, or integration.")

    # CDS P2A / T2A ribosomal skip — intrinsic cleavage biology
    if rule == "INT-CDS-2A-01" and stim and len(prod_list) >= 2:
        a, b = prod_list[0], prod_list[1]
        return (f"The **{stim['name']}** ribosomal-skip peptide splits the ORF "
                f"co-translationally, releasing **{a['name']}**{_span(a)} and "
                f"**{b['name']}**{_span(b)} as two independent polypeptides.")

    # Cassette-coupled genetic_production of one of two 2A products
    if rule == "INT-CDS-2A-02" and stim and tpl:
        text = (f"**{stim['name']}**{_span(stim)} drives expression of "
                f"**{tpl['name']}**{_span(tpl)} — one of two independent "
                f"proteins produced from a single 2A-linked ORF; both products "
                f"share this promoter and its polyA signal")
        if mod:
            text += f" (**{mod['name']}**{_span(mod)})"
        return text + "."

    # Insulator block
    if rule == "INT-INS-01" and mod:
        return (f"**{mod['name']}**{_span(mod)} brackets an insulator-bounded "
                f"expression region, suppressing enhancer reach across the boundary.")


    # Bacterial selection — promoter drives selection CDS (no wrapper module)
    if rule == "INT-BSEL-01" and stim and tpl:
        return (f"**{stim['name']}**{_span(stim)} drives expression of "
                f"**{tpl['name']}**{_span(tpl)} (bacterial selection).")

    # Bacterial selection — terminator regulates selection CDS (peer interaction)
    if rule == "INT-BSEL-TERM-01" and tpl:
        mod_p = mod or (ix.get("participants", []) or [None])[0]
        if mod_p:
            return (f"**{mod_p['name']}**{_span(mod_p)} terminates transcription "
                    f"of **{tpl['name']}**{_span(tpl)} (bacterial selection).")
    # Fallback
    label = SBO_LABELS.get(ix.get("sbo_term", ""), ix.get("interaction_type", "interaction"))
    names = ", ".join(p.get("name", "?") for p in ix.get("participants", []))
    return f"{label.capitalize()} involving {names}."


# --------------------------------------------------------------------------- #
# Top-line summary
# --------------------------------------------------------------------------- #


def _summary_paragraph(interactions: List[Dict[str, Any]]) -> str:
    if not interactions:
        return ("No within-module functional interactions were inferred — the plasmid's "
                "detected modules don't carry enough submodule information yet to map "
                "explicit relationships.")

    # Count *distinct* modules per source type. Two interactions come from
    # the same physical module when they share source_module AND any
    # participant span overlaps — exact match OR range containment in
    # either direction (so an interaction whose protein span (5110,9217)
    # still counts as the same module as one with the enclosing cds_module
    # span (5110,9949)).
    from collections import defaultdict

    def _overlap_or_contain(a, b):
        """True if any (start,end) in set a overlaps any (start,end) in b."""
        for sa, ea in a:
            if sa is None or ea is None:
                continue
            for sb, eb in b:
                if sb is None or eb is None:
                    continue
                if sa == sb and ea == eb:
                    return True
                if sa >= sb and ea <= eb:   # a contained in b
                    return True
                if sb >= sa and eb <= ea:   # b contained in a
                    return True
                if sa < eb and sb < ea:     # any partial intersection
                    return True
        return False

    groups: Dict[str, List[set]] = defaultdict(list)  # src -> list of span-sets
    for ix in interactions:
        src = ix.get("source_module", "other")
        spans = {
            (p.get("start"), p.get("end"))
            for p in ix.get("participants", [])
            if p.get("start") is not None
        }
        matched = False
        for g in groups[src]:
            if _overlap_or_contain(g, spans):
                g.update(spans)
                matched = True
                break
        if not matched:
            groups[src].append(spans)
    by_source: Dict[str, int] = {k: len(v) for k, v in groups.items()}

    # Human-readable module-type labels
    pretty = {
        "mammalian_pol2_expression_cassette": "mammalian Pol II expression cassette",
        "mammalian_lentiviral_expression_cassette": "lentiviral expression cassette",
        "lac_alpha_blue_white_module": "lacZα blue/white screening module",
        "floxed_region": "floxed (loxP-flanked) cassette",
        "lsl_cassette": "LSL (loxP-STOP-loxP) cassette",
        "frt_flanked_cassette": "FRT-flanked cassette",
        "gateway_entry_cassette": "Gateway entry cassette",
        "gateway_dest_cassette": "Gateway destination cassette",
        "gateway_recombination": "Gateway recombination site",
        "integrase_landing_pad": "phage integrase landing pad",
        "insulated_expression_block": "insulator-bounded expression block",
        "cds_orf": "coding-ORF internal relationship",
    }

    # Roll up duplicates into counts, preserving declaration order.
    ordered: List[str] = []
    seen: Dict[str, int] = {}
    for src in by_source:
        label = pretty.get(src, src.replace("_", " "))
        seen[label] = seen.get(label, 0) + by_source[src]
        if label not in ordered:
            ordered.append(label)

    parts = [
        f"{seen[l]} {l}{'s' if seen[l] != 1 and not l.endswith('s') else ''}"
        for l in ordered
    ]
    joined = ", ".join(parts)
    return (f"This plasmid carries {joined}. "
            f"The interactions below describe how the parts inside each module "
            f"cooperate to carry out the module's function.")


# --------------------------------------------------------------------------- #
# Public entry point
# --------------------------------------------------------------------------- #


def describe_interactions(
    interactions: List[Dict[str, Any]],
    *,
    plasmid_name: Optional[str] = None,
) -> Dict[str, Any]:
    """Build a structured description for the UI.

    Returns dict with keys: summary, bullets, counts, markdown, by_module.
    """
    summary = _summary_paragraph(interactions)

    # Group bullets by source_module so the UI can render a nested, cassette-
    # level narrative if it wants to.
    by_module: Dict[str, List[str]] = {}
    for ix in interactions:
        by_module.setdefault(ix.get("source_module", "other"), []).append(
            f"- {_describe_one(ix)}"
        )

    # Flat bullets list (kept for back-compat with the existing sidebar)
    bullets: List[str] = []
    for _mod, group in by_module.items():
        bullets.extend(group)

    counts: Dict[str, int] = {}
    for ix in interactions:
        key = SBO_LABELS.get(ix.get("sbo_term", ""), ix.get("interaction_type", "other"))
        counts[key] = counts.get(key, 0) + 1

    # Build a cassette-section markdown
    md_parts = []
    if plasmid_name:
        md_parts.append(f"## {plasmid_name}")
    md_parts.append(summary)
    for src, group in by_module.items():
        label = src.replace("_", " ")
        md_parts.append(f"\n**{label}**\n" + "\n".join(group))
    md = "\n\n".join(md_parts)

    return {
        "summary": summary,
        "bullets": bullets,
        "by_module": by_module,
        "counts": counts,
        "markdown": md,
    }


__all__ = ["describe_interactions"]
