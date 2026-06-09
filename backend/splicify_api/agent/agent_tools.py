"""
Tool implementations for the AIPlasmidDesign agent (v2).

Adds: score_sanger_primer, analyze_design_intent, verify_assembly,
compare_to_choice. Existing tools unchanged in behaviour but
simulate_assembly now registers the assembled product as a new
attachment so downstream tools can reference it.
"""
from __future__ import annotations

import hashlib
import logging
import re
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

logger = logging.getLogger("agent.tools")


# --- attachment registry -------------------------------------------------

@dataclass
class Attachment:
    attachment_id: str
    name: str
    sequence: str
    circular: bool = True
    role: str = "inventory"  # "target" | "inventory" | "product"


@dataclass
class MCQChoice:
    letter: str
    text: str  # may be DNA (long) or arbitrary prose


@dataclass
class AttachmentRegistry:
    items: Dict[str, Attachment] = field(default_factory=dict)
    choices: Dict[str, MCQChoice] = field(default_factory=dict)
    _next_product_id: int = 1

    def add(self, att: Attachment) -> None:
        self.items[att.attachment_id] = att

    def get(self, attachment_id: str) -> Optional[Attachment]:
        return self.items.get(attachment_id)

    def register_product(self, name: str, sequence: str, circular: bool = True) -> str:
        aid = f"att_product_{self._next_product_id}"
        self._next_product_id += 1
        self.add(Attachment(attachment_id=aid, name=name, sequence=sequence,
                            circular=circular, role="product"))
        return aid

    def add_choice(self, letter: str, text: str) -> None:
        self.choices[letter.upper()] = MCQChoice(letter=letter.upper(), text=text or "")

    def public_summary(self) -> List[Dict[str, Any]]:
        return [
            {"attachment_id": a.attachment_id, "name": a.name,
             "length_bp": len(a.sequence),
             "topology": "circular" if a.circular else "linear",
             "role": a.role}
            for a in self.items.values()
        ]


# --- helpers --------------------------------------------------------------

def _clean_seq(s: str) -> str:
    return re.sub(r"[^ACGTNacgtn]", "", s).upper()


def extract_seq_from_genbank(gb_text: str) -> str:
    if "ORIGIN" not in gb_text:
        return _clean_seq(gb_text)
    after = gb_text.split("ORIGIN", 1)[1].split("//", 1)[0]
    return _clean_seq(after)


def extract_name_from_genbank(gb_text: str, fallback: str) -> str:
    m = re.search(r"^LOCUS\s+(\S+)", gb_text, re.MULTILINE)
    return m.group(1) if m else fallback


def _revcomp(s: str) -> str:
    comp = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    return "".join(comp.get(b, "N") for b in reversed(s))


def _trim_annotation(ann: Dict[str, Any], max_features: int = 60) -> Dict[str, Any]:
    out: Dict[str, Any] = {"length_bp": ann.get("sequence_length") or ann.get("length")}
    feats = ann.get("annotations") or ann.get("features") or []
    out["features"] = [
        {"name": f.get("Feature") or f.get("name"),
         "type": f.get("Type") or f.get("type"),
         "start": f.get("qstart") if f.get("qstart") is not None else f.get("start"),
         "end": f.get("qend") if f.get("qend") is not None else f.get("end"),
         "strand": f.get("sframe") or f.get("strand")}
        for f in feats[:max_features]
    ]
    if "modules" in ann:
        out["modules"] = [
            {"name": m.get("name"), "type": m.get("type"),
             "start": m.get("start"), "end": m.get("end")}
            for m in (ann.get("modules") or [])[:30]
        ]
    if "cloning_features" in ann:
        # cloning_features evolved from list[dict] -> dict with a "features"
        # list. Tolerate both shapes so genomic-kind annotations (which the
        # v2 router runs on linear gene records) don't trip a dict[slice]
        # KeyError on the slicing below.
        cf_value = ann.get("cloning_features") or []
        if isinstance(cf_value, dict):
            cf_value = cf_value.get("features", []) or []
        cf_clean = []
        for cf in cf_value[:60]:
            if isinstance(cf, dict):
                cf_clean.append({k: v for k, v in cf.items()
                                 if k not in ("sequence", "context_sequence", "raw")})
        out["cloning_features"] = cf_clean
    return out


def _strip_sequences(obj: Any) -> Any:
    SEQ_KEYS = {"sequence", "seq", "genbank", "gb", "raw_sequence",
                "fragment_sequence", "product_sequence", "context_sequence"}
    if isinstance(obj, dict):
        return {k: ("[redacted]" if k in SEQ_KEYS else _strip_sequences(v))
                for k, v in obj.items()}
    if isinstance(obj, list):
        return [_strip_sequences(x) for x in obj]
    if isinstance(obj, str):
        if len(obj) >= 60 and re.fullmatch(r"[ACGTNacgtn\s]+", obj or ""):
            return f"[redacted DNA, {len(obj.strip())} chars]"
    return obj


def _gb_for_attachment(att: Attachment) -> str:
    seq = att.sequence.lower()
    topo = "circular" if att.circular else "linear"
    lines = [f"LOCUS       {att.name[:16]:<16} {len(seq):>8} bp    DNA     {topo:<8} SYN 01-JAN-2026",
             f"DEFINITION  {att.name}", "FEATURES             Location/Qualifiers",
             f"     source          1..{len(seq)}",
             '                     /organism="synthetic DNA construct"',
             '                     /mol_type="other DNA"', "ORIGIN"]
    for start in range(0, len(seq), 60):
        chunk = seq[start:start + 60]
        groups = [chunk[i:i + 10] for i in range(0, len(chunk), 10)]
        lines.append(f"{start+1:>9} " + " ".join(groups))
    lines.append("//")
    return "\n".join(lines) + "\n"


def _resolve_feature_position(annotation: Dict[str, Any], feature_name: str) -> Optional[int]:
    """Locate the midpoint of a named feature in an annotation result."""
    feats = annotation.get("annotations") or annotation.get("features") or []
    needle = (feature_name or "").lower()
    if not needle:
        return None
    best = None
    for f in feats:
        n = (f.get("Feature") or f.get("name") or "").lower()
        if needle in n or n in needle:
            s = f.get("qstart") if f.get("qstart") is not None else f.get("start")
            e = f.get("qend") if f.get("qend") is not None else f.get("end")
            if s is not None and e is not None:
                best = (int(s) + int(e)) // 2
                break
    if best is None and (annotation.get("modules") or []):
        for m in annotation["modules"]:
            n = (m.get("name") or "").lower()
            if needle in n or n in needle:
                if m.get("start") is not None and m.get("end") is not None:
                    best = (int(m["start"]) + int(m["end"])) // 2
                    break
    return best


# --- tools ----------------------------------------------------------------

async def tool_annotate_attachment(args: Dict[str, Any], reg: AttachmentRegistry) -> Dict[str, Any]:
    from ..annotation_cache import annotate_cached
    aid = args.get("attachment_id")
    att = reg.get(aid) if aid else None
    if not att:
        return {"error": f"unknown attachment_id: {aid!r}. "
                         f"Available: {[a.attachment_id for a in reg.items.values()]}"}
    try:
        ann = await annotate_cached(att.sequence, circular=att.circular, depth="full")
    except Exception as e:
        logger.exception("annotate_cached failed")
        return {"error": f"annotation failed: {type(e).__name__}: {e}"}
    out = _trim_annotation(ann)
    out["attachment_id"] = aid
    out["name"] = att.name
    return out


async def _annotate_full(att: Attachment) -> Optional[Dict[str, Any]]:
    from ..annotation_cache import annotate_cached
    try:
        return await annotate_cached(att.sequence, circular=att.circular, depth="full")
    except Exception as e:
        logger.exception("annotate failed for %s: %s", att.attachment_id, e)
        return None


async def tool_simulate_assembly(args: Dict[str, Any], reg: AttachmentRegistry) -> Dict[str, Any]:
    from ..chat import chat as chat_endpoint
    from fastapi import UploadFile
    from io import BytesIO
    from starlette.datastructures import Headers

    instruction = args.get("instruction", "")
    target_id = args.get("target_attachment_id")
    inv_ids: List[str] = args.get("inventory_attachment_ids") or []

    def _upload(name: str, gb_text: str) -> UploadFile:
        return UploadFile(filename=name, file=BytesIO(gb_text.encode()),
                          headers=Headers({"content-type": "chemical/seq-na-genbank"}))

    target_att = reg.get(target_id) if target_id else None
    file = _upload(f"{target_att.name}.gb", _gb_for_attachment(target_att)) if target_att else None

    inv_atts = [a for a in (reg.get(i) for i in inv_ids) if a is not None]
    inv_files = [_upload(f"{a.name}.gb", _gb_for_attachment(a)) for a in inv_atts] or None

    try:
        result = await chat_endpoint(
            message=instruction, session_id="agent",
            include_ai_explanation="true", describe_plasmid_intent="false",
            file=file, inventory_files=inv_files,
        )
    except Exception as e:
        logger.exception("simulate_assembly failed")
        return {"error": f"chat dispatch failed: {type(e).__name__}: {e}"}

    if hasattr(result, "body"):
        import json as _json
        result = _json.loads(result.body)

    # Find an emitted GenBank file in the response and register it as a product.
    # output_builders.py emits files as {"fileName": "..._assembled.gb", "dataBase64": "..."}.
    # Older code paths may use lowercase "filename" / "content"; we accept both.
    product_aid = None
    files_raw = (result or {}).get("files")
    files = files_raw if isinstance(files_raw, list) else []
    import base64 as _b64

    def _extract_gb_text(f: dict) -> str:
        for b64_key in ("dataBase64", "data_base64", "content_b64", "b64"):
            v = f.get(b64_key)
            if isinstance(v, str) and v:
                try:
                    return _b64.b64decode(v).decode("utf-8", errors="ignore")
                except Exception:
                    continue
        for txt_key in ("content", "data", "text"):
            v = f.get(txt_key)
            if isinstance(v, str) and v:
                return v
        return ""

    files_emitted_names: list = []
    for f in files:
        if not isinstance(f, dict):
            continue
        fname = f.get("fileName") or f.get("filename") or f.get("name") or ""
        files_emitted_names.append(fname)
        # Prioritise actual assembled-product files.
        if fname and ("assembled" not in fname.lower()) and ("product" not in fname.lower()):
            continue
        gb_text = _extract_gb_text(f)
        if gb_text and "ORIGIN" in gb_text:
            seq = extract_seq_from_genbank(gb_text)
            if seq and len(seq) >= 100:
                product_aid = reg.register_product(
                    name=(fname or "assembly_product").replace(".gb", ""),
                    sequence=seq, circular=True,
                )
                break
    # Fallback: if no assembled.gb but some other GenBank-shaped file exists, take the longest one.
    if product_aid is None:
        best_seq = ""
        best_name = ""
        for f in files:
            if not isinstance(f, dict):
                continue
            gb_text = _extract_gb_text(f)
            if gb_text and "ORIGIN" in gb_text:
                seq = extract_seq_from_genbank(gb_text)
                if len(seq) > len(best_seq):
                    best_seq = seq
                    best_name = f.get("fileName") or f.get("filename") or "product"
        if best_seq and len(best_seq) >= 100:
            product_aid = reg.register_product(
                name=best_name.replace(".gb", ""),
                sequence=best_seq, circular=True,
            )

    return _strip_sequences({
        "reply":                  (result or {}).get("reply", ""),
        "intent":                 (result or {}).get("intent"),
        "product_attachment_id":  product_aid,
        "product_registered":     bool(product_aid),
        "files_emitted":          files_emitted_names,
        "predesign_verdict":      ((result or {}).get("predesign_context") or {}).get("predesign_evaluation"),
    })


async def tool_digest_plasmid(args: Dict[str, Any], reg: AttachmentRegistry) -> Dict[str, Any]:
    aid = args.get("attachment_id")
    att = reg.get(aid) if aid else None
    if not att:
        return {"error": f"unknown attachment_id: {aid!r}"}
    enzymes: List[str] = args.get("enzymes") or []
    if not enzymes:
        return {"error": "at least one enzyme required"}
    try:
        from Bio.Restriction import RestrictionBatch
        from Bio.Seq import Seq
    except ImportError:
        return {"error": "BioPython not installed"}
    try:
        batch = RestrictionBatch(enzymes)
    except Exception as e:
        return {"error": f"unknown enzyme: {e}"}

    s = Seq(att.sequence)
    cut_map = batch.search(s, linear=not att.circular)
    by_enzyme = {str(e): sorted(p) for e, p in cut_map.items()}
    cuts = sorted({p for ps in cut_map.values() for p in ps})

    if att.circular and cuts:
        rotated = cuts + [cuts[0] + len(att.sequence)]
        fragments = [rotated[i+1] - rotated[i] for i in range(len(cuts))]
    elif cuts:
        fragments = ([cuts[0]] +
                     [cuts[i+1] - cuts[i] for i in range(len(cuts) - 1)] +
                     [len(att.sequence) - cuts[-1]])
    else:
        fragments = [len(att.sequence)]
    fragments.sort(reverse=True)

    return {
        "attachment_id": aid, "name": att.name, "length_bp": len(att.sequence),
        "cut_positions_by_enzyme": by_enzyme, "n_cuts_total": len(cuts),
        "fragment_lengths_bp": fragments,
    }


async def tool_find_primer_binding(args: Dict[str, Any], reg: AttachmentRegistry) -> Dict[str, Any]:
    aid = args.get("template_attachment_id")
    att = reg.get(aid) if aid else None
    if not att:
        return {"error": f"unknown template_attachment_id: {aid!r}"}
    primers: List[str] = args.get("primers") or []
    if not primers:
        return {"error": "primers list is empty"}

    template = att.sequence
    haystack = template + (template if att.circular else "")
    hits = []
    for p in primers:
        p_clean = _clean_seq(p)
        if not p_clean:
            continue
        anneal = p_clean[-20:] if len(p_clean) >= 20 else p_clean
        fwd = haystack.find(anneal) % len(template) if anneal in haystack else -1
        rc = _revcomp(anneal)
        rev = haystack.find(rc) % len(template) if rc in haystack else -1
        hits.append({
            "primer_length": len(p_clean), "primer_3prime_anneal": anneal,
            "fwd_pos": fwd if fwd >= 0 else None,
            "rev_pos": rev if rev >= 0 else None,
            "binds": (fwd >= 0) or (rev >= 0),
        })
    return {"attachment_id": aid, "name": att.name, "hits": hits}


async def tool_score_sanger_primer(args: Dict[str, Any], reg: AttachmentRegistry) -> Dict[str, Any]:
    aid = args.get("template_attachment_id")
    att = reg.get(aid) if aid else None
    if not att:
        return {"error": f"unknown template_attachment_id: {aid!r}"}
    primers: List[str] = args.get("primers") or []
    if not primers:
        return {"error": "primers list is empty"}

    target_pos = args.get("target_position")
    target_name = args.get("target_feature_name")

    if target_pos is None and target_name:
        ann = await _annotate_full(att)
        if ann is None:
            return {"error": "annotation failed; cannot resolve target_feature_name"}
        target_pos = _resolve_feature_position(ann, target_name)
        if target_pos is None:
            available = [(f.get("Feature") or f.get("name"))
                         for f in (ann.get("annotations") or ann.get("features") or [])][:30]
            return {"error": f"feature {target_name!r} not found",
                    "available_features": available}

    if target_pos is None:
        return {"error": "must provide either target_feature_name or target_position"}

    try:
        from ..sanger_scoring import score_sanger_primer as _scorer
    except Exception as e:
        return {"error": f"sanger_scoring not available: {type(e).__name__}: {e}"}

    # Pass each primer as both forward and reverse — caller doesn't always know.
    primer_specs = []
    for p in primers:
        p_clean = _clean_seq(p)
        primer_specs.append({"sequence": p_clean, "direction": "forward",
                             "name": f"{p_clean[:8]}...fwd"})
        primer_specs.append({"sequence": p_clean, "direction": "reverse",
                             "name": f"{p_clean[:8]}...rev"})
    try:
        scored = _scorer(att.sequence, primer_specs, int(target_pos))
    except Exception as e:
        logger.exception("score_sanger_primer raised")
        return {"error": f"scoring failed: {type(e).__name__}: {e}"}

    # Aggregate per-primer best (fwd vs rev) so the agent sees one row per primer.
    best_by_primer: Dict[str, Dict[str, Any]] = {}
    for spec, s in zip(primer_specs, scored):
        key = spec["sequence"]
        if key not in best_by_primer or s.get("overall_score", 0) > best_by_primer[key].get("overall_score", -1):
            best_by_primer[key] = {
                "primer":         key,
                "best_direction": spec["direction"],
                "overall_score":  s.get("overall_score"),
                "rating":         s.get("rating"),
                "warnings":       s.get("warnings"),
            }
    return {
        "attachment_id": aid, "target_position": target_pos,
        "scored": list(best_by_primer.values()),
    }


async def tool_analyze_design_intent(args: Dict[str, Any], reg: AttachmentRegistry) -> Dict[str, Any]:
    user_message = args.get("user_message") or ""
    if not user_message.strip():
        return {"error": "user_message is empty"}
    try:
        from ..intent import parse_intent
        from ..target_from_inventory_router import analyze_design_intent
    except Exception as e:
        return {"error": f"intent module not available: {type(e).__name__}: {e}"}
    try:
        intent_result = await parse_intent(
            user_message, has_target=False, has_inventory=False,
            seq_count=0, redacted_message=user_message,
        )
        completeness = intent_result.get("design_completeness")
        if completeness is None:
            completeness = analyze_design_intent(intent_result, user_message)
    except Exception as e:
        logger.exception("analyze_design_intent failed")
        return {"error": f"{type(e).__name__}: {e}"}

    return _strip_sequences({
        "intent":              intent_result.get("intent"),
        "kb_resolved_parts":   [
            i.get("feature_name") or i.get("name")
            for i in (intent_result.get("kb_resolved") or {}).get("identified") or []
        ],
        "design_completeness": completeness,
    })


async def tool_verify_assembly(args: Dict[str, Any], reg: AttachmentRegistry) -> Dict[str, Any]:
    aid = args.get("attachment_id")
    att = reg.get(aid) if aid else None
    if not att:
        return {"error": f"unknown attachment_id: {aid!r}"}
    ann = await _annotate_full(att)
    if ann is None:
        return {"error": "annotation failed"}
    try:
        from ..target_from_inventory_router import verify_target_design
    except Exception as e:
        return {"error": f"verify_target_design not available: {type(e).__name__}: {e}"}

    annotations = ann.get("annotations") or ann.get("features") or []
    modules = ann.get("modules") or []
    interactions = ann.get("interactions") or []
    try:
        verdict = verify_target_design(
            target_sequence=att.sequence,
            target_annotations=annotations,
            target_modules=modules,
            target_interactions=interactions,
            target_name=att.name,
        )
    except Exception as e:
        logger.exception("verify_target_design raised")
        return {"error": f"verify failed: {type(e).__name__}: {e}"}
    return _strip_sequences({"attachment_id": aid, "name": att.name, **verdict})


async def tool_compare_to_choice(args: Dict[str, Any], reg: AttachmentRegistry) -> Dict[str, Any]:
    aid = args.get("attachment_id")
    letter = (args.get("choice_letter") or "").upper()
    att = reg.get(aid) if aid else None
    if not att:
        return {"error": f"unknown attachment_id: {aid!r}"}
    if not letter or letter not in reg.choices:
        return {"error": f"unknown choice_letter: {letter!r}",
                "available_choices": sorted(reg.choices.keys())}
    choice = reg.choices[letter]
    choice_seq = _clean_seq(choice.text)
    if not choice_seq:
        return {"error": "this choice does not look like DNA",
                "choice_text_preview": (choice.text or "")[:120]}

    att_seq = att.sequence
    length_matches = len(att_seq) == len(choice_seq)

    def _hash(s: str) -> str:
        return hashlib.sha1(s.encode()).hexdigest()[:16]

    if not length_matches:
        return {"match": False, "length_matches": False,
                "attachment_length_bp": len(att_seq),
                "choice_length_bp": len(choice_seq)}

    # Try exact match in any rotation, both strands.
    n = len(att_seq)
    candidates = [att_seq, _revcomp(att_seq)]
    target_hash = _hash(choice_seq)
    for cand in candidates:
        doubled = cand + cand
        for off in range(n):
            window = doubled[off:off + n]
            if _hash(window) == target_hash:
                return {"match": True, "length_matches": True,
                        "best_rotation_offset": off,
                        "strand": "forward" if cand is att_seq else "reverse_complement"}

    # No exact match — give a quick mismatch estimate using the best rotation
    # of the forward strand (cheap bag-of-kmers comparison).
    def _kmer_set(s: str, k: int = 12) -> set:
        return {s[i:i+k] for i in range(0, len(s) - k + 1, max(1, len(s) // 200))}
    a_kmers = _kmer_set(att_seq)
    c_kmers = _kmer_set(choice_seq)
    overlap = len(a_kmers & c_kmers) / max(1, len(a_kmers | c_kmers))
    return {"match": False, "length_matches": True,
            "kmer_jaccard": round(overlap, 3),
            "interpretation": "high jaccard (>0.9) suggests near-match; <0.5 suggests different sequence"}


# --- deterministic in-silico Golden Gate ---------------------------------

_TYPE_IIS_RECOGNITION = {
    "Esp3I":  "CGTCTC",
    "BsmBI":  "CGTCTC",
    "BsaI":   "GGTCTC",
    "BsaI-HFv2": "GGTCTC",
    "BbsI":   "GAAGAC",
    "BbsI-HF": "GAAGAC",
    "SapI":   "GCTCTTC",
    "PaqCI":  "CACCTGC",
    "AarI":   "CACCTGC",
}

# Each of these enzymes cuts N1 nt away on the top strand, leaving a 4-nt
# 5' overhang (SapI is 1 + 3-nt overhang). For Esp3I/BsmBI/BsaI/BbsI:
# 5'...CGTCTC N1 NNNN...3'
# 3'...GCAGAG N5 NNNN...5'  → top cut after the spacer; 4-nt overhang.
_ENZYME_PROFILE = {
    "Esp3I":  {"site": "CGTCTC", "spacer": 1, "overhang": 4},
    "BsmBI":  {"site": "CGTCTC", "spacer": 1, "overhang": 4},
    "BsaI":   {"site": "GGTCTC", "spacer": 1, "overhang": 4},
    "BsaI-HFv2": {"site": "GGTCTC", "spacer": 1, "overhang": 4},
    "BbsI":   {"site": "GAAGAC", "spacer": 2, "overhang": 4},
    "BbsI-HF": {"site": "GAAGAC", "spacer": 2, "overhang": 4},
    "SapI":   {"site": "GCTCTTC", "spacer": 1, "overhang": 3},
}


def _find_overhang_fragment(seq: str, enzyme: str, circular: bool) -> Optional[Dict[str, Any]]:
    """For a Type-IIS digest of a single circular plasmid (which has TWO
    inward-facing recognition sites flanking the insert), compute the
    insert fragment + its left/right 4-nt overhangs.

    Returns {sequence, left_overhang, right_overhang, is_full_plasmid}.
    Returns None on ambiguity (>2 sites) or single-site (uncuttable).
    """
    profile = _ENZYME_PROFILE.get(enzyme)
    if profile is None:
        return None
    site = profile["site"]
    spacer = profile["spacer"]
    overhang_len = profile["overhang"]
    rc_site = _revcomp(site)

    seq = seq.upper()
    haystack = seq + seq if circular else seq

    # Find every recognition occurrence (forward and reverse).
    fwd = []
    i = 0
    while True:
        i = haystack.find(site, i)
        if i < 0:
            break
        # Cut position on the top strand: after the spacer
        cut = i + len(site) + spacer
        fwd.append({"strand": "+", "site_start": i, "cut": cut})
        i += 1
    rev = []
    i = 0
    while True:
        i = haystack.find(rc_site, i)
        if i < 0:
            break
        # - strand recognition (rc_site at top position i): the TOP-strand
        # cut releasing the 4-nt 5' overhang is at i - spacer - overhang_len.
        cut = i - spacer - overhang_len
        if cut >= 0:
            rev.append({"strand": "-", "site_start": i, "cut": cut})
        i += 1

    # Deduplicate by % seq length
    n = len(seq)

    def _norm(c: int) -> int:
        return c % n if circular else c

    sites = []
    for s in fwd + rev:
        c = _norm(s["cut"])
        sites.append({**s, "cut_norm": c})

    # We want a clean 2-cut release (one + and one -) in a circular plasmid.
    plus_cuts = sorted({s["cut_norm"] for s in sites if s["strand"] == "+"})
    minus_cuts = sorted({s["cut_norm"] for s in sites if s["strand"] == "-"})
    if not plus_cuts or not minus_cuts:
        return None

    # Take the inner-facing pair: + cut comes first walking forward,
    # - cut second; the segment between them on the + strand is the insert.
    # For a single insert per part, exactly one + + one - inside it.
    if len(plus_cuts) > 1 or len(minus_cuts) > 1:
        # Could be a multi-fragment input — for simplicity reject.
        return None

    p_cut = plus_cuts[0]
    m_cut = minus_cuts[0]

    # Build doubled sequence for clean rotation handling.
    doubled = seq + seq

    # Determine fragment span on + strand: cut sites point inward.
    # For + recognition at position i, the cut produces: ...spacer | overhang...
    # Insert end starts at p_cut (inclusive of overhang on the LEFT of insert).
    # For - recognition at j, cut on top strand is at m_cut, and the
    # overhang occupies [m_cut, m_cut + overhang_len) on the top strand,
    # belonging to the insert's RIGHT end.
    #
    # Walk from p_cut to m_cut + overhang_len modulo n.
    if m_cut >= p_cut:
        start = p_cut
        end = m_cut + overhang_len
        if end <= n:
            insert = seq[start:end]
        else:
            insert = doubled[start:end]
    else:
        # Wraps around origin
        start = p_cut
        end = n + m_cut + overhang_len
        insert = doubled[start:end]

    if len(insert) <= overhang_len * 2:
        return None
    left_oh = insert[:overhang_len]
    right_oh = insert[-overhang_len:]
    return {
        "sequence": insert,
        "left_overhang": left_oh,
        "right_overhang": right_oh,
        "length_bp": len(insert),
    }


async def tool_golden_gate_assemble(
    args: Dict[str, Any], reg: AttachmentRegistry,
) -> Dict[str, Any]:
    """Deterministic in-silico Golden Gate assembly: digest each input
    plasmid with the type-IIS enzyme, find compatible 4-nt overhangs,
    ligate in a unique cyclic order, register the product as a new
    attachment.
    """
    enzyme = (args.get("enzyme") or "").strip()
    if enzyme not in _ENZYME_PROFILE:
        return {"error": f"unknown / unsupported enzyme: {enzyme!r}",
                "supported": sorted(_ENZYME_PROFILE.keys())}

    aids: List[str] = args.get("attachment_ids") or []
    if len(aids) < 2:
        return {"error": "need at least 2 attachment_ids"}
    parts = []
    for aid in aids:
        att = reg.get(aid)
        if not att:
            return {"error": f"unknown attachment_id: {aid!r}"}
        frag = _find_overhang_fragment(att.sequence, enzyme, circular=att.circular)
        if frag is None:
            return {"error": f"could not find a unique 2-cut Type-IIS release in {att.name} ({aid}) with {enzyme}"}
        parts.append({"aid": aid, "name": att.name, **frag})

    # Build directed graph: for each part, its right_overhang must equal
    # the next part's left_overhang.
    by_left: Dict[str, List[int]] = {}
    for i, p in enumerate(parts):
        by_left.setdefault(p["left_overhang"], []).append(i)

    # Try to walk a unique cycle starting from each part.
    n = len(parts)
    for start in range(n):
        order = [start]
        used = {start}
        ok = True
        for _ in range(n - 1):
            cur = parts[order[-1]]
            nxts = [j for j in by_left.get(cur["right_overhang"], []) if j not in used]
            if len(nxts) != 1:
                ok = False
                break
            order.append(nxts[0])
            used.add(nxts[0])
        if not ok:
            continue
        # Final closure: last part's right_overhang must equal first part's left_overhang
        if parts[order[-1]]["right_overhang"] != parts[order[0]]["left_overhang"]:
            continue
        # Concatenate inserts, dropping the duplicated overhang at each junction.
        ohlen = _ENZYME_PROFILE[enzyme]["overhang"]
        product = parts[order[0]]["sequence"]
        for k in range(1, n):
            nxt = parts[order[k]]["sequence"]
            product += nxt[ohlen:]
        # Drop the trailing overhang (it duplicates the first overhang to close the circle)
        product = product[:-ohlen]
        # Register the assembled product
        product_name = "_".join(parts[i]["name"] for i in order) + "_GG"
        product_aid = reg.register_product(name=product_name[:60],
                                           sequence=product, circular=True)
        return {
            "product_attachment_id": product_aid,
            "length_bp": len(product),
            "enzyme": enzyme,
            "fragment_order": [parts[i]["name"] for i in order],
            "junction_overhangs": [parts[i]["right_overhang"] for i in order],
            "feasible": True,
        }

    return {
        "feasible": False,
        "reason": "no unique cyclic ligation order — overhangs ambiguous or incompatible",
        "fragments": [{"name": p["name"], "length_bp": p["length_bp"],
                       "left_overhang": p["left_overhang"],
                       "right_overhang": p["right_overhang"]}
                      for p in parts],
    }


async def tool_route_workflow(
    args, reg
):
    """Score every cloning workflow against an attached target + inventory
    using `target_from_inventory_router.route`. Returns the chosen winner
    plus the full report list (workflow, feasible, score, work_estimate,
    success_estimate, rationale) so the agent can pick the right method
    when the prompt does not name one explicitly."""
    target_id = args.get("target_attachment_id")
    inv_ids = args.get("inventory_attachment_ids") or []
    if not target_id:
        return {"error": "target_attachment_id is required"}
    target_att = reg.get(target_id)
    if not target_att:
        return {"error": f"unknown target_attachment_id: {target_id!r}"}
    inv_atts = []
    for aid in inv_ids:
        a = reg.get(aid)
        if a is None:
            return {"error": f"unknown inventory_attachment_id: {aid!r}"}
        inv_atts.append(a)
    try:
        from ..target_from_inventory_router import (
            annotate_one, route, build_audit_markdown,
        )
    except Exception as e:
        return {"error": f"router unavailable: {type(e).__name__}: {e}"}

    target_ctx = await annotate_one(target_att.name, None, target_att.sequence)
    inv_ctxs = [await annotate_one(a.name, None, a.sequence) for a in inv_atts]
    chosen, reports = route(target_ctx, inv_ctxs)

    summarized = []
    for r in reports:
        summarized.append({
            "workflow":         r.workflow,
            "feasible":         r.feasible,
            "score":            round(r.score, 3),
            "work_estimate":    r.work_estimate,
            "success_estimate": round(r.success_estimate, 3),
            "rationale":        (r.rationale or "")[:280],
        })
    summarized.sort(key=lambda x: (not x["feasible"], -x["score"]))
    return _strip_sequences({
        "chosen_workflow": (chosen.workflow if chosen else None),
        "chosen_score":    (round(chosen.score, 3) if chosen else None),
        "chosen_rationale": ((chosen.rationale or "")[:400] if chosen else None),
        "reports":         summarized,
        "audit_markdown":  build_audit_markdown(reports, chosen)[:1500],
    })


async def tool_lookup_kb_part(args, reg):
    """Search AIPlasmidDesign's knowledge base for a feature/part by name
    (CDS, fluorescent protein, common cloning element, etc.). Optionally,
    if `attachment_id` is supplied, align the KB sequence against that
    attachment and report what fraction of the KB part is present.
    """
    query = (args.get("name") or "").strip()
    if not query:
        return {"error": "name is required"}

    try:
        from ..predesign.knowledge_base import get_knowledge_base
    except Exception as e:
        return {"error": f"knowledge_base unavailable: {type(e).__name__}: {e}"}

    try:
        kb = get_knowledge_base()
    except Exception as e:
        return {"error": f"kb load failed: {type(e).__name__}: {e}"}
    try:
        hits = await kb.search_feature(query)
    except Exception as e:
        return {"error": f"kb search failed: {type(e).__name__}: {e}"}

    if not hits:
        return {"name": query, "matches": [], "note": "no KB matches; try a synonym (e.g. 'MS2 coat protein' for MCP)"}

    # Re-rank: prefer feature_motifs and feature_reference (curated GenoLIB
    # short motifs + canonical cloning parts) over swissprot/fpbase noise.
    # When the query is short (<= 5 chars, typical of motif acronyms like
    # P2A, T2A, NLS, FLAG), bias even harder toward motifs because
    # substring matches in protein DBs return false positives (e.g.
    # "Phosphatase 2A" for "P2A").
    _DB_PRIORITY = {"feature_motifs": 0, "feature_reference": 1,
                     "feature_protein": 2, "fpbase": 3, "swissprot": 4}
    _short_query = len(query) <= 5
    def _rank(h):
        db = h.get("database") or "other"
        base = _DB_PRIORITY.get(db, 9)
        if _short_query and db == "feature_motifs":
            base -= 100
        # Prefer exact-name matches.
        name = (h.get("name") or "").lower()
        exact = name == query.lower() or name.startswith(query.lower() + "_")
        return (base, 0 if exact else 1,
                -float(h.get("confidence") or 0.0), name)
    hits = sorted(hits, key=_rank)

    # Optionally align the top hit against an attachment to report fraction-present.
    attachment_id = args.get("attachment_id")
    fraction_payload = None
    if attachment_id:
        att = reg.get(attachment_id)
        if att is None:
            return {"error": f"unknown attachment_id: {attachment_id!r}"}
        top = hits[0]
        kb_seq = (top.get("sequence") or "").upper()
        kb_seq = "".join(c for c in kb_seq if c in "ACGTN")
        att_seq = att.sequence.upper()
        att_doubled = att_seq + (att_seq if att.circular else "")

        def _revcomp(s):
            comp = {"A":"T","T":"A","G":"C","C":"G","N":"N"}
            return "".join(comp.get(b, "N") for b in reversed(s))

        # Find longest contiguous KB-prefix that occurs anywhere on either strand of the attachment.
        # Use sliding 30-mer probes to locate likely region(s), then extend.
        def _largest_contig(kb_seq, hay):
            best = 0
            n = len(kb_seq)
            # Probe with a 30-nt window from the start, end, middle thirds.
            probe_len = min(30, n)
            for start in (0, n // 4, n // 2, 3 * n // 4, max(0, n - probe_len)):
                p = kb_seq[start:start + probe_len]
                if not p or p not in hay:
                    continue
                hay_idx = hay.find(p)
                # Extend left and right
                left = start
                hay_left = hay_idx
                while left > 0 and hay_left > 0 and kb_seq[left - 1] == hay[hay_left - 1]:
                    left -= 1; hay_left -= 1
                right = start + probe_len
                hay_right = hay_idx + probe_len
                while right < n and hay_right < len(hay) and kb_seq[right] == hay[hay_right]:
                    right += 1; hay_right += 1
                if right - left > best:
                    best = right - left
            return best

        rc = _revcomp(kb_seq)
        contig_fwd = _largest_contig(kb_seq, att_doubled)
        contig_rev = _largest_contig(rc, att_doubled)
        contig = max(contig_fwd, contig_rev)
        n = len(kb_seq) or 1
        frac = contig / n

        if frac >= 0.95:
            qualitative = "essentially the entire KB part is present"
        elif frac >= 0.55:
            qualitative = f"roughly {frac:.0%} of the KB part is present (about 2/3)" if 0.55 <= frac <= 0.78 else f"roughly {frac:.0%} of the KB part is present"
        elif frac >= 0.25:
            qualitative = f"roughly {frac:.0%} of the KB part is present (about 1/3)" if 0.25 <= frac <= 0.45 else f"roughly {frac:.0%} of the KB part is present"
        elif frac >= 0.05:
            qualitative = "only a small fragment / first few bases"
        else:
            qualitative = "essentially none of the KB part is present"

        fraction_payload = {
            "kb_part_length_bp": n,
            "longest_contiguous_match_bp": contig,
            "fraction_present": round(frac, 3),
            "best_strand": "forward" if contig_fwd >= contig_rev else "reverse_complement",
            "qualitative": qualitative,
        }

    # Resolve DNA for each hit. The KB stores some features as protein
    # only (e.g. mCherry, Cas9 from the GenoLIB protein tier); for those
    # we call feature_dna_resolver.get_feature_dna() which back-translates
    # via dnachisel and caches under feature_cds_cache.fna.
    register = bool(args.get("register_attachment", True))
    organism = (args.get("organism") or "h_sapiens").strip()
    try:
        from ..feature_dna_resolver import get_feature_dna
    except Exception:
        get_feature_dna = None  # type: ignore[assignment]

    def _looks_like_dna(s: str) -> bool:
        s = (s or "").upper()
        if len(s) < 30:
            return False
        nt = sum(1 for c in s if c in "ACGTN")
        return nt / max(len(s), 1) >= 0.95

    summarized = []
    for idx, h in enumerate(hits[:5]):
        seq = (h.get("sequence") or "")
        provenance = "direct" if _looks_like_dna(seq) else None
        dna = seq if provenance == "direct" else ""

        # Back-translate if we have a protein hit and the resolver is available.
        if not dna and get_feature_dna is not None:
            try:
                resolved = get_feature_dna(
                    h.get("id") or h.get("sseqid") or h.get("name") or "",
                    protein_sequence=seq if not _looks_like_dna(seq) else None,
                    json_representative_sequence=h.get("representative_sequence"),
                    organism=organism,
                    allow_backtranslation=True,
                )
                if resolved and resolved.sequence:
                    dna = resolved.sequence
                    provenance = resolved.provenance
            except Exception as e:
                logger.warning("kb back-translate failed for %r: %s", h.get("name"), e)

        row = {
            "id":          h.get("id"),
            "name":        h.get("name"),
            "database":    h.get("database"),
            "length_bp":   len(dna or seq),
            "confidence":  round(float(h.get("confidence") or 0.0), 3),
            "description": (h.get("description") or "")[:160],
            "provenance":  provenance,
        }

        # Register the TOP hit as an attachment so downstream tools
        # (graft_parts) can refer to it by attachment_id instead of
        # needing the LLM to thread a raw sequence back. Threshold is
        # low enough to capture short motifs (NLS=21 bp, FLAG=24 bp,
        # Kozak=8 bp) — anything shorter is usually a sub-motif we
        # don't want as a stand-alone part.
        if register and idx == 0 and dna and len(dna) >= 8:
            att_name = f"kb_{(h.get('name') or query).replace(' ', '_')[:40]}"
            aid = reg.register_product(name=att_name, sequence=dna,
                                        circular=False)
            row["attachment_id"] = aid

        summarized.append(row)

    out = {"query": query, "matches": summarized}
    if fraction_payload is not None:
        out["alignment_to_attachment"] = {"attachment_id": attachment_id, **fraction_payload}
    return _strip_sequences(out)


# --- dispatch table -------------------------------------------------------

TOOL_HANDLERS = {
    "annotate_attachment":       tool_annotate_attachment,
    "simulate_assembly":         tool_simulate_assembly,
    "digest_plasmid":            tool_digest_plasmid,
    "find_primer_binding_sites": tool_find_primer_binding,
    "score_sanger_primer":       tool_score_sanger_primer,
    "analyze_design_intent":     tool_analyze_design_intent,
    "verify_assembly":           tool_verify_assembly,
    "compare_to_choice":         tool_compare_to_choice,
    "golden_gate_assemble":      tool_golden_gate_assemble,
    "route_workflow":            tool_route_workflow,
    "lookup_kb_part":            tool_lookup_kb_part,
}


async def dispatch_tool(name: str, args: Dict[str, Any],
                        reg: AttachmentRegistry) -> Dict[str, Any]:
    handler = TOOL_HANDLERS.get(name)
    if handler is None:
        return {"error": f"unknown tool: {name}"}
    result = await handler(args, reg)
    return _strip_sequences(result)
