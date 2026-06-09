"""CRISPR + primer design tool wrappers for the agent_v2 main agent.

Layered on top of v1's `guide_designer.design_guides`,
`pegrna_designer.design_pegrnas`, and `pcr.design_primers`. Each wrapper:

  - Resolves attachment_id -> registry.get(aid).sequence
  - Forwards the call to v1 with sensible defaults
  - Returns a digested envelope that excludes raw template DNA but
    preserves the actual designed sequences (spacers, primers, full pegRNAs)
    since those ARE the answer the agent needs to reply with.

Dispatched through the same `dispatch_with_emitters` chain used by the
output emitters, so the main agent gets one uniform handler table.
"""
from __future__ import annotations
import asyncio
import base64

from typing import Any, Optional

import agent_v2  # noqa: F401 - triggers path shim


# Max wall-clock time for a single primer3 call. primer3 is a synchronous
# C extension — if it stalls on an unsatisfiable region (rare but real),
# calling it on the event-loop thread freezes the entire orchestrator.
# Wrapping the call in asyncio.to_thread + asyncio.wait_for fixes both
# problems: the sync call runs in a worker thread (event loop stays
# responsive), and the timeout caps each request so the LLM gets a clean
# ok=False envelope instead of a hung tool spinner that systemd has to
# SIGKILL.
PRIMER3_TIMEOUT_SECONDS = 20.0


# Sonnet 4.6 sometimes passes the entire response of an earlier tool call
# (e.g. resolve_feature_position -> {plasmid_position, codon, ...}) as the
# argument to a later tool that expects a scalar int. Plain int(v) raises
# TypeError on dicts / lists. This helper extracts a sensible integer from
# the common bad shapes so the pipeline survives the LLM mistake.
_INT_KEYS = ("plasmid_position", "position", "start", "value", "left_pos", "right_pos")


def _arg_to_int(v, default=None):
    """Best-effort int conversion that tolerates the shapes Claude actually emits.

    - int / float / bool       -> int(v)
    - numeric str              -> int(v)
    - tuple / list (non-empty) -> recursively convert first element
    - dict                     -> the first int-like value at a known key
    - None / unrecoverable     -> default (may be None)
    """
    if v is None:
        return default
    if isinstance(v, bool):
        return int(v)
    if isinstance(v, (int, float)):
        return int(v)
    if isinstance(v, str):
        s = v.strip()
        if not s:
            return default
        try:
            return int(s)
        except (TypeError, ValueError):
            try:
                return int(float(s))
            except (TypeError, ValueError):
                return default
    if isinstance(v, (list, tuple)):
        if not v:
            return default
        return _arg_to_int(v[0], default=default)
    if isinstance(v, dict):
        for k in _INT_KEYS:
            if k in v:
                inner = _arg_to_int(v[k], default=None)
                if inner is not None:
                    return inner
        # Fall back to any single int-typed value inside the dict.
        for val in v.values():
            if isinstance(val, (int, float)) and not isinstance(val, bool):
                return int(val)
        return default
    return default


# ---------------------------------------------------------------------------
# Cas9 / Cas12a sgRNA design
# ---------------------------------------------------------------------------
async def design_guides_tool(
    args: dict[str, Any],
    registry: Any,
    *,
    output_dir: Optional[str] = None,  # unused; kept for dispatch signature
) -> dict[str, Any]:
    aid = args.get("attachment_id")
    att = registry.get(aid) if aid else None
    if att is None:
        return {"ok": False, "error": f"unknown attachment_id: {aid!r}"}

    L = len(att.sequence)
    rs = _arg_to_int(args.get("region_start"), default=1)
    re_ = _arg_to_int(args.get("region_end"), default=L)
    if rs is None or re_ is None or rs < 1 or re_ < rs or re_ > L:
        return {"ok": False, "error": f"invalid region {rs!r}..{re_!r} for {L}-bp plasmid"}

    from splicify_api.guide_designer import design_guides as _v1
    raw = _v1(
        sequence=att.sequence,
        region_start=rs,
        region_end=re_,
        pam=args.get("pam", "NGG"),
        guide_length=_arg_to_int(args.get("guide_length"), default=20),
        pam_position=args.get("pam_position", "3prime"),
        max_guides=_arg_to_int(args.get("max_guides"), default=5),
        min_score=float(args.get("min_score") or 0.0),
        score_method=args.get("score_method", "doench2014"),
    )
    if isinstance(raw, dict) and not raw.get("ok", True):
        raw.setdefault("attachment_id", aid)
        return raw

    digested = []
    for g in raw.get("guides", []) or []:
        digested.append({
            "name": g.get("name"),
            "spacer": g.get("spacer"),
            "pam": g.get("pam"),
            "start": g.get("start"),
            "end": g.get("end"),
            "direction": g.get("direction"),
            "score": g.get("score"),
            "score_method": g.get("score_method"),
            "gc_fraction": g.get("gc_fraction"),
            "n_offtargets": g.get("n_offtargets"),
        })
    return {
        "ok": True,
        "attachment_id": aid,
        "guides": digested,
        "summary": raw.get("summary", {}),
    }


# ---------------------------------------------------------------------------
# pegRNA design (prime editing)
# ---------------------------------------------------------------------------
async def design_pegrnas_tool(
    args: dict[str, Any],
    registry: Any,
    *,
    output_dir: Optional[str] = None,
) -> dict[str, Any]:
    aid = args.get("attachment_id")
    att = registry.get(aid) if aid else None
    if att is None:
        return {"ok": False, "error": f"unknown attachment_id: {aid!r}"}

    edit_start = _arg_to_int(args.get("edit_start"))
    edit_end = _arg_to_int(args.get("edit_end"))
    if edit_start is None or edit_end is None:
        return {"ok": False,
                 "error": ("edit_start and edit_end required as 1-indexed plasmid "
                            "integers. If you previously called resolve_feature_position, "
                            "pass its `plasmid_position` field — not the whole result object.")}

    from splicify_api.pegrna_designer import design_pegrnas as _v1
    raw = _v1(
        sequence=att.sequence,
        edit_start_1based=edit_start,
        edit_end_1based=edit_end,
        alt=args.get("alt", "") or "",
        edit_type=args.get("edit_type", "substitution"),
        n_results=_arg_to_int(args.get("n_results"), default=3),
        use_pe3=bool(args.get("use_pe3", True)),
    )
    if isinstance(raw, dict) and not raw.get("ok", True):
        raw.setdefault("attachment_id", aid)
        return raw

    digested = []
    for p in raw.get("pegrnas", []) or []:
        digested.append({
            "rank": p.get("rank"),
            "name": p.get("name"),
            "predicted_efficiency": p.get("predicted_efficiency"),
            "spacer": p.get("spacer"),
            "pam": p.get("pam"),
            "spacer_start": p.get("spacer_start"),
            "spacer_end": p.get("spacer_end"),
            "direction": p.get("direction"),
            "cas9_score": p.get("cas9_score"),
            "rtt": p.get("rtt"),
            "rtt_length": p.get("rtt_length"),
            "pbs": p.get("pbs"),
            "pbs_length": p.get("pbs_length"),
            "scaffold": p.get("scaffold"),
            "full_pegrna": p.get("full_pegrna"),
            "full_pegrna_length": p.get("full_pegrna_length"),
            "is_dpam": p.get("is_dpam"),
            "is_pe3b": p.get("is_pe3b"),
            "ngrna": p.get("ngrna"),
            "edit_type": p.get("edit_type"),
            "edit_ref": p.get("edit_ref"),
            "edit_alt": p.get("edit_alt"),
            "edit_start_1based": p.get("edit_start_1based"),
            "edit_end_1based": p.get("edit_end_1based"),
        })
    return {
        "ok": True,
        "attachment_id": aid,
        "pegrnas": digested,
        "summary": raw.get("summary", {}),
    }


# ---------------------------------------------------------------------------
# Primer design (fragment / sanger / illumina)
# ---------------------------------------------------------------------------
async def design_primers_tool(
    args: dict[str, Any],
    registry: Any,
    *,
    output_dir: Optional[str] = None,
) -> dict[str, Any]:
    aid = args.get("attachment_id")
    att = registry.get(aid) if aid else None
    if att is None:
        return {"ok": False, "error": f"unknown attachment_id: {aid!r}"}

    L = len(att.sequence)
    rs = _arg_to_int(args.get("region_start"), default=1)
    re_ = _arg_to_int(args.get("region_end"), default=L)
    if rs is None or re_ is None or rs < 1 or re_ < rs or re_ > L:
        return {"ok": False, "error": f"invalid region {rs!r}..{re_!r} for {L}-bp plasmid"}

    template = att.sequence[rs - 1:re_]
    template_len = len(template)

    excl_start_plasmid = _arg_to_int(args.get("excluded_start"))
    excl_end_plasmid = _arg_to_int(args.get("excluded_end"))
    excluded_start_tpl: Optional[int] = None
    excluded_length: Optional[int] = None
    if excl_start_plasmid is not None and excl_end_plasmid is not None:
        excluded_start_tpl = excl_start_plasmid - rs
        excluded_length = excl_end_plasmid - excl_start_plasmid + 1
        if excluded_start_tpl < 0 or excluded_start_tpl + excluded_length > template_len:
            return {"ok": False, "error": "excluded region outside the requested plasmid region"}

    application = args.get("application", "sanger")
    default_min, default_max = {
        "fragment": (100, 300),
        "sanger":   (250, 500),
        "illumina": (150, 290),
    }.get(application, (200, 500))

    from splicify_api.pcr import PrimerRequest, design_primers as _v1
    req = PrimerRequest(
        fragments_in=template,
        excluded_start=excluded_start_tpl,
        excluded_length=excluded_length,
        product_size_min=_arg_to_int(args.get("product_size_min"), default=default_min),
        product_size_max=_arg_to_int(args.get("product_size_max"), default=default_max),
        primer_opt_tm=float(args.get("primer_opt_tm") or 60.0),
        num_return=_arg_to_int(args.get("num_return"), default=5),
        application=application,
    )
    # Run primer3 in a worker thread so the sync C call cannot block the
    # event loop, and cap it with a wall-clock timeout so a hung call
    # fails fast with a clean envelope instead of stalling the entire
    # orchestrator (and forcing systemd to SIGKILL on shutdown).
    try:
        raw = await asyncio.wait_for(
            asyncio.to_thread(_v1, req),
            timeout=PRIMER3_TIMEOUT_SECONDS,
        )
    except asyncio.TimeoutError:
        return {
            "ok": False,
            "attachment_id": aid,
            "region_start_plasmid": rs,
            "region_end_plasmid": re_,
            "application": application,
            "error": (
                f"primer3 timed out after {PRIMER3_TIMEOUT_SECONDS:.0f}s "
                f"on region {rs}..{re_} ({template_len} bp, excluded "
                f"{excluded_start_tpl}+{excluded_length}). The constraints "
                f"are likely unsatisfiable — try widening region_start/"
                f"region_end, loosening product_size_min/max, or removing "
                f"the excluded region."
            ),
        }
    except Exception as e:
        return {
            "ok": False,
            "attachment_id": aid,
            "region_start_plasmid": rs,
            "region_end_plasmid": re_,
            "application": application,
            "error": f"primer3 raised {type(e).__name__}: {e}",
        }

    def _pos_to_int(v):
        # primer3 returns left_pos / right_pos as (start, length) tuples
        # (serialised as 2-element lists across JSON). Accept either shape:
        # tuple/list -> first element; bare int -> use as-is. Anything else
        # -> None (defensive; we never want to int() a dict).
        if v is None:
            return None
        if isinstance(v, (list, tuple)) and v:
            try:
                return int(v[0])
            except (TypeError, ValueError):
                return None
        if isinstance(v, (int, float)):
            return int(v)
        return None

    if isinstance(raw, dict):
        if not raw.get("ok", True):
            raw.setdefault("attachment_id", aid)
            return raw
        lp = _pos_to_int(raw.get("left_pos"))
        rp = _pos_to_int(raw.get("right_pos"))
        # 1-indexed plasmid region coords; primer3 returns template-relative
        # 0-indexed positions. plasmid_pos_0idx = (rs - 1) + template_pos.
        if lp is not None:
            raw["left_pos_plasmid"] = lp + rs - 1
        if rp is not None:
            raw["right_pos_plasmid"] = rp + rs - 1
        raw["attachment_id"] = aid
        raw["region_start_plasmid"] = rs
        raw["region_end_plasmid"] = re_
    return raw


# ---------------------------------------------------------------------------
# External genomic-record lookup (NCBI Gene -> RefSeqGene .gb -> register)
# ---------------------------------------------------------------------------
async def find_genomic_record_tool(
    args: dict[str, Any],
    registry: Any,
    *,
    output_dir: Optional[str] = None,
) -> dict[str, Any]:
    """Resolve a gene symbol + organism to a registered genomic .gb attachment.

    Designed for the CRISPR / pegRNA pipeline when the user references a
    gene the agent doesn't have a local .gb for. Downloads the RefSeqGene
    (NG_*) when available — that's the genomic record with introns +
    flanking the pegRNA designer needs — and falls back to RefSeq mRNA
    (NM_*) when no RefSeqGene exists.
    """
    gene_symbol = (args.get("gene_symbol") or "").strip()
    organism = (args.get("organism") or "Homo sapiens").strip() or "Homo sapiens"
    if not gene_symbol:
        return {"ok": False, "error": "gene_symbol required (e.g. 'CGAS')"}

    from splicify_api.external_search import search_ncbi_gene, fetch_ncbi_genbank
    from agent_v2.file_kind import classify_genbank
    from agent_v2.genomic_annotator import annotate_genomic_gb
    from agent_v2 import attachment_kinds

    # Belt-and-suspenders: search_ncbi_gene + fetch_ncbi_genbank both already
    # catch httpx.HTTPError and return None (see PR #80
    # fix/ncbi-rate-limit-graceful). Wrap here too so ANY uncaught exception
    # (network blip, asyncio cancellation, unexpected upstream shape)
    # becomes a clean ok=False envelope rather than a cryptic
    # "ConnectTimeout: " bubbling up to the orchestrator's catch-all and
    # appearing in the SSE error event.
    try:
        hit = await search_ncbi_gene(gene_symbol, organism)
    except Exception as e:
        return {"ok": False,
                 "error": (f"NCBI gene search failed: {type(e).__name__}: {e}. "
                            "NCBI may be rate-limiting; retry in a few seconds, "
                            "or upload the .gb manually."),
                 "gene_symbol": gene_symbol, "organism": organism}
    if hit is None:
        return {
            "ok": False,
            "error": (f"no NCBI gene record found for {gene_symbol!r} in {organism!r}. "
                       "Either the gene symbol / organism is wrong, or NCBI returned "
                       "no hits / rate-limited the lookup. Retry, or upload the .gb manually."),
            "gene_symbol": gene_symbol, "organism": organism,
        }

    try:
        gb_text = await fetch_ncbi_genbank(
            hit.accession,
            seq_start=hit.seq_start, seq_stop=hit.seq_stop,
        )
    except Exception as e:
        return {"ok": False,
                 "error": (f"NCBI .gb fetch failed for {hit.accession!r}: "
                            f"{type(e).__name__}: {e}. NCBI may be rate-limiting; "
                            "retry in a few seconds, or upload the .gb manually."),
                 "gene_symbol": gene_symbol, "organism": organism,
                 "accession": hit.accession}
    if not gb_text:
        return {
            "ok": False,
            "error": (f"NCBI returned a hit ({hit.accession}) but the .gb download "
                       "failed or was empty. Ask the user to upload the .gb manually."),
            "gene_symbol": gene_symbol, "organism": organism,
            "accession": hit.accession, "title": hit.title,
        }

    # Pull the sequence + register as a non-circular genomic attachment.
    try:
        from splicify_api.agent.agent_tools import extract_seq_from_genbank
        seq = extract_seq_from_genbank(gb_text)
    except Exception as e:
        return {"ok": False,
                 "error": f"sequence extraction failed: {type(e).__name__}: {e}",
                 "accession": hit.accession}
    if not seq or len(seq) < 50:
        return {"ok": False,
                 "error": f"NCBI .gb for {hit.accession} had no usable sequence",
                 "accession": hit.accession}

    classification = classify_genbank(gb_text)
    aid = registry.register_product(
        hit.accession or gene_symbol, seq,
        circular=False,   # always linear for NCBI gene records
    )
    attachment_kinds.stash_kind(aid, classification, gb_text=gb_text)

    # Eagerly build + stash the GenomicAnnotation so resolve_feature_position
    # doesn't re-annotate on its first call. The digest below reports
    # transcripts + CDS protein length to the LLM in this same response.
    try:
        ann = annotate_genomic_gb(gb_text)
        attachment_kinds.set_annotation(aid, ann)
        gene_keys = list(ann.genes.keys())
        transcript_ids = list(ann.transcripts.keys())
        cds_count = sum(1 for f in ann.features if f.type == "CDS")
        protein_lens = [len(f.translation or "") for f in ann.features
                          if f.type == "CDS" and f.translation]
        primary_cds_aa = max(protein_lens) if protein_lens else 0
    except Exception:
        gene_keys, transcript_ids, cds_count, primary_cds_aa = [], [], 0, 0

    # Surface the retrieved .gb to the frontend mid-stream so the user
    # sees the gene visualisation as soon as it's fetched, instead of
    # waiting for the full pegRNA-design pipeline to finish. The file
    # envelope rides along on the tool result and is forwarded via the
    # on_tool_event hook in main_agent.
    safe_name = (gene_symbol or hit.accession or "genomic").strip().replace(" ", "_")
    safe_name = "".join(c for c in safe_name if c.isalnum() or c in "._-")[:60]
    file_envelope = {
        "fileName": f"{safe_name}.gb" if safe_name else "genomic.gb",
        "dataBase64": base64.b64encode(gb_text.encode("utf-8")).decode("ascii"),
    }

    notes = {
        "refseqgene": "RefSeqGene record — genomic with introns + flanking.",
        "chromosomal_slice": (
            "Chromosomal slice from the gene's primary-assembly placement "
            "(NC_*). Includes introns + exons + the user-configured flanking "
            "window on each side, so intron/exon-boundary guides design "
            "correctly. Coordinates returned by resolve_feature_position "
            "and design_* tools are on this slice (1-indexed within seq_start "
            "to seq_stop), not on the parent chromosome."
        ),
        "mrna": (
            "mRNA fallback — the record is spliced (no introns), so intronic "
            "or splice-junction edits are out of scope on this attachment."
        ),
    }
    return {
        "ok": True,
        "file": file_envelope,
        "attachment_id": aid,
        "source": "ncbi",
        "db_source": hit.db_source,     # "refseqgene" | "chromosomal_slice" | "mrna"
        "gene_symbol": gene_symbol,
        "organism": organism,
        "accession": hit.accession,
        "ncbi_gene_id": hit.gene_id,
        "title": hit.title,
        "length_bp": len(seq),
        "kind": classification.kind,
        "kind_confidence": classification.confidence,
        "genes": gene_keys[:10],
        "transcripts": transcript_ids[:10],
        "n_cds": cds_count,
        "primary_cds_aa": primary_cds_aa,
        # chromosomal-slice spatial context (None for the other two db_source values)
        "chromosome_accession": hit.accession if hit.db_source == "chromosomal_slice" else None,
        "slice_seq_start": hit.seq_start,
        "slice_seq_stop": hit.seq_stop,
        "flanking_bp": hit.flanking_bp,
        "gene_chr_start": hit.gene_chr_start,
        "gene_chr_stop": hit.gene_chr_stop,
        "note": notes.get(hit.db_source, ""),
    }
