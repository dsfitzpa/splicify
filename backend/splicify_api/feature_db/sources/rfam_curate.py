"""Rfam curation — regex keep-list filter over the full Rfam 14.10 CM
file, then cmfetch+cmpress to produce an index suitable for
``cmscan --rfam --cut_ga`` in feature_annotator.py.

Spec categories (LLM_ANNOTATION_WORKFLOW.md, line 82) — viral cis-
acting elements, ribozymes, group I/II introns, riboswitches, aptamers
and fluorogenic RNAs, structural RNAs, bacterial regulatory sRNAs,
miRNA/lncRNA families, CRISPR repeats/cr/tracrRNA, att sites, UTR cis
elements (SECIS / Histone 3' / IRE / polyA).

Pipeline:
    external_src/Rfam.cm                    (full, 8,340 families)
    external_src/Rfam.clanin                (full clan map)
    -> parse NAME/DESC for each family
    -> match against KEEP_PATTERNS
    -> emit keep_ids.txt
    -> cmfetch -f -> rfam.cm
    -> filter clanin to clans with >=2 remaining members
       -> rfam.clanin
    -> cmpress rfam.cm
    -> rm Rfam.cm Rfam.cm.gz Rfam.clanin

Run:
    python -m splicify_api.feature_db.sources.rfam_curate
"""
from __future__ import annotations

import re
import shutil
import subprocess
import sys
from pathlib import Path

_SCRIPT_DIR = Path(__file__).resolve().parents[1]
_SPLICIFY_API = _SCRIPT_DIR.parent
_FDD = _SPLICIFY_API / "feature_db_data"
_EXT = _FDD / "external_src"
_BLAST = _FDD / "BLAST_dbs"


# ---------------------------------------------------------------------------
# Keep-list patterns. Case-insensitive regex; matches NAME or DESC.
# Categories follow the spec wording verbatim where practical.
# ---------------------------------------------------------------------------
KEEP_PATTERNS: list[str] = [
    # ---- Viral cis-acting elements ------------------------------------------
    r"\bHIV.{0,3}(psi|RRE|TAR|FE|SLI{1,3}|DIS|SL[0-9]+)\b",
    r"\bWPRE\b", r"\bcPPT\b", r"\bCTS\b",
    r"\bHCV\b", r"\bEMCV\b", r"\bFMDV\b", r"\bPolio\b", r"\bCrPV\b",
    r"\bIRES\b",
    r"\bFlavi", r"\bCorona", r"\bEntero",
    r"\bCRE\b.*virus", r"virus.*\bCRE\b",
    r"\bHEV\b", r"\bSARS\b", r"\b3.UTR.*virus", r"\bRSV\b",
    r"\bVirus_5", r"\bVirus_3",
    # ---- Ribozymes, group I/II introns --------------------------------------
    r"ribozyme", r"HDV", r"hammerhead", r"twister",
    r"\bGroup[_ ]?[III]+", r"intron.*group",
    r"hairpin.{0,3}ribozyme", r"\bVS_ribo",
    # ---- Riboswitches -------------------------------------------------------
    r"\bTPP\b", r"\bFMN\b", r"\bSAM\b", r"purine", r"guanine",
    r"\bB12\b", r"cobalamin", r"fluoride", r"\bMg\b", r"c.di.GMP",
    r"riboswitch",
    r"glycine.*riboswitch", r"lysine.*riboswitch", r"theo",
    # ---- Aptamers + fluorogenic RNAs ----------------------------------------
    r"Spinach", r"Broccoli", r"Mango", r"Pepper",
    r"theophylline", r"tetracycline", r"neomycin", r"streptomycin",
    r"aptamer",
    # ---- Structural RNAs ----------------------------------------------------
    r"\btRNA\b",
    r"\b5S_rRNA\b", r"\b5\.8S_rRNA\b", r"\b18S_rRNA\b", r"\b28S_rRNA\b",
    r"\bSSU_rRNA", r"\bLSU_rRNA",
    r"\b7SK\b", r"\b7SL\b",
    r"\bU[0-9]+\b",                  # U1-U12, U6atac
    r"\bU[0-9]+atac\b",
    r"\bVault\b",
    r"\bY[_ ]?RNA\b",
    r"RNase.?P", r"RNase.?MRP",
    r"telomerase", r"\btmRNA\b", r"\bSRP\b",
    # ---- Bacterial regulatory sRNAs -----------------------------------------
    r"\bDsrA\b", r"\bRyhB\b", r"\bMic[A-F]\b", r"\bOxyS\b",
    r"\bCsrB\b", r"\bCsrC\b", r"\b6S\b", r"Spot.?42",
    r"\bCopA\b", r"\bDicF\b", r"\bPrfA\b", r"\bQrr\b", r"T[_ ]box",
    # ---- miRNA + lncRNA -----------------------------------------------------
    r"\blet-7", r"\bmir[-_]", r"\bmiR[-_]",
    r"\bXist\b", r"\bMALAT1\b", r"\bNEAT1\b", r"\bHOTAIR\b",
    r"long.?non.?coding", r"\blincRNA\b",
    # ---- att sites ----------------------------------------------------------
    r"\batt[BPLR][1-9]?\b", r"\battB\b", r"\battP\b", r"\battL\b", r"\battR\b",
    # ---- CRISPR repeats / cr / tracrRNA -------------------------------------
    r"CRISPR.*repeat", r"\bcrRNA\b", r"\btracrRNA\b",
    # ---- UTR cis elements ---------------------------------------------------
    r"SECIS", r"Histone.?3.UTR", r"Histone.?stem.?loop",
    r"\bIRE\b", r"iron.responsive",
    r"polyA", r"poly.?A.signal",
    r"\bAUF1\b", r"\bAU.rich\b",
    r"\bK.turn", r"\bC.D.box\b", r"\bH.ACA",
]
_KEEP_RE = re.compile("|".join(KEEP_PATTERNS), re.IGNORECASE)


def parse_cm_families(cm_path: Path) -> list[dict]:
    """Stream Rfam.cm and emit one dict per family with NAME / ACC / DESC.

    Only reads header lines; sequence emission models are skipped.
    """
    families: list[dict] = []
    cur: dict[str, str] = {}
    with cm_path.open("r") as fh:
        for line in fh:
            if line.startswith("INFERNAL"):
                if cur:
                    families.append(cur)
                cur = {}
            elif line.startswith("NAME "):
                cur["name"] = line[5:].strip()
            elif line.startswith("ACC "):
                cur["acc"] = line[4:].strip()
            elif line.startswith("DESC "):
                cur["desc"] = line[5:].strip()
            elif line.startswith("CM"):
                # CM section header — we have everything we need
                if cur:
                    families.append(cur)
                cur = {}
    return families


def _is_keep(fam: dict) -> bool:
    text = f"{fam.get('name','')} {fam.get('desc','')}"
    return bool(_KEEP_RE.search(text))


def filter_clanin(src_clanin: Path, keep_names: set[str], out_clanin: Path) -> int:
    """Drop clan members not in keep_names; keep only clans with >=2
    remaining members (matches spec)."""
    kept_clans = 0
    with src_clanin.open("r") as src, out_clanin.open("w") as dst:
        for line in src:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            clan_id, *members = parts
            kept_members = [m for m in members if m in keep_names]
            if len(kept_members) >= 2:
                dst.write("\t".join([clan_id] + kept_members) + "\n")
                kept_clans += 1
    return kept_clans


def main() -> int:
    src_cm = _EXT / "Rfam.cm"
    src_clanin = _EXT / "Rfam.clanin"
    if not src_cm.exists() or not src_clanin.exists():
        print(f"FAIL: Rfam.cm or Rfam.clanin missing in {_EXT}", file=sys.stderr)
        return 1

    print(f"[rfam] parsing {src_cm} (this takes a few s)")
    families = parse_cm_families(src_cm)
    print(f"[rfam] {len(families)} families in source")

    keep = [f for f in families if _is_keep(f)]
    keep_names = {f["name"] for f in keep}
    print(f"[rfam] keep-list match: {len(keep)} families")

    # Write keep-ids file (one NAME per line) for cmfetch.
    keep_ids_path = _FDD / "rfam_keep_ids.txt"
    keep_ids_path.write_text("\n".join(sorted(keep_names)) + "\n")
    print(f"[rfam] {keep_ids_path} -> {len(keep_names)} names")

    # cmfetch -f writes a multi-CM file from a list of names.
    _BLAST.mkdir(exist_ok=True)
    out_cm = _BLAST / "rfam.cm"
    if out_cm.exists():
        out_cm.unlink()
    for suffix in (".i1f", ".i1i", ".i1m", ".i1p"):
        p = out_cm.with_suffix(out_cm.suffix + suffix)
        if p.exists():
            p.unlink()

    print(f"[rfam] cmfetch -> {out_cm}")
    rc = subprocess.run(
        ["cmfetch", "-f", "-o", str(out_cm), str(src_cm), str(keep_ids_path)],
        capture_output=True, text=True,
    )
    if rc.returncode != 0:
        print(f"[rfam] cmfetch FAIL: {rc.stderr.strip()[:500]}")
        return 1
    print(f"[rfam] cmfetch OK ({out_cm.stat().st_size // (1024*1024)} MB)")

    # cmpress builds .i1f .i1i .i1m .i1p indices
    print(f"[rfam] cmpress {out_cm}")
    rc = subprocess.run(
        ["cmpress", "-F", str(out_cm)],
        capture_output=True, text=True,
    )
    if rc.returncode != 0:
        print(f"[rfam] cmpress FAIL: {rc.stderr.strip()[:500]}")
        return 1
    print(f"[rfam] cmpress OK")

    out_clanin = _BLAST / "rfam.clanin"
    kept_clans = filter_clanin(src_clanin, keep_names, out_clanin)
    print(f"[rfam] rfam.clanin: {kept_clans} clans (>=2 remaining members)")

    # Delete originals (user-requested)
    print(f"[rfam] deleting originals")
    for p in (src_cm, _EXT / "Rfam.cm.gz", src_clanin):
        if p.exists():
            sz_mb = p.stat().st_size // (1024 * 1024)
            p.unlink()
            print(f"[rfam]   removed {p.name} ({sz_mb} MB)")

    print(f"[rfam] DONE — curated to {len(keep_names)} families")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
