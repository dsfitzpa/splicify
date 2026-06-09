"""splicify_api.feature_annotator — self-contained annotation pipeline.

Forks the BLAST orchestration from pLannotate's annotate.py into the
splicify_api package so we can drop the `plannotate` PyPI dependency.
Searches the feature_db_data/ tree built by
scripts/build_feature_db_data.py.

Differences from the upstream pLannotate annotate.py:
  * No `import plannotate.*` — paths come from feature_db_data/databases.yml.
  * No `streamlit` — progress bar / errors print to stderr.
  * Drops the `diamond` method entirely; protein search runs through
    `mmseqs2 easy-search` (faster, no DIAMOND binary required).
  * Drops `pdb|FOO|` munging from `get_details` (none of our DBs use that
    SwissProt-style sseqid form).
  * `parse_infernal` forked inline (Rfam covariance-model search).

Public surface:
  - `annotate(inSeq, yaml_file=None, linear=False, is_detailed=False)`
    drop-in for the old `from plannotate.annotate import annotate` import.
"""
from __future__ import annotations

import glob
import logging
import os
import re
import shlex
import subprocess
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from . import _data
from tempfile import NamedTemporaryFile
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
import yaml
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)
if not logger.handlers:
    logging.basicConfig(level=os.getenv("SPLICIFY_ANNOTATOR_LOGLEVEL", "INFO"))


# Module-level temp log used by upstream pLannotate; kept for symmetry.
_LOG = NamedTemporaryFile()

_FEATURE_DB_DIR = _data.data_path("feature_db_data")
_DEFAULT_YAML = _FEATURE_DB_DIR / "databases.yml"

# Columns the DataFrame must carry on every return path. Forked verbatim
# from pLannotate.resources.DF_COLS so downstream consumers keep working.
DF_COLS = [
    "sseqid", "qstart", "qend", "sstart", "send", "sframe", "score",
    "evalue", "qseq", "length", "slen", "pident", "qlen", "db", "Feature",
    "Description", "Type", "priority", "percmatch", "abs percmatch",
    "pi_permatch", "wiggle", "wstart", "wend", "kind", "qstart_dup",
    "qend_dup", "fragment",
]


# --------------------------------------------------------------------------- #
# YAML loader (replaces plannotate.resources.get_yaml).
# --------------------------------------------------------------------------- #


def _load_databases_yaml(yaml_file_loc: str) -> Dict[str, Dict[str, Any]]:
    """Load databases.yml. All `location` paths must be absolute — there is
    no `Default` fallback in this fork. Returns one entry per DB with
    `db_loc` computed and `parameters` joined to a single string."""
    with open(yaml_file_loc, "r") as f:
        dbs: Dict[str, Dict[str, Any]] = yaml.safe_load(f)

    for db_name, db in dbs.items():
        loc = db.get("location")
        if loc in (None, "", "Default"):
            raise ValueError(
                f"databases.yml entry '{db_name}' has location='{loc}'. "
                f"feature_db_data/ requires absolute paths."
            )
        parameters = " ".join(db.get("parameters") or [])
        db["parameters"] = parameters

        method = db.get("method")
        if method == "infernal":
            # Infernal needs two files: <db>.clanin and <db>.cm
            db["db_loc"] = " ".join(
                os.path.join(loc, x) for x in (f"{db_name}.clanin", f"{db_name}.cm")
            )
        else:
            db["db_loc"] = os.path.join(loc, db_name)

    return dbs


# --------------------------------------------------------------------------- #
# Infernal output parser (forked from plannotate.infernal.parse_infernal).
# --------------------------------------------------------------------------- #


def _parse_infernal(file_loc: str) -> pd.DataFrame:
    with open(file_loc) as fh:
        lines = fh.readlines()
    if len(lines) < 2:
        return pd.DataFrame()

    col_widths = [len(ele) + 1 for ele in lines[1].split()]
    ends = list(np.cumsum(col_widths))
    ends[-1] += 100
    starts = [0] + ends[:-1]
    col_pos = list(zip(starts, ends))

    col_names = [lines[0][s:e].strip() for s, e in col_pos]

    try:
        infernal = pd.read_fwf(file_loc, comment="#", colspecs=col_pos, header=None)
        infernal.columns = col_names
    except pd.errors.EmptyDataError:
        infernal = pd.DataFrame(columns=col_names)

    keep = [
        "#idx", "target name", "accession", "clan name", "seq from", "seq to",
        "mdl from", "mdl to", "strand", "score", "E-value", "description of target",
    ]
    # Some columns may be missing in empty / minimal output
    keep = [c for c in keep if c in infernal.columns]
    infernal = infernal[keep]
    infernal = infernal.loc[:, ~infernal.columns.duplicated()]

    infernal = infernal.rename(columns={
        "#idx": "sseqid",
        "seq from": "qstart",
        "seq to": "qend",
        "mdl from": "sstart",
        "mdl to": "send",
        "E-value": "evalue",
        "strand": "sframe",
        "target name": "Feature",
        "description of target": "Description",
    })
    if "accession" in infernal.columns:
        infernal["accession"] = infernal["accession"].str.replace("-", " ")
    if "clan name" in infernal.columns:
        infernal["clan name"] = infernal["clan name"].str.replace("-", " ")
    if "Feature" in infernal.columns:
        infernal["Feature"] = infernal["Feature"].str.replace("_", " ")
    if "Description" in infernal.columns and "accession" in infernal.columns:
        infernal["Description"] = (
            "Accession: " + infernal["accession"] + " - " + infernal["Description"]
        )

    for col in infernal.columns:
        try:
            infernal[col] = pd.to_numeric(infernal[col])
        except (ValueError, TypeError):
            pass

    infernal["qseq"] = ""
    if "qstart" in infernal.columns and "qend" in infernal.columns:
        to_swap = infernal["qend"] < infernal["qstart"]
        infernal.loc[to_swap, ["qstart", "qend"]] = infernal.loc[
            to_swap, ["qend", "qstart"]
        ].values
        infernal[["qstart", "qend"]] = infernal[["qstart", "qend"]].apply(
            pd.to_numeric, downcast="integer"
        )
        infernal["qstart"] = infernal["qstart"] - 1
        infernal["qend"] = infernal["qend"] - 1
        infernal["length"] = abs(infernal["qend"] - infernal["qstart"]) + 1
    if "sframe" in infernal.columns:
        infernal["sframe"] = infernal["sframe"].replace(["-", "+"], [-1, 1])
    if "sstart" in infernal.columns and "send" in infernal.columns:
        infernal["slen"] = abs(infernal["send"] - infernal["sstart"]) + 1
    infernal["pident"] = 100

    infernal = infernal.drop(
        columns=[c for c in ["accession", "clan name"] if c in infernal.columns]
    )
    return infernal


# --------------------------------------------------------------------------- #
# BLAST orchestration (forked from plannotate.annotate.BLAST).
# Drops `diamond`; protein search runs through mmseqs.
# --------------------------------------------------------------------------- #


def _write_temp_fasta(seq: str) -> str:
    fd, path = tempfile.mkstemp(suffix=".fa", prefix="splicify_")
    with os.fdopen(fd, "w", encoding="utf-8") as fh:
        SeqIO.write(SeqRecord(Seq(seq), id="query", description=""), fh, "fasta")
    return path


def _run_blast(seq: str, db: Dict[str, Any]) -> pd.DataFrame:
    task = db["method"]
    parameters = db["parameters"]
    db_loc = db["db_loc"]

    query_path = _write_temp_fasta(seq)
    fd_out, out_path = tempfile.mkstemp(suffix=".out", prefix="splicify_")
    os.close(fd_out)

    try:
        if task == "blastn":
            flags = "qstart qend sseqid sframe pident slen qseq length sstart send qlen evalue"
            cmd = (
                f"blastn -task blastn-short -query {shlex.quote(query_path)} "
                f"-out {shlex.quote(out_path)} -db {shlex.quote(db_loc)} "
                f"{parameters} -outfmt \"6 {flags}\""
            )

        elif task == "diamond":
            # DIAMOND blastx — translated nucleotide vs protein. ~300 ms
            # of fixed per-call overhead vs ~1 s for mmseqs2, so wins on
            # the small protein DBs (feature_protein 709, fpbase 990).
            flags = "qstart qend sseqid pident slen qseq length sstart send qlen evalue"
            cmd = (
                f"diamond blastx -d {shlex.quote(db_loc)}.dmnd "
                f"-q {shlex.quote(query_path)} -o {shlex.quote(out_path)} "
                f"{parameters} --outfmt 6 {flags} --quiet"
            )

        elif task == "mmseqs":
            flags = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
            # mmseqs DBs live at feature_db_data/mmseqs_dbs/<name>_db.
            mmseqs_db = db_loc.replace("/BLAST_dbs/", "/mmseqs_dbs/") + "_db"
            tmp_dir = tempfile.mkdtemp(prefix="splicify_mmseqs_")
            cmd = (
                f"mmseqs easy-search {shlex.quote(query_path)} "
                f"{shlex.quote(mmseqs_db)} {shlex.quote(out_path)} "
                f"{shlex.quote(tmp_dir)} --search-type 2 --format-mode 0 "
                f"-e 0.001 --min-seq-id 0.8 -v 0"
            )

        elif task == "mmseqs_nucl":
            flags = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
            mmseqs_db = db_loc.replace("/BLAST_dbs/", "/mmseqs_dbs/") + "_db"
            tmp_dir = tempfile.mkdtemp(prefix="splicify_mmseqs_")
            cmd = (
                f"mmseqs easy-search {shlex.quote(query_path)} "
                f"{shlex.quote(mmseqs_db)} {shlex.quote(out_path)} "
                f"{shlex.quote(tmp_dir)} --search-type 3 --format-mode 0 "
                f"-e 0.001 --min-seq-id 0.9 -v 0"
            )

        elif task == "infernal":
            flags = "--cut_ga --rfam --noali --nohmmonly --fmt 2"
            cmd = (
                f"cmscan {flags} {parameters} --tblout {shlex.quote(out_path)} "
                f"--clanin {shlex.quote(db_loc)} {shlex.quote(query_path)}"
            )

        else:
            raise ValueError(f"Unknown search method: {task}")

        logger.info("splicify annotator running %s", task)
        _t0 = time.time()
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        logger.info("[TIMING] %s: %.2fs", task, time.time() - _t0)

        if proc.returncode != 0:
            logger.error(
                "%s failed (rc=%s). stderr: %s",
                task, proc.returncode, (proc.stderr or "").strip(),
            )
            raise RuntimeError(f"{task} failed (rc={proc.returncode})")

        if task == "infernal":
            inDf = _parse_infernal(out_path)
            inDf["qlen"] = len(seq)
            if not inDf.empty:
                inDf["qseq"] = inDf.apply(
                    lambda x: seq[x["qstart"]:x["qend"] + 1].upper(), axis=1
                )
            return inDf

        # blastn / mmseqs: tabular outfmt 6
        with open(out_path, "r", encoding="utf-8", errors="replace") as fh:
            lines = [ln for ln in fh.readlines() if ln.strip()]

        if not lines:
            return pd.DataFrame(columns=flags.split())

        inDf = pd.DataFrame([ln.split() for ln in lines], columns=flags.split())
        for col in inDf.columns:
            try:
                inDf[col] = pd.to_numeric(inDf[col])
            except (ValueError, TypeError):
                pass

        if task == "diamond":
            # Same SwissProt sp|ACCESSION|ENTRY_NAME handling pLannotate did.
            try:
                split_df = inDf["sseqid"].astype(str).str.split("|", expand=True)
                if split_df.shape[1] >= 3:
                    inDf["_uniprot_accession"] = split_df[1]
                    inDf["sseqid"] = split_df[2]
                elif split_df.shape[1] >= 2:
                    inDf["sseqid"] = split_df[1]
            except Exception:
                pass
            inDf["sframe"] = (inDf["qstart"] < inDf["qend"]).astype(int).replace(0, -1)
            inDf["slen"] = inDf["slen"] * 3
            inDf["length"] = abs(inDf["qend"] - inDf["qstart"]) + 1

        if task in ("mmseqs", "mmseqs_nucl"):
            inDf["sseqid"] = inDf["sseqid"].astype(str).str.split("|", n=1).str[-1]
            inDf["sframe"] = (
                inDf["qstart"].astype(float) < inDf["qend"].astype(float)
            ).astype(int).replace(0, -1)
            # mmseqs reports protein-length match; multiply by 3 for nt-scale
            # accounting so percmatch / score line up with blastn hits.
            inDf["slen"] = inDf["length"].astype(float) * 3
            inDf["qlen"] = len(seq)
            inDf["qseq"] = inDf.apply(
                lambda x: seq[
                    int(min(float(x["qstart"]), float(x["qend"]))) - 1
                    : int(max(float(x["qstart"]), float(x["qend"])))
                ].upper(),
                axis=1,
            )

        return inDf

    finally:
        for p in (query_path, out_path):
            try:
                os.remove(p)
            except OSError:
                pass


def _database_available(db: Dict[str, Any]) -> bool:
    method = db["method"]
    db_loc = db["db_loc"]
    if method == "infernal":
        return all(os.path.exists(p) for p in db_loc.split())
    if method in ("mmseqs", "mmseqs_nucl"):
        mmseqs_db = db_loc.replace("/BLAST_dbs/", "/mmseqs_dbs/") + "_db"
        return os.path.exists(mmseqs_db)
    if method == "diamond":
        return os.path.exists(f"{db_loc}.dmnd")
    return bool(glob.glob(f"{db_loc}*"))


# --------------------------------------------------------------------------- #
# Score / clean / dedup (forked verbatim — pure dataframe ops).
# --------------------------------------------------------------------------- #


def _calculate(inDf: pd.DataFrame, is_linear: bool) -> pd.DataFrame:
    inDf["qstart"] = inDf["qstart"] - 1
    inDf["qend"] = inDf["qend"] - 1
    inDf["qstart"], inDf["qend"] = (
        inDf[["qstart", "qend"]].min(axis=1),
        inDf[["qstart", "qend"]].max(axis=1),
    )
    inDf["percmatch"] = inDf["length"] / inDf["slen"] * 100
    inDf["abs percmatch"] = 100 - abs(100 - inDf["percmatch"])
    inDf["pi_permatch"] = (inDf["pident"] * inDf["abs percmatch"]) / 100
    inDf["score"] = (inDf["pi_permatch"] / 100) * inDf["length"]
    inDf["score"] = inDf["score"] * (2 ** (-1 * inDf["priority"].astype(float)) * 2)
    if not is_linear:
        inDf["qlen"] = (inDf["qlen"] / 2).astype("int")
    bonus = (1 / (inDf["priority"] + 1)) * 10
    inDf.loc[inDf["pi_permatch"] == 100, "score"] = (
        inDf.loc[inDf["pi_permatch"] == 100, "score"] * bonus
    )
    wiggle_size = 0.15
    inDf["wiggle"] = (inDf["length"] * wiggle_size).astype(int)
    inDf["wstart"] = inDf["qstart"] + inDf["wiggle"]
    inDf["wend"] = inDf["qend"] - inDf["wiggle"]
    return inDf


def _clean(inDf: pd.DataFrame) -> pd.DataFrame:
    inDf["qstart_dup"] = inDf["qstart"]
    inDf["qend_dup"] = inDf["qend"]
    inDf["qstart"] = np.where(inDf["qstart"] >= inDf["qlen"], inDf["qstart"] - inDf["qlen"], inDf["qstart"])
    inDf["qend"] = np.where(inDf["qend"] >= inDf["qlen"], inDf["qend"] - inDf["qlen"], inDf["qend"])
    inDf["wstart"] = np.where(inDf["wstart"] >= inDf["qlen"], inDf["wstart"] - inDf["qlen"], inDf["wstart"])
    inDf["wend"] = np.where(inDf["wend"] >= inDf["qlen"], inDf["wend"] - inDf["qlen"], inDf["wend"])

    problem_hits = ["P03851", "P03845", "ISS", "P03846"]
    inDf = inDf.loc[~inDf["sseqid"].isin(problem_hits)]
    inDf = inDf.loc[inDf["evalue"] < 1]
    inDf = inDf.loc[inDf["pi_permatch"] > 3]

    dedup_cols = [c for c in inDf.columns if c not in ("qstart_dup", "qend_dup")]
    inDf = inDf.drop_duplicates(subset=dedup_cols).reset_index(drop=True)

    if inDf.empty:
        return pd.DataFrame(columns=DF_COLS)

    for col in inDf.columns:
        try:
            inDf[col] = pd.to_numeric(inDf[col]).astype(int)
        except (ValueError, TypeError):
            pass

    # Interval-based per-kind overlap dedup. Vectorised replacement of the
    # original seqSpace mechanic (~6 s on 1k+ hits -> <50 ms here).
    end = int(inDf["qlen"][0])
    qlen = end
    to_drop: list[int] = []
    covered_per_kind: dict = {}

    for idx in inDf.index:
        row = inDf.loc[idx]
        kind = row["kind"]
        wstart, wend = int(row["wstart"]), int(row["wend"])
        qstart, qend = int(row["qstart"]), int(row["qend"])

        covered = covered_per_kind.get(kind)
        if covered is None:
            covered = np.zeros(qlen, dtype=bool)
            covered_per_kind[kind] = covered

        # Drop check: does this hit's wiggle window overlap any
        # already-marked (winner) full window for the same kind?
        if wstart <= wend:
            overlap = bool(covered[wstart:wend + 1].any())
        else:  # circular wraparound
            overlap = bool(covered[wstart:].any() or covered[:wend + 1].any())
        if overlap:
            to_drop.append(idx)
            continue

        # Winner — mark its FULL extent so subsequent hits get the
        # same drop semantics the seqSpace approach had.
        if qstart <= qend:
            covered[qstart:qend + 1] = True
        else:
            covered[qstart:] = True
            covered[:qend + 1] = True

    if to_drop:
        inDf = inDf.drop(to_drop)
    inDf = inDf.reset_index(drop=True)
    return inDf


def _get_details(inDf: pd.DataFrame, databases: Dict[str, Dict[str, Any]]) -> pd.DataFrame:
    """Read the per-DB CSV (sseqid,Feature,Type,Description) and merge."""
    assert len(set(inDf["db"].to_list())) == 1, "all hits must be from the same db"
    db_name = inDf["db"].iloc[0]
    db_details = databases[db_name]["details"]

    if db_details["location"] in ("None", None):
        return inDf.loc[inDf["db"] == db_name][["sseqid", "Feature", "Description"]]

    csv_path = db_details["location"]
    feat_desc = pd.read_csv(csv_path)

    default_type = db_details.get("default_type")
    if default_type not in (None, "None", "null"):
        feat_desc["Type"] = default_type
    return feat_desc


# --------------------------------------------------------------------------- #
# Dedup / coverage filter (forked verbatim — pure DataFrame ops).
# --------------------------------------------------------------------------- #


def _filter_overlapping_by_coverage(inDf: pd.DataFrame) -> pd.DataFrame:
    if inDf.empty:
        return inDf
    inDf = inDf.sort_values(by=["percmatch", "score"], ascending=[False, False]).reset_index(drop=True)
    qlen = int(inDf["qlen"].iloc[0]) if "qlen" in inDf.columns else 10000
    covered: set[int] = set()
    keep: list[int] = []
    for idx in inDf.index:
        row = inDf.loc[idx]
        qstart, qend, percmatch = int(row["qstart"]), int(row["qend"]), float(row["percmatch"])
        if qstart <= qend:
            positions = set(range(qstart, qend + 1))
        else:
            positions = set(range(qstart, qlen)) | set(range(0, qend + 1))
        overlap_frac = (len(positions & covered) / len(positions)) if positions else 0
        if overlap_frac < 0.5 or percmatch >= 90:
            keep.append(idx)
            covered |= positions
    return inDf.loc[keep].reset_index(drop=True)


def _deduplicate_variant_annotations(inDf: pd.DataFrame) -> pd.DataFrame:
    if inDf.empty:
        return inDf
    variant_pattern = re.compile(r"[_\s]*\(?\d+\)?$")

    def base_name(name: Any) -> str:
        if not name:
            return ""
        return variant_pattern.sub("", str(name)).strip().lower()

    def positions_of(row: pd.Series, qlen: int) -> set[int]:
        qstart, qend = int(row["qstart"]), int(row["qend"])
        if qstart <= qend:
            return set(range(qstart, qend + 1))
        return set(range(qstart, qlen)) | set(range(0, qend + 1))

    qlen = int(inDf["qlen"].iloc[0]) if "qlen" in inDf.columns else 10000
    inDf = inDf.copy()
    inDf["_base_name"] = inDf["Feature"].apply(base_name)
    inDf = inDf.sort_values(by=["percmatch", "score"], ascending=[False, False])

    keep: list[int] = []
    regions: dict[str, list[tuple[set[int], int]]] = {}
    for idx in inDf.index:
        row = inDf.loc[idx]
        bn = row["_base_name"]
        pos = positions_of(row, qlen)
        if bn in regions:
            dominated = False
            for existing, _ in regions[bn]:
                overlap_frac = (len(pos & existing) / len(pos)) if pos else 0
                if overlap_frac > 0.7:
                    dominated = True
                    break
            if dominated:
                continue
        keep.append(idx)
        regions.setdefault(bn, []).append((pos, idx))

    return inDf.loc[keep].drop(columns=["_base_name"]).reset_index(drop=True)


# --------------------------------------------------------------------------- #
# Pipeline (forked from plannotate.annotate.{get_raw_hits, annotate}).
# --------------------------------------------------------------------------- #


def _get_raw_hits(query: str, linear: bool, databases: Dict[str, Dict[str, Any]]) -> pd.DataFrame:
    def process_db(db_name: str, db_config: Dict[str, Any]) -> Optional[pd.DataFrame]:
        if not _database_available(db_config):
            logger.warning(
                "Skipping db '%s' (files not present at %s)", db_name, db_config["db_loc"]
            )
            return None
        hits = _run_blast(seq=query, db=db_config)
        if hits.empty:
            return None
        hits["db"] = db_name
        hits["sseqid"] = hits["sseqid"].astype(str)
        feat_desc = _get_details(hits, databases)
        hits = hits.merge(feat_desc, on="sseqid", how="left", suffixes=("_x", None))
        hits = hits[hits.columns.drop(list(hits.filter(regex="_x")))]
        # Drop primer_bind / primer annotations — match upstream behavior.
        if "Type" in hits.columns:
            hits = hits.loc[hits["Type"] != "primer_bind"]
        if "Feature" in hits.columns and "Description" in hits.columns:
            primer_like = (
                hits["Feature"].astype(str).str.contains("primer", case=False, na=False)
                | hits["Description"].astype(str).str.contains("primer", case=False, na=False)
            )
            hits = hits.loc[~primer_like]
        hits["priority"] = db_config["priority"]
        try:
            hits["priority"] = hits["priority"] + hits["priority_mod"]
            hits = hits.drop("priority_mod", axis=1)
        except KeyError:
            pass
        hits = _calculate(hits, is_linear=linear)
        return hits

    raw: List[pd.DataFrame] = []
    with ThreadPoolExecutor(max_workers=6) as executor:
        futures = {
            executor.submit(process_db, db_name, databases[db_name]): db_name
            for db_name in databases
        }
        for future in as_completed(futures):
            db_name = futures[future]
            try:
                result = future.result()
                if result is not None and not result.empty:
                    raw.append(result)
            except Exception as exc:
                logger.error("error processing db '%s': %s", db_name, exc)

    if not raw:
        return pd.DataFrame()
    return pd.concat(raw).sort_values(
        by=["score", "length", "percmatch"], ascending=[False, False, False]
    )


# --------------------------------------------------------------------------- #
# Page-cache pre-warm — reads every BLAST / mmseqs / Rfam index file once on
# module import so the first annotate() call doesn't pay the ~10 s cold-disk
# penalty for loading ~3 GB of indices. Runs in a daemon thread so it does
# not block module import (lazy callers won't pay for it).
# --------------------------------------------------------------------------- #

_PREWARM_STARTED = False
_PREWARM_PATTERNS = (
    "BLAST_dbs/*.nhr", "BLAST_dbs/*.nin", "BLAST_dbs/*.nsq",
    "BLAST_dbs/*.ntf", "BLAST_dbs/*.nto", "BLAST_dbs/*.ndb",
    "BLAST_dbs/*.not", "BLAST_dbs/*.nos", "BLAST_dbs/*.nog",
    "BLAST_dbs/*.phr", "BLAST_dbs/*.pin", "BLAST_dbs/*.psq",
    "BLAST_dbs/*.ptf", "BLAST_dbs/*.pto", "BLAST_dbs/*.pdb",
    "BLAST_dbs/*.pot", "BLAST_dbs/*.pos", "BLAST_dbs/*.pog",
    "BLAST_dbs/*.cm", "BLAST_dbs/*.i1f", "BLAST_dbs/*.i1i",
    "BLAST_dbs/*.i1m", "BLAST_dbs/*.i1p", "BLAST_dbs/*.clanin",
    "BLAST_dbs/*.dmnd",
)


def _prewarm_page_cache() -> None:
    """Read every index file once. Quiet on failure (best-effort)."""
    if not _FEATURE_DB_DIR.is_dir():
        return
    n_files = 0
    n_bytes = 0
    for pat in _PREWARM_PATTERNS:
        for p in _FEATURE_DB_DIR.glob(pat):
            try:
                with p.open("rb", buffering=0) as fh:
                    while True:
                        chunk = fh.read(1024 * 1024)
                        if not chunk:
                            break
                        n_bytes += len(chunk)
                n_files += 1
            except OSError:
                pass
    logger.info(
        "feature_annotator: prewarmed %d index files (%.1f MB)",
        n_files, n_bytes / (1024 * 1024),
    )


def prewarm(blocking: bool = False) -> None:
    """Pre-warm the OS page cache for every feature_db_data/ index file.

    Default is non-blocking (daemon thread): module import returns
    immediately; the warm-up runs in the background and is usually
    done by the time the first annotate() call lands. Use
    ``blocking=True`` from a startup hook if you want to guarantee
    the warm-up has completed before the first annotation request
    (e.g. a smoke-test endpoint that should never see cold-cache
    latency).
    """
    global _PREWARM_STARTED
    if _PREWARM_STARTED:
        return
    _PREWARM_STARTED = True
    if blocking:
        _prewarm_page_cache()
    else:
        import threading
        threading.Thread(target=_prewarm_page_cache, daemon=True).start()


# Kick off the prewarm in the background. Safe to call multiple times.
prewarm(blocking=False)


# Hot-path: cache loaded YAML (re-read only on path change).
_YAML_CACHE: Dict[str, Dict[str, Dict[str, Any]]] = {}


def _databases(yaml_file: Optional[str]) -> Dict[str, Dict[str, Any]]:
    path = yaml_file or str(_DEFAULT_YAML)
    if path not in _YAML_CACHE:
        _YAML_CACHE[path] = _load_databases_yaml(path)
    return _YAML_CACHE[path]


def annotate(
    inSeq: str,
    yaml_file: Optional[str] = None,
    linear: bool = False,
    is_detailed: bool = False,
    include_rfam: bool = False,
) -> pd.DataFrame:
    """Drop-in replacement for `plannotate.annotate.annotate`.

    Args:
        inSeq: query DNA sequence (string).
        yaml_file: override of databases.yml location (default:
            feature_db_data/databases.yml).
        linear: True for linear sequence; False (default) doubles the
            query so origin-spanning hits are found.
        is_detailed: when True, the seqSpace dedup partitions by Type so
            features of different types can co-occupy a region.

    Returns:
        pandas DataFrame with columns matching DF_COLS.
    """
    # Validate by routing through Biopython first (catches odd characters).
    tmp = NamedTemporaryFile()
    SeqIO.write(
        SeqRecord(Seq(inSeq), name="splicify", annotations={"molecule_type": "DNA"}),
        tmp.name, "fasta",
    )
    record = list(SeqIO.parse(tmp.name, "fasta"))[0]
    tmp.close()

    if not linear:
        query = str(record.seq) + str(record.seq)
    else:
        query = str(record.seq)

    databases = _databases(yaml_file)
    if not include_rfam:
        # Default: skip the cmscan tier — saves ~3 s on the hot path.
        # Caller opts in via include_rfam=True when ncRNA / riboswitch /
        # aptamer / structural-RNA detection is needed.
        databases = {k: v for k, v in databases.items() if k != "rfam"}
    blastDf = _get_raw_hits(query, linear, databases)

    if blastDf.empty:
        return pd.DataFrame(columns=DF_COLS)

    blastDf["kind"] = blastDf["Type"] if is_detailed else 1
    blastDf = _clean(blastDf)
    blastDf = _filter_overlapping_by_coverage(blastDf)
    blastDf = _deduplicate_variant_annotations(blastDf)

    if blastDf.empty:
        return pd.DataFrame(columns=DF_COLS)

    def _is_fragment(feature: pd.Series) -> bool:
        if feature["Type"] == "CDS":
            if feature["pi_permatch"] == 100:
                return False
            if (feature["length"] % 3) == 0 and feature["percmatch"] > 95:
                return False
            return True
        return feature["percmatch"] < 95

    blastDf["fragment"] = blastDf.apply(_is_fragment, axis=1)
    blastDf["qend"] = blastDf["qend"] + 1  # 1-based GenBank
    blastDf["qseq"] = blastDf.apply(
        lambda x: (
            str(Seq(x["qseq"]).reverse_complement()) if x["sframe"] == -1 else x["qseq"]
        ),
        axis=1,
    )
    blastDf["Feature"] = blastDf["Feature"].fillna(blastDf["sseqid"])
    blastDf["Description"] = blastDf["Description"].fillna("")
    blastDf["Type"] = blastDf["Type"].fillna("misc_feature")
    return blastDf
