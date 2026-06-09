"""Inventory of PlasmidIndex objects + name/feature-driven resolution.

The interpreter agent works on either a single named plasmid or fans
out across every plasmid in the registry. resolve_plasmid() is the
canonical way to map a user-supplied free-text reference (filename,
gene name, "the one with EGFP and PuroR") to one or more concrete
indexes — or to surface that nothing matches.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from difflib import SequenceMatcher
from typing import Any, Callable, Optional

from agent_v2.interpreter.plasmid_index import PlasmidIndex


# Stop-words stripped from free-text plasmid references before scoring
# against feature names. Keeps the matcher focused on nouns.
_STOP_WORDS = {
    "the", "a", "an", "of", "with", "that", "which", "has", "have", "and",
    "or", "for", "on", "in", "to", "from", "by", "is", "are", "was", "were",
    "plasmid", "plasmids", "vector", "construct", "one", "where", "what",
    "how", "any", "all", "this", "these", "those", "uploaded", "showing",
    "carry", "carrying", "uses", "use", "contains", "containing",
}


def _normalise_name(s: str) -> str:
    s = (s or "").lower().strip()
    s = re.sub(r"[\s_\-\.]+", " ", s)
    return s


def _tokens(s: str) -> list[str]:
    return [t for t in re.findall(r"[a-z0-9α]+", (s or "").lower())
            if t not in _STOP_WORDS and len(t) >= 2]


@dataclass
class PlasmidRegistry:
    """Holds N PlasmidIndex objects. Methods either operate on one
    (by plasmid_id) or fan out across the inventory and aggregate.
    """
    indexes: dict[str, PlasmidIndex] = field(default_factory=dict)

    # ──────────────────────────────────────────────────────────────────
    def register(self, plasmid_id: str, envelope: dict[str, Any], name: Optional[str] = None) -> PlasmidIndex:
        idx = PlasmidIndex.from_envelope(plasmid_id, envelope, name=name)
        self.indexes[plasmid_id] = idx
        return idx

    def get(self, plasmid_id: str) -> Optional[PlasmidIndex]:
        return self.indexes.get(plasmid_id)

    def all(self) -> list[PlasmidIndex]:
        return list(self.indexes.values())

    def n(self) -> int:
        return len(self.indexes)

    # ──────────────────────────────────────────────────────────────────
    # Plasmid resolution
    # ──────────────────────────────────────────────────────────────────
    def resolve_plasmid(self, query: str) -> dict[str, Any]:
        """Map a free-text reference to one or more concrete plasmid
        indexes. Strategy:

          1. Exact / substring match against plasmid name OR plasmid_id.
          2. Distinguishing-annotation match: tokenise the query, score
             each plasmid by how many query tokens appear in its feature
             names / module names / gene names.
          3. Fuzzy SequenceMatcher fallback against names.

        Returns:
          {
            "ok": bool,
            "matches": [{plasmid_id, name, score, reason}],
            "method": "exact" | "substring" | "feature_overlap" | "fuzzy" | "none",
            "rejected": [{plasmid_id, name, reason}],  # surfaced when ok=False
          }
        """
        if not self.indexes:
            return {"ok": False, "matches": [], "method": "none",
                    "reason": "Registry is empty — no plasmids uploaded."}

        q_raw = (query or "").strip()
        q = _normalise_name(q_raw)
        if not q:
            return {"ok": False, "matches": [], "method": "none",
                    "reason": "Empty query."}

        # 1. Exact / substring on name or plasmid_id.
        exact, substring = [], []
        for idx in self.indexes.values():
            nm = _normalise_name(idx.name or "")
            pid = _normalise_name(idx.plasmid_id)
            if q == nm or q == pid:
                exact.append(idx)
            elif (q in nm and nm) or (q in pid and pid):
                substring.append(idx)
        if exact:
            return {
                "ok": True, "method": "exact",
                "matches": [{"plasmid_id": i.plasmid_id, "name": i.name,
                              "score": 1.0, "reason": "exact name match"}
                             for i in exact],
            }
        if substring:
            return {
                "ok": True, "method": "substring",
                "matches": [{"plasmid_id": i.plasmid_id, "name": i.name,
                              "score": 0.9, "reason": "name substring match"}
                             for i in substring],
            }

        # 2. Distinguishing-annotation scoring.
        toks = _tokens(q_raw)
        if toks:
            scored = []
            for idx in self.indexes.values():
                pool = self._feature_token_pool(idx)
                hits = sorted({t for t in toks if any(t in p for p in pool)})
                if hits:
                    score = len(hits) / len(toks)
                    scored.append({
                        "plasmid_id": idx.plasmid_id,
                        "name": idx.name,
                        "score": round(score, 3),
                        "reason": f"matches {len(hits)}/{len(toks)} distinguishing tokens: {hits}",
                    })
            scored.sort(key=lambda x: x["score"], reverse=True)
            # Only return scored matches with at least one hit AND the top
            # candidate strictly above any tied runner-up by 0.05+.
            if scored:
                # Filter to the top tier — within 0.05 of the highest score.
                top = scored[0]["score"]
                tier = [m for m in scored if m["score"] >= top - 0.05]
                if top >= 0.5 or len(toks) <= 2:
                    return {"ok": True, "method": "feature_overlap", "matches": tier,
                            "rejected": [m for m in scored if m not in tier][:5]}

        # 3. Fuzzy name matcher.
        fuzz = []
        for idx in self.indexes.values():
            nm = _normalise_name(idx.name or idx.plasmid_id)
            score = SequenceMatcher(None, q, nm).ratio()
            if score >= 0.6:
                fuzz.append({"plasmid_id": idx.plasmid_id, "name": idx.name,
                             "score": round(score, 3),
                             "reason": f"fuzzy name similarity {score:.2f}"})
        if fuzz:
            fuzz.sort(key=lambda x: x["score"], reverse=True)
            return {"ok": True, "method": "fuzzy", "matches": fuzz[:3]}

        return {
            "ok": False, "method": "none", "matches": [],
            "reason": (
                "No plasmid in the registry matched the query by name, "
                "distinguishing features, or fuzzy similarity."
            ),
            "registry_names": [i.name or i.plasmid_id for i in self.indexes.values()],
        }

    def _feature_token_pool(self, idx: PlasmidIndex) -> set[str]:
        """All searchable tokens from a plasmid's features, module names,
        and KB metadata. Used by distinguishing-annotation matching."""
        pool: set[str] = set()
        for a in idx.annotations():
            for k in ("name", "sseqid"):
                pool.update(_tokens(str(a.get(k, ""))))
            kb = a.get("kb_data") or {}
            for k in ("gene_name", "protein_name", "feature_name", "entry_name"):
                pool.update(_tokens(str(kb.get(k, ""))))
        for m in idx.modules():
            pool.update(_tokens(str(m.get("name", ""))))
            pool.update(_tokens(str(m.get("module_type", ""))))
            pool.update(_tokens(str(m.get("rule_id", ""))))
        return pool

    # ──────────────────────────────────────────────────────────────────
    # Fan-out
    # ──────────────────────────────────────────────────────────────────
    def fan_out(self, method_name: str, *, plasmid_id: Optional[str] = None, **kwargs) -> dict[str, Any]:
        """Call a PlasmidIndex method on a single plasmid (when
        plasmid_id is provided) or every plasmid (when None). Returns
        {ok, n_plasmids_searched, results: [...], no_results_in: [pid,...]}.
        Each result is the raw return value of the underlying method.
        """
        if plasmid_id is not None:
            idx = self.get(plasmid_id)
            if idx is None:
                return {"ok": False, "n_plasmids_searched": 0,
                        "reason": f"No plasmid registered with id={plasmid_id}",
                        "results": []}
            fn: Callable = getattr(idx, method_name)
            res = fn(**kwargs)
            return {
                "ok": True, "n_plasmids_searched": 1,
                "scope": "single",
                "results": _wrap_results(res, idx.plasmid_id),
                "no_results_in": [] if _has_results(res) else [idx.plasmid_id],
            }

        results: list[dict[str, Any]] = []
        empty: list[str] = []
        for idx in self.indexes.values():
            fn = getattr(idx, method_name)
            res = fn(**kwargs)
            wrapped = _wrap_results(res, idx.plasmid_id)
            if wrapped:
                results.extend(wrapped)
            if not _has_results(res):
                empty.append(idx.plasmid_id)
        return {
            "ok": True, "n_plasmids_searched": len(self.indexes),
            "scope": "inventory",
            "results": results,
            "no_results_in": empty,
        }


def _has_results(res: Any) -> bool:
    if res is None:
        return False
    if isinstance(res, list):
        return len(res) > 0
    if isinstance(res, dict):
        # AA lookup returns {"ok": True, ...} on success or {"ok": False, ...}
        # on out-of-range. Treat both as "had a result"; only None is empty.
        return True
    return True


def _wrap_results(res: Any, plasmid_id: str) -> list[dict[str, Any]]:
    """Coerce a method's return value into a flat list, tagging each
    item with the originating plasmid_id."""
    if res is None:
        return []
    if isinstance(res, list):
        return [{**(r if isinstance(r, dict) else {"value": r}), "plasmid_id": plasmid_id}
                for r in res]
    if isinstance(res, dict):
        return [{**res, "plasmid_id": plasmid_id}]
    return [{"value": res, "plasmid_id": plasmid_id}]
