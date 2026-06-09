"""
Diff-first meta strategy routing for choosing cloning operators.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from enum import Enum
from typing import Any, Dict, List, Optional


class ModuleChange(str, Enum):
    UNCHANGED = "UNCHANGED"
    MODIFIED = "MODIFIED"
    NEW = "NEW"


class EditType(str, Enum):
    NO_CHANGE = "NO_CHANGE"
    SDM_SINGLE_EDIT = "SDM_SINGLE_EDIT"
    SINGLE_SWAP = "SINGLE_SWAP"
    MULTI_EDIT = "MULTI_EDIT"


@dataclass
class DiffResult:
    anchor_source: Optional[str]
    anchor_coverage: int
    edit_type: EditType
    routing_recommendation: str
    module_diffs: List[Dict[str, Any]]


class DiffRouter:
    """Analyze desired design vs anchor plasmid modules and recommend route."""

    def analyze(
        self,
        desired_modules: List[Dict[str, Any]],
        anchor_modules: Optional[List[Dict[str, Any]]] = None,
    ) -> Dict[str, Any]:
        if not desired_modules:
            return {
                "anchor_source": None,
                "anchor_coverage": 0,
                "edit_type": EditType.NO_CHANGE,
                "routing_recommendation": "no_cloning_needed",
                "module_diffs": [],
            }

        inferred_source, inferred_anchor = self._infer_anchor(desired_modules)
        anchor = anchor_modules if anchor_modules is not None else inferred_anchor
        anchor_lookup = {self._module_key(m): m for m in anchor}

        module_diffs: List[Dict[str, Any]] = []
        modified_records: List[Dict[str, Any]] = []

        for idx, mod in enumerate(desired_modules):
            key = self._module_key(mod)
            desired_seq = (mod.get("sequence") or "").upper()
            anchor_mod = anchor_lookup.get(key)
            name = mod.get("canonical_id") or mod.get("description") or f"module_{idx}"

            if anchor_mod is None:
                module_diffs.append({
                    "module_index": idx,
                    "module_name": name,
                    "status": ModuleChange.NEW,
                    "old_length": 0,
                    "new_length": len(desired_seq),
                })
                continue

            anchor_seq = (anchor_mod.get("sequence") or "").upper()
            if desired_seq == anchor_seq:
                module_diffs.append({
                    "module_index": idx,
                    "module_name": name,
                    "status": ModuleChange.UNCHANGED,
                    "old_length": len(anchor_seq),
                    "new_length": len(desired_seq),
                })
                continue

            edit = self._simple_diff(anchor_seq, desired_seq)
            row = {
                "module_index": idx,
                "module_name": name,
                "status": ModuleChange.MODIFIED,
                "old_length": len(anchor_seq),
                "new_length": len(desired_seq),
                "edit": edit,
            }
            module_diffs.append(row)
            modified_records.append(row)

        edit_type = self._classify_edit_type(module_diffs, modified_records)
        routing = self._recommend(edit_type)

        result = DiffResult(
            anchor_source=inferred_source,
            anchor_coverage=len(anchor),
            edit_type=edit_type,
            routing_recommendation=routing,
            module_diffs=module_diffs,
        )

        return {
            "anchor_source": result.anchor_source,
            "anchor_coverage": result.anchor_coverage,
            "edit_type": result.edit_type,
            "routing_recommendation": result.routing_recommendation,
            "module_diffs": result.module_diffs,
        }

    def _infer_anchor(self, desired_modules: List[Dict[str, Any]]) -> tuple[Optional[str], List[Dict[str, Any]]]:
        sources = [m.get("source") for m in desired_modules if m.get("source")]
        if not sources:
            return None, []
        source, _ = Counter(sources).most_common(1)[0]
        return source, [m for m in desired_modules if m.get("source") == source]

    def _module_key(self, module: Dict[str, Any]) -> str:
        return (
            module.get("canonical_id")
            or module.get("description")
            or module.get("role")
            or ""
        ).strip().lower()

    def _simple_diff(self, old_seq: str, new_seq: str) -> Dict[str, Any]:
        if not old_seq and new_seq:
            return {"kind": "insertion", "position": 0, "old_seq": "", "new_seq": new_seq}
        if old_seq and not new_seq:
            return {"kind": "deletion", "position": 0, "old_seq": old_seq, "new_seq": ""}

        # longest common prefix/suffix heuristic for a single local edit
        prefix = 0
        max_prefix = min(len(old_seq), len(new_seq))
        while prefix < max_prefix and old_seq[prefix] == new_seq[prefix]:
            prefix += 1

        suffix = 0
        max_suffix = min(len(old_seq) - prefix, len(new_seq) - prefix)
        while suffix < max_suffix and old_seq[-(suffix + 1)] == new_seq[-(suffix + 1)]:
            suffix += 1

        old_mid_end = len(old_seq) - suffix if suffix else len(old_seq)
        new_mid_end = len(new_seq) - suffix if suffix else len(new_seq)

        old_mid = old_seq[prefix:old_mid_end]
        new_mid = new_seq[prefix:new_mid_end]

        if len(old_mid) == len(new_mid):
            kind = "substitution"
        elif len(old_mid) < len(new_mid):
            kind = "insertion"
        else:
            kind = "deletion"

        return {
            "kind": kind,
            "position": prefix,
            "old_seq": old_mid,
            "new_seq": new_mid,
            "delta_bp": len(new_mid) - len(old_mid),
        }

    def _classify_edit_type(
        self,
        module_diffs: List[Dict[str, Any]],
        modified_records: List[Dict[str, Any]],
    ) -> EditType:
        modified_count = sum(1 for row in module_diffs if row["status"] == ModuleChange.MODIFIED)
        new_count = sum(1 for row in module_diffs if row["status"] == ModuleChange.NEW)

        if modified_count == 0 and new_count == 0:
            return EditType.NO_CHANGE

        if modified_count == 1 and new_count == 0:
            edit = modified_records[0].get("edit", {})
            old_len = len(edit.get("old_seq", ""))
            new_len = len(edit.get("new_seq", ""))
            if max(old_len, new_len) <= 60:
                return EditType.SDM_SINGLE_EDIT

        if modified_count + new_count == 1:
            return EditType.SINGLE_SWAP

        return EditType.MULTI_EDIT

    def _recommend(self, edit_type: EditType) -> str:
        if edit_type == EditType.NO_CHANGE:
            return "no_cloning_needed"
        if edit_type == EditType.SDM_SINGLE_EDIT:
            return "sdm_q5"
        if edit_type == EditType.SINGLE_SWAP:
            return "restriction_or_gibson"
        return "gibson_hifi"
