"""
Visualization object builders for the AI Plasmid Design chat endpoint.
Converts raw endpoint responses into frontend-ready viz objects.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional


def build_pcr_viz(result: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    """
    Build PCR visualization from a /design-primers response.

    Annotation coordinates use the same convention as the frontend:
      - left primer:  [start, start+len)  direction=1
      - right primer: [start_3prime-len+1, start_3prime+1)  direction=-1
      - excluded region: direction=0
    """
    if not result:
        return None

    template = result.get("fragments_in", "")
    left_pos = result.get("left_pos")   # {"start": int, "len": int}
    right_pos = result.get("right_pos") # {"start_3prime": int, "len": int}
    excluded = result.get("excluded_region")  # {"start": int, "length": int}

    annotations: List[Dict[str, Any]] = []

    if left_pos:
        l_len = left_pos["len"]
        l_start = left_pos["start"]
        annotations.append({
            "name": f"Left Primer ({l_len} bp)",
            "start": l_start,
            "end": l_start + l_len,
            "direction": 1,
            "sequence": result.get("left_primer"),
            "tm": result.get("left_tm"),
        })

    if right_pos:
        r3 = right_pos["start_3prime"]
        r_len = right_pos["len"]
        annotations.append({
            "name": f"Right Primer ({r_len} bp)",
            "start": r3 - r_len + 1,
            "end": r3 + 1,
            "direction": -1,
            "sequence": result.get("right_primer"),
            "tm": result.get("right_tm"),
        })

    if excluded and excluded.get("length"):
        excl_start = excluded["start"]
        excl_len = excluded["length"]
        annotations.append({
            "name": f"Excluded Region ({excl_len} bp)",
            "start": excl_start,
            "end": excl_start + excl_len,
            "direction": 0,
        })

    return {
        "type": "pcr",
        "sequence": template,
        "annotations": annotations,
        "product_size": result.get("product_size"),
    }


def build_batch_pcr_viz_list(result: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Build a list of PCR viz objects from a /batch-design-primers response."""
    if not result or "results" not in result:
        return []
    viz_list = []
    for r in result["results"]:
        v = build_pcr_viz(r)
        if v is not None:
            v["meta"] = {
                "template_index": r.get("template_index", 0),
                "template_name": r.get("template_name", ""),
            }
            viz_list.append(v)
    return viz_list


def build_gibson_viz(result: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    """
    Extract the viz object from a /design-gibson-primers response and augment
    it with primer binding site annotations derived from primers_by_fragment.

    The endpoint builds fragment + overlap annotations but omits primer
    annotations (they were removed to avoid duplication with the n8n viz node).
    We rebuild them here from the anneal sequences and fragment positions.
    """
    if not result:
        return None

    viz = result.get("viz")
    if not viz:
        return None

    # Build a name→position map from the existing fragment annotations
    existing_anns = viz.get("annotations", [])
    frag_pos: Dict[str, Dict[str, Any]] = {
        a["name"]: a
        for a in existing_anns
        if "name" in a and "start" in a and "end" in a
        # Only plain fragment entries (no overlap/primer type tags)
        and "type" not in a
    }

    primer_anns: List[Dict[str, Any]] = []
    seq_len = len(result.get("viz", {}).get("sequence", ""))

    for p in result.get("primers_by_fragment", []):
        if not p.get("needs_primers"):
            continue

        frag_name = p.get("fragment", "")
        frag = frag_pos.get(frag_name)
        if not frag:
            continue

        frag_start: int = frag["start"]
        frag_end: int = frag["end"]

        fwd_anneal = p.get("forward_anneal_seq", "") or ""
        rev_anneal = p.get("reverse_anneal_seq", "") or ""
        fwd_primer = p.get("forward_primer", "")
        rev_primer = p.get("reverse_primer", "")
        fwd_extension = p.get("forward_extension_seq", "") or ""
        rev_extension = p.get("reverse_extension_seq", "") or ""
        fwd_tm = p.get("forward_anneal_tm")
        rev_tm = p.get("reverse_anneal_tm")

        # Forward primer: extension (from prev frag) + anneal at START of fragment
        if fwd_anneal:
            ann_len = len(fwd_anneal)
            ext_len = len(fwd_extension)
            total_len = len(fwd_primer)
            name = f"{frag_name} FWD primer ({total_len}bp)"
            ext_start = frag_start - ext_len
            if ext_start >= 0 or seq_len == 0:
                primer_anns.append({
                    "name": name,
                    "start": max(0, ext_start),
                    "end": frag_start + ann_len,
                    "direction": 1,
                    "sequence": fwd_primer,
                    "tm": fwd_tm,
                    "type": "primer",
                })
            else:
                # Extension wraps around origin — emit two same-named annotations
                primer_anns.append({
                    "name": name,
                    "start": seq_len + ext_start,  # ext_start is negative
                    "end": seq_len,
                    "direction": 1,
                    "sequence": fwd_primer,
                    "tm": fwd_tm,
                    "type": "primer",
                })
                primer_anns.append({
                    "name": name,
                    "start": 0,
                    "end": frag_start + ann_len,
                    "direction": 1,
                    "sequence": fwd_primer,
                    "tm": fwd_tm,
                    "type": "primer",
                })

        # Reverse primer: anneal at END of fragment + extension (into next frag)
        if rev_anneal:
            ann_len = len(rev_anneal)
            ext_len = len(rev_extension)
            total_len = len(rev_primer)
            name = f"{frag_name} REV primer ({total_len}bp)"
            ext_end = frag_end + ext_len
            if ext_end <= seq_len or seq_len == 0:
                primer_anns.append({
                    "name": name,
                    "start": frag_end - ann_len,
                    "end": ext_end,
                    "direction": -1,
                    "sequence": rev_primer,
                    "tm": rev_tm,
                    "type": "primer",
                })
            else:
                # Extension wraps beyond end — emit two same-named annotations
                primer_anns.append({
                    "name": name,
                    "start": frag_end - ann_len,
                    "end": seq_len,
                    "direction": -1,
                    "sequence": rev_primer,
                    "tm": rev_tm,
                    "type": "primer",
                })
                primer_anns.append({
                    "name": name,
                    "start": 0,
                    "end": ext_end - seq_len,
                    "direction": -1,
                    "sequence": rev_primer,
                    "tm": rev_tm,
                    "type": "primer",
                })

    return {
        **viz,
        "annotations": existing_anns + primer_anns,
    }


def build_inv_gib_viz(result: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    """
    Build a viz object from an /inv-gib response.
    Returns type "gibson" so the frontend renders it with the Gibson viewer.
    """
    if not result:
        return None

    fragments = result.get("fragments_in", [])
    seq = result.get("target_sequence", "")
    if not fragments and not seq:
        return None

    seq_len = len(seq) if seq else 0
    annotations: List[Dict[str, Any]] = []
    for frag in fragments:
        if not isinstance(frag, dict):
            continue
        start = frag.get("target_start", frag.get("start", 0))
        raw_end = frag.get("target_end", frag.get("end", start + frag.get("length_bp", 0)))
        frag_seq = frag.get("sequence", "")
        end = start + len(frag_seq) if frag_seq else raw_end
        frag_name = frag.get("name", "Fragment")
        direction = 1 if frag.get("source_orientation", "+") == "+" else -1
        is_synth = "Synthesis" in frag_name
        ann_base = {
            "direction": direction,
            "source": frag.get("source_inventory"),
            "type": "synthesis_gap" if is_synth else "inventory",
        }
        if seq_len > 0 and end > seq_len:
            # Fragment wraps around origin — emit two same-named annotations
            # so mergeOriginCrossingAnnotations in the frontend merges them
            annotations.append({**ann_base, "name": frag_name, "start": start, "end": seq_len})
            wrap_end = end - seq_len
            if wrap_end > 0:
                annotations.append({**ann_base, "name": frag_name, "start": 0, "end": wrap_end})
        else:
            annotations.append({**ann_base, "name": frag_name, "start": start, "end": end})

    return {
        "type": "gibson",
        "sequence": seq,
        "annotations": annotations,
    }

def build_sdm_viz(plan_dict: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    """
    Build SDM visualization from an SDMBuildPlan dict.
    
    Shows the mutated plasmid with annotations for:
    - Mutation site (deletion/insertion/substitution region)
    - Forward primer annealing region
    - Reverse primer annealing region
    - Mutation payload region
    
    Adapts the PCR visualization pattern for SDM.
    """
    if not plan_dict:
        return None
    
    primer_design = plan_dict.get("primer_design")
    if not primer_design:
        return None
    
    # Use mutated sequence for visualization
    sequence = plan_dict.get("mutated_sequence", "")
    if not sequence:
        return None
    
    edit_start = primer_design.get("edit_start", 0)
    old_seq = primer_design.get("old_sequence", "")
    new_seq = primer_design.get("new_sequence", "")
    mutation_type = plan_dict.get("mutation_type", "substitution")
    primer_strategy = plan_dict.get("primer_strategy", "back_to_back")
    
    fwd_primer = primer_design.get("forward_primer", "")
    rev_primer = primer_design.get("reverse_primer", "")
    fwd_tm = primer_design.get("forward_tm")
    rev_tm = primer_design.get("reverse_tm")
    
    fwd_anneal_seq = plan_dict.get("fwd_anneal_seq", "")
    rev_anneal_seq = plan_dict.get("rev_anneal_seq", "")
    
    annotations: List[Dict[str, Any]] = []
    
    # Mutation site annotation (in mutated sequence coordinates)
    if mutation_type == "deletion":
        # Deletion: show as a point marker where sequence was removed
        annotations.append({
            "name": f"Deletion site ({len(old_seq)} bp removed)",
            "start": edit_start,
            "end": edit_start + 1,  # Minimal span for visibility
            "direction": 0,
            "type": "sdm_mutation",
            "color": "#EF4444",  # Red
        })
    elif mutation_type == "insertion":
        # Insertion: show the inserted region
        annotations.append({
            "name": f"Insertion ({len(new_seq)} bp)",
            "start": edit_start,
            "end": edit_start + len(new_seq),
            "direction": 0,
            "type": "sdm_mutation",
            "color": "#10B981",  # Green
        })
    else:  # substitution
        annotations.append({
            "name": f"Substitution ({len(old_seq)}->{len(new_seq)} bp)",
            "start": edit_start,
            "end": edit_start + len(new_seq),
            "direction": 0,
            "type": "sdm_mutation",
            "color": "#F59E0B",  # Amber
        })
    
    # Primer annotations based on strategy
    if primer_strategy == "back_to_back":
        # Forward primer: starts at mutation, extends downstream
        fwd_anneal_start = edit_start + len(new_seq)
        fwd_anneal_end = fwd_anneal_start + len(fwd_anneal_seq)
        annotations.append({
            "name": f"Fwd Primer ({len(fwd_primer)} bp)",
            "start": edit_start,
            "end": fwd_anneal_end,
            "direction": 1,
            "sequence": fwd_primer,
            "tm": fwd_tm,
            "type": "primer",
        })
        
        # Reverse primer: anneals upstream of mutation
        rev_anneal_end = edit_start
        rev_anneal_start = max(0, rev_anneal_end - len(rev_anneal_seq))
        annotations.append({
            "name": f"Rev Primer ({len(rev_primer)} bp)",
            "start": rev_anneal_start,
            "end": rev_anneal_end,
            "direction": -1,
            "sequence": rev_primer,
            "tm": rev_tm,
            "type": "primer",
        })
        
    elif primer_strategy == "single_primer":
        # Mutagenic primer spans the mutation with flanking regions
        # fwd_anneal_seq contains both upstream and downstream flanks
        flank_len = len(fwd_anneal_seq) // 2
        primer_start = max(0, edit_start - flank_len)
        primer_end = edit_start + len(new_seq) + flank_len
        annotations.append({
            "name": f"Mutagenic Primer ({len(fwd_primer)} bp)",
            "start": primer_start,
            "end": min(primer_end, len(sequence)),
            "direction": 1,
            "sequence": fwd_primer,
            "tm": fwd_tm,
            "type": "primer",
        })
        
        # Reverse primer is downstream
        # Position it after the mutation region
        rev_start = edit_start + len(new_seq) + 150  # Approximate position
        if rev_start < len(sequence):
            annotations.append({
                "name": f"Rev Primer ({len(rev_primer)} bp)",
                "start": rev_start,
                "end": min(rev_start + len(rev_anneal_seq), len(sequence)),
                "direction": -1,
                "sequence": rev_primer,
                "tm": rev_tm,
                "type": "primer",
            })
            
    elif primer_strategy == "overlapping":
        # Both primers share the mutation region
        # Forward: upstream_anneal + first half of mutation
        fwd_start = max(0, edit_start - len(fwd_anneal_seq))
        fwd_end = edit_start + len(new_seq) // 2
        annotations.append({
            "name": f"Fwd Primer ({len(fwd_primer)} bp)",
            "start": fwd_start,
            "end": fwd_end,
            "direction": 1,
            "sequence": fwd_primer,
            "tm": fwd_tm,
            "type": "primer",
        })
        
        # Reverse: second half of mutation + downstream_anneal
        rev_start = edit_start + len(new_seq) // 2
        rev_end = edit_start + len(new_seq) + len(rev_anneal_seq)
        annotations.append({
            "name": f"Rev Primer ({len(rev_primer)} bp)",
            "start": rev_start,
            "end": min(rev_end, len(sequence)),
            "direction": -1,
            "sequence": rev_primer,
            "tm": rev_tm,
            "type": "primer",
        })
        
        # Show overlap region
        if plan_dict.get("overlap_seq"):
            annotations.append({
                "name": f"Primer Overlap ({len(new_seq)} bp)",
                "start": edit_start,
                "end": edit_start + len(new_seq),
                "direction": 0,
                "type": "overlap",
                "color": "#8B5CF6",  # Purple
            })
    
    return {
        "type": "design",
        "sequence": sequence,
        "annotations": annotations,
        "mutation_type": mutation_type,
        "primer_strategy": primer_strategy,
        "old_sequence": old_seq,
        "new_sequence": new_seq,
    }
