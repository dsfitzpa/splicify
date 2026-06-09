"""
Rule-based module detector for discrete boundary modules and high-confidence patterns.

This detector handles modules that LLM-based parsing struggles with:
- LTR/ITR-bounded viral payloads (lentiviral, AAV)
- T-DNA modules (Agrobacterium binary vectors)
- Flank:cargo patterns (Cre/loxP, FLP/FRT, Gateway, transposons)
- Standalone modules with discrete boundaries and high weights (≥0.9)
"""

import csv
import re
from typing import List, Dict, Any, Optional, Tuple
from pathlib import Path



# FRT (Flp Recombination Target) site sequences for motif scanning
FRT_MOTIFS = {
    'FRT': 'GAAGTTCCTATTCTCTAGAAAGTATAGGAACTTC',  # Minimal FRT (34bp)
    'FRT_full': 'GAAGTTCCTATTCCGAAGTTCCTATTCTCTAGAAAGTATAGGAACTTC',  # Full FRT with spacer (48bp)
    'FRT3': 'GAAGTTCCTATACTTTCTAGAGAATAGGAACTTC',  # FRT3 variant
    'FRT5': 'GAAGTTCCTATTCTTCAAATAGTATAGGAACTTC',  # FRT5 variant
}


class RuleBasedModuleDetector:
    """Detects modules using rule-based pattern matching."""

    def __init__(self, heuristics_csv_path: str):
        """
        Initialize detector with heuristics CSV.

        Args:
            heuristics_csv_path: Path to module_heuristics.csv
        """
        self.heuristics_csv_path = heuristics_csv_path
        self.rules = self._load_rules()

    def _load_rules(self) -> List[Dict[str, Any]]:
        """Load and parse heuristics CSV."""
        rules = []
        csv_path = Path(self.heuristics_csv_path)

        if not csv_path.exists():
            return rules

        with open(csv_path, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                # Parse weight
                try:
                    row['weight'] = float(row['weight'])
                except (ValueError, KeyError):
                    row['weight'] = 0.0

                # Parse features list
                if 'features' in row:
                    row['features'] = [f.strip() for f in row['features'].split(';')]

                rules.append(row)

        return rules

    def detect_modules(self, features: List[Dict[str, Any]], sequence: str = None) -> List[Dict[str, Any]]:
        """
        Detect all rule-based modules from features.

        Args:
            features: List of pLannotate features with keys:
                - name: Feature name
                - start: Start position
                - end: End position
                - strand: +1 or -1
                - type: Feature type
                - kb_class: (optional) Knowledge base feature class

        Returns:
            List of detected modules with structure:
                - module_type: Type of module
                - start: Start position
                - end: End position
                - name: Module name
                - features: List of feature indices used
                - weight: Confidence weight
                - detection_method: "rule_based"
        """
        modules = []

        # Detect LTR-bounded payloads (lentiviral)
        ltr_modules = self._detect_ltr_payload(features)
        modules.extend(ltr_modules)

        # Detect ITR-bounded payloads (AAV)
        itr_modules = self._detect_itr_payload(features)
        modules.extend(itr_modules)

        # Detect T-DNA modules (Agrobacterium binary vectors)
        tdna_modules = self._detect_tdna_module(features, sequence)
        modules.extend(tdna_modules)

        # Detect flank:cargo patterns
        flank_modules = self._detect_flank_cargo_modules(features)
        modules.extend(flank_modules)

        # Detect standalone high-weight modules
        standalone_modules = self._detect_standalone_modules(features)
        modules.extend(standalone_modules)

        # Detect replication origins by KB class
        ori_modules = self._detect_ori_by_kb_class(features)
        modules.extend(ori_modules)

        # Detect Gateway cloning cassettes
        gateway_modules = self._detect_gateway_cassettes(features)
        modules.extend(gateway_modules)

        # Detect floxed regions (loxP-flanked cassettes)
        floxed_modules = self._detect_floxed_regions(features)
        modules.extend(floxed_modules)

        # Detect FRT-flanked regions (Flp recombination)
        frt_modules = self._detect_frt_flanked_regions(features, sequence)
        modules.extend(frt_modules)

        # Detect transposon ITR-flanked regions
        transposon_modules = self._detect_transposon_flanked_regions(features)
        modules.extend(transposon_modules)

        # Detect baculovirus homology regions
        baculo_modules = self._detect_baculovirus_homology_regions(features)
        print(f'[DEBUG detect_modules] Baculovirus modules found: {len(baculo_modules)}')
        for m in baculo_modules:
            print(f'  - {m["name"]} ({m["start"]}-{m["end"]})')
        modules.extend(baculo_modules)

        # BAC-R-01: mini-F replicon (repE + sopABC)
        modules.extend(self._detect_bac_f_replicon(features))

        # MR-EBV-01: EBV oriP + EBNA1 episomal module
        modules.extend(self._detect_ebv_episomal(features))

        # MR-SV-01: SV40 promoter + SV40 ori replication+early module
        modules.extend(self._detect_sv40_replication_module(features))

        # REC-GW-03: Gateway destination cassette (attR1 + ccdB + CmR + attR2)
        modules.extend(self._detect_gateway_dest_cassette_composite(features))

        # MCS-04: Dual-phage IVT cloning cassette (T7 + MCS + T3/SP6)
        modules.extend(self._detect_ivt_cloning_cassette(features))

        # INS-01/02: Paired cHS4 / β-globin insulator bracket
        modules.extend(self._detect_insulated_expression_block(features))

        # REC-INT-01: Phage integrase landing pad (attP/attB)
        modules.extend(self._detect_integrase_landing_pad(features))

        # CSEL-01..05: Counter-selection CDS (ccdB, sacB, DTA, barnase, codA)
        modules.extend(self._detect_counter_selection_cds(features))

        # REG-TET-01: Tet regulator cassette (tTA / rtTA / Tet-On 3G)
        modules.extend(self._detect_tet_regulator_cassette(features))

        # REG-AID-01: AID degron system (OsTIR1 + AID)
        modules.extend(self._detect_aid_degron_system(features))

        # REG-DIM-01: FKBP/FRB chemical dimerization
        modules.extend(self._detect_fkbp_frb_dimerization(features))

        # LAC-BW-01..05: lacZα blue/white screening modules (variant-specific)
        modules.extend(self._detect_lac_blue_white_module(features, sequence))
        modules.extend(self._detect_lac_alpha_disrupted_module(features, sequence))

        # POL3-01: Pol III guide-RNA expression cassette (U6/H1/7SK + gRNA scaffold)
        modules.extend(self._detect_pol3_expression_cassette(features, sequence))

        # SEL-BAC-01: Bacterial selection cassette (AmpR promoter + AmpR/KanR/CmR/TcR)
        modules.extend(self._detect_bacterial_selection_cassette(features))

        # SEL-MAM-01: Mammalian selection cassette (weak promoter + PuroR/NeoR/HygR/BSD/ZeoR + polyA)
        modules.extend(self._detect_mammalian_selection_cassette(features))

        return modules

    @staticmethod
    def _name_matches(feat: Dict[str, Any], patterns) -> bool:
        """Case-insensitive substring match of any pattern against feature name.

        Underscores in the incoming name are normalised to spaces so patterns
        like `'t7 promoter'` match pLannotate's `'T7_promoter'` output. Pattern
        strings are assumed to already be in normalised form.
        """
        raw = (feat.get('name') or '').lower()
        name = raw.replace('_', ' ')
        if isinstance(patterns, str):
            patterns = [patterns]
        return any(p.lower() in name for p in patterns)

    @staticmethod
    def _find_first(features, patterns):
        for i, f in enumerate(features):
            if RuleBasedModuleDetector._name_matches(f, patterns):
                return i, f
        return None, None

    @staticmethod
    def _find_all(features, patterns):
        return [(i, f) for i, f in enumerate(features)
                if RuleBasedModuleDetector._name_matches(f, patterns)]


    # Type IIs enzymes commonly used in sgRNA Golden Gate (BsmBI / BbsI /
    # BsaI / SapI / Esp3I / AarI). Recognition site only (cut offset not
    # needed for site presence/coords).
    _TYPE_IIS_FOR_GG = {
        "BsmBI": "CGTCTC",
        "Esp3I": "CGTCTC",   # isoschizomer
        "BbsI":  "GAAGAC",
        "BsaI":  "GGTCTC",
        "SapI":  "GCTCTTC",
        "AarI":  "CACCTGC",
    }

    @staticmethod
    def _rc(seq):
        comp = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
        return seq.translate(comp)[::-1]

    def _detect_sgrna_golden_gate(self, features, sequence, prom_feat, scaff_feat):
        """Check whether a Stuffer + matching flanking Type IIs sites
        live between prom_feat and scaff_feat. Returns a dict with the
        Stuffer span + enzyme + flanking site coords, or None if no
        such cassette is detectable."""
        if not sequence:
            return None
        # Find Stuffer-like features
        stuffer_pats = ("stuffer", "filler", "dropout")
        stuffer = None
        for f in features:
            nm = (f.get("name") or "").lower()
            if any(p in nm for p in stuffer_pats):
                if (prom_feat["end"] <= f["start"] <= scaff_feat["start"]
                        or scaff_feat["end"] <= f["start"] <= prom_feat["start"]):
                    stuffer = f
                    break
        if stuffer is None:
            return None
        s_start, s_end = int(stuffer["start"]), int(stuffer["end"])

        # For each Type IIs enzyme, scan the sequence for sites within 20 bp
        # of EITHER end of the Stuffer. Require a hit on each side (forward
        # or reverse-complement) for the same enzyme.
        seq_upper = sequence.upper()
        rc_lookup = {name: self._rc(site) for name, site in self._TYPE_IIS_FOR_GG.items()}
        window = 20

        for enzyme, site in self._TYPE_IIS_FOR_GG.items():
            sites_fwd = []
            sites_rev = []
            start_idx = 0
            while True:
                hit = seq_upper.find(site, start_idx)
                if hit < 0:
                    break
                sites_fwd.append(hit)
                start_idx = hit + 1
            start_idx = 0
            rc_site = rc_lookup[enzyme]
            while True:
                hit = seq_upper.find(rc_site, start_idx)
                if hit < 0:
                    break
                sites_rev.append(hit)
                start_idx = hit + 1
            all_sites = sites_fwd + sites_rev
            if not all_sites:
                continue
            # 5' flank: site within 20 bp upstream of stuffer start
            near_5 = [h for h in all_sites if s_start - window <= h <= s_start + 5]
            # 3' flank: site within 20 bp downstream of stuffer end
            near_3 = [h for h in all_sites if s_end - 5 <= h <= s_end + window]
            if near_5 and near_3:
                return {
                    "enzyme": enzyme,
                    "recognition_site": site,
                    "stuffer_name": stuffer.get("name", "Stuffer"),
                    "stuffer_start": s_start,
                    "stuffer_end": s_end,
                    "site_5prime": min(near_5),
                    "site_3prime": min(near_3),
                    "notes": (
                        f"sgRNA Golden Gate cassette: {enzyme} sites flank "
                        f"the {stuffer.get('name', 'Stuffer')} ("
                        f"5' site at {min(near_5)}, 3' site at {min(near_3)}). "
                        f"Cut + ligate a 20 bp guide oligo to replace the Stuffer."
                    ),
                }
        return None

    def _detect_pol3_expression_cassette(self, features: List[Dict[str, Any]], sequence: str = None) -> List[Dict[str, Any]]:
        """POL3-01: Pol III guide-RNA expression cassette.

        Pattern: [U6 | H1 | 7SK promoter] → (≤500 bp gap) → [gRNA scaffold | tracrRNA | sgRNA | pegRNA].
        Emits: guide_expression_cassette with submodule list [pol3_promoter, sgrna_scaffold].
        Weight: 0.9
        """
        modules = []
        promoter_pats = ['u6 promoter', 'h1 promoter', '7sk promoter', 'hu6', 'mu6']
        scaffold_pats = ['grna scaffold', 'sgrna scaffold', 'tracrrna', 'pegrna scaffold',
                         'sgrna', 'grna_scaffold']

        promoters = self._find_all(features, promoter_pats)
        scaffolds = self._find_all(features, scaffold_pats)
        if not promoters or not scaffolds:
            return modules

        for p_idx, p_feat in promoters:
            p_end = p_feat.get('end', 0)
            p_strand = p_feat.get('strand', 1)
            # pick the closest downstream scaffold on same strand, within 5000 bp
            best = None
            for s_idx, s_feat in scaffolds:
                if s_feat.get('strand', 1) != p_strand:
                    continue
                s_start = s_feat.get('start', 0)
                gap = s_start - p_end
                if p_strand == 1 and 0 <= gap <= 5000:
                    if best is None or gap < best[0]:
                        best = (gap, s_idx, s_feat)
                elif p_strand == -1:
                    # reverse-strand: scaffold should sit upstream of promoter
                    rev_gap = p_feat.get('start', 0) - s_feat.get('end', 0)
                    if 0 <= rev_gap <= 5000 and (best is None or rev_gap < best[0]):
                        best = (rev_gap, s_idx, s_feat)
            if best is None:
                continue
            _, s_idx, s_feat = best
            start = min(p_feat['start'], s_feat['start'])
            end = max(p_feat['end'], s_feat['end'])

            submods = [
                {
                    'module_type': 'pol3_promoter',
                    'start': p_feat['start'], 'end': p_feat['end'],
                    'strand': p_strand,
                    'name': p_feat.get('name', 'pol3_promoter'),
                },
                {
                    'module_type': 'sgrna_scaffold',
                    'start': s_feat['start'], 'end': s_feat['end'],
                    'strand': s_feat.get('strand', 1),
                    'name': s_feat.get('name', 'sgRNA scaffold'),
                },
            ]

            # Check for sgRNA Golden Gate cloning: Stuffer between promoter and
            # scaffold + Type IIs sites within 20 bp of the Stuffer's ends.
            gg_meta = self._detect_sgrna_golden_gate(features, sequence,
                                                     p_feat, s_feat) if sequence else None
            if gg_meta:
                submods.append({
                    "module_type": "stuffer_module",
                    "start": gg_meta["stuffer_start"],
                    "end": gg_meta["stuffer_end"],
                    "strand": p_strand,
                    "name": gg_meta.get("stuffer_name", "Stuffer"),
                })

            modules.append({
                'module_type': 'guide_expression_cassette',
                'start': start,
                'end': end,
                'strand': p_strand,
                'name': f"Pol III guide cassette ({p_feat.get('name', 'U6')})",
                'features': [p_idx, s_idx],
                'submodules': submods,
                'weight': 0.9,
                'detection_method': 'rule_based',
                'rule_id': 'POL3-GG-01' if gg_meta else 'POL3-01',
                'golden_gate': gg_meta,
                'notes': f"Pol III promoter {p_feat.get('name')} + scaffold {s_feat.get('name')} within 5000 bp",
            })
        return modules

    def _detect_bacterial_selection_cassette(self, features: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """SEL-BAC-01: Bacterial selection cassette.

        Pattern: [AmpR promoter | cat promoter | bla promoter | KanR promoter | tet promoter]
                 → [AmpR | KanR | CmR | TcR | SmR | ZeoR (bacterial)]
                 [→ terminator]
        Weight: 0.9
        """
        modules = []
        promoter_pats = ['ampr promoter', 'bla promoter', 'cat promoter', 'kanr promoter',
                         'kan promoter', 'tet promoter', 'tetr promoter', 'sm promoter']
        cds_pats = ['ampr', 'bla', 'kanr', 'neor/kanr', 'cmr', 'cat', 'tcr', 'tetr', 'smr',
                    'aadA', 'gmr']

        promoters = self._find_all(features, promoter_pats)
        cdss = self._find_all(features, cds_pats)
        if not promoters or not cdss:
            return modules

        used_cds = set()
        for p_idx, p_feat in promoters:
            p_strand = p_feat.get('strand', 1)
            p_start, p_end = p_feat['start'], p_feat['end']
            best = None
            for c_idx, c_feat in cdss:
                if c_idx in used_cds:
                    continue
                if c_feat.get('strand', 1) != p_strand:
                    continue
                if p_strand == 1:
                    gap = c_feat['start'] - p_end
                else:
                    gap = p_start - c_feat['end']
                # bacterial selection: CDS is immediately downstream on same strand, ≤200 bp gap
                if 0 <= gap <= 200 and (best is None or gap < best[0]):
                    best = (gap, c_idx, c_feat)
            if best is None:
                continue
            _, c_idx, c_feat = best
            used_cds.add(c_idx)

            start = min(p_start, c_feat['start'])
            end = max(p_end, c_feat['end'])

            submods = [
                {'module_type': 'bacterial_promoter',
                 'start': p_start, 'end': p_end, 'strand': p_strand,
                 'name': p_feat.get('name', 'bacterial_promoter')},
                {'module_type': 'resistance_cds',
                 'start': c_feat['start'], 'end': c_feat['end'],
                 'strand': c_feat.get('strand', 1),
                 'name': c_feat.get('name', 'resistance_cds')},
            ]

            modules.append({
                'module_type': 'bacterial_selection_cassette',
                'start': start,
                'end': end,
                'strand': p_strand,
                'name': f"Bacterial selection ({c_feat.get('name')})",
                'features': [p_idx, c_idx],
                'submodules': submods,
                'weight': 0.9,
                'detection_method': 'rule_based',
                'rule_id': 'SEL-BAC-01',
                'notes': f"{p_feat.get('name')} → {c_feat.get('name')} (gap {best[0]} bp)",
            })
        return modules

    def _detect_mammalian_selection_cassette(self, features: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """SEL-MAM-01: Mammalian selection cassette.

        Pattern: [PGK | SV40 | TK | EF-1α core | EM7 | RSV | CMV] → [PuroR | NeoR | HygR | BSD | ZeoR | BleoR]
                 → [polyA]
        Weight: 0.9
        """
        modules = []
        promoter_pats = ['pgk promoter', 'sv40 promoter', 'hsv tk promoter', 'tk promoter',
                         'ef-1α core', 'ef1 core', 'ef-1a core', 'em7 promoter', 'em7',
                         'rsv promoter']
        cds_pats = ['puror', 'puro', 'neor', 'neo', 'hygr', 'hyg', 'bsd', 'blast',
                    'zeor', 'ble', 'bleor', 'dhfr']
        polya_pats = ['poly(a)', 'polya', 'poly(a) signal', 'sv40 poly', 'bgh poly',
                      'hgh poly', 'β-globin poly', 'rb_glob_pa']

        promoters = self._find_all(features, promoter_pats)
        cdss = self._find_all(features, cds_pats)
        polyas = self._find_all(features, polya_pats)
        if not promoters or not cdss:
            return modules

        used_cds = set()
        for p_idx, p_feat in promoters:
            p_strand = p_feat.get('strand', 1)
            p_start, p_end = p_feat['start'], p_feat['end']
            best = None
            for c_idx, c_feat in cdss:
                if c_idx in used_cds:
                    continue
                if c_feat.get('strand', 1) != p_strand:
                    continue
                if p_strand == 1:
                    gap = c_feat['start'] - p_end
                else:
                    gap = p_start - c_feat['end']
                if 0 <= gap <= 5000 and (best is None or gap < best[0]):
                    best = (gap, c_idx, c_feat)
            if best is None:
                continue
            _, c_idx, c_feat = best
            used_cds.add(c_idx)

            # find polyA downstream of CDS (same strand, ≤1000 bp)
            polya_idx = None
            polya_feat = None
            best_pa = None
            for pa_idx, pa_feat in polyas:
                if pa_feat.get('strand', 1) != p_strand:
                    continue
                if p_strand == 1:
                    pa_gap = pa_feat['start'] - c_feat['end']
                else:
                    pa_gap = c_feat['start'] - pa_feat['end']
                if 0 <= pa_gap <= 1000 and (best_pa is None or pa_gap < best_pa[0]):
                    best_pa = (pa_gap, pa_idx, pa_feat)
            if best_pa is not None:
                _, polya_idx, polya_feat = best_pa

            start = min(p_start, c_feat['start'])
            end = max(p_end, c_feat['end'])
            if polya_feat is not None:
                end = max(end, polya_feat['end'])
                start = min(start, polya_feat['start'])

            submods = [
                {'module_type': 'mammalian_promoter',
                 'start': p_start, 'end': p_end, 'strand': p_strand,
                 'name': p_feat.get('name', 'mammalian_promoter')},
                {'module_type': 'resistance_cds',
                 'start': c_feat['start'], 'end': c_feat['end'],
                 'strand': c_feat.get('strand', 1),
                 'name': c_feat.get('name', 'resistance_cds')},
            ]
            if polya_feat is not None:
                submods.append({
                    'module_type': 'polyA_signal',
                    'start': polya_feat['start'], 'end': polya_feat['end'],
                    'strand': polya_feat.get('strand', 1),
                    'name': polya_feat.get('name', 'polyA'),
                })

            feats_idx = [p_idx, c_idx]
            if polya_idx is not None:
                feats_idx.append(polya_idx)

            modules.append({
                'module_type': 'mammalian_selection_cassette',
                'start': start,
                'end': end,
                'strand': p_strand,
                'name': f"Mammalian selection ({c_feat.get('name')})",
                'features': feats_idx,
                'submodules': submods,
                'weight': 0.9 if polya_feat is not None else 0.85,
                'detection_method': 'rule_based',
                'rule_id': 'SEL-MAM-01',
                'notes': f"{p_feat.get('name')} → {c_feat.get('name')}" + (f" → {polya_feat.get('name')}" if polya_feat is not None else ''),
            })
        return modules

    def _detect_bac_f_replicon(self, features):
        """BAC-R-01: mini-F replicon = repE + sopA + sopB + sopC (optionally ori2) adjacent."""
        repE = self._find_all(features, ['repE'])
        sopA = self._find_all(features, ['sopA'])
        sopB = self._find_all(features, ['sopB'])
        sopC = self._find_all(features, ['sopC'])
        ori2 = self._find_all(features, ['ori2'])
        # Require repE and at least two of sopA/B/C
        have_sop = sum(1 for s in (sopA, sopB, sopC) if s) >= 2
        if not repE or not have_sop:
            return []
        feat_idxs = [i for group in (repE, sopA, sopB, sopC, ori2) for i, _ in group]
        feats_hit = [features[i] for i in feat_idxs]
        start = min(f['start'] for f in feats_hit)
        end = max(f['end'] for f in feats_hit)
        # Reject if span is unreasonably large (>15 kb → parts are on opposite sides of a circular plasmid)
        if end - start > 15000:
            return []
        return [{
            'module_type': 'bac_f_replicon',
            'start': start, 'end': end,
            'strand': 1,
            'name': 'mini-F replicon (repE + sopABC)',
            'features': sorted(set(feat_idxs)),
            'weight': 0.98,
            'detection_method': 'rule_based',
            'rule_id': 'BAC-R-01',
            'notes': 'repE initiator + sopABC partition system'
        }]

    def _detect_ebv_episomal(self, features):
        """MR-EBV-01: oriP + EBNA1 give EBV-maintained episome in human cells."""
        orip = self._find_all(features, ['oriP', 'mini-oriP'])
        ebna = self._find_all(features, ['EBNA1', 'EBNA-1'])
        if not orip or not ebna:
            return []
        feat_idxs = [i for grp in (orip, ebna) for i, _ in grp]
        feats_hit = [features[i] for i in feat_idxs]
        start = min(f['start'] for f in feats_hit)
        end = max(f['end'] for f in feats_hit)
        if end - start > 20000:
            # oriP and EBNA1 may be on opposite sides; still emit two separate sub-hits? Keep composite anyway.
            pass
        return [{
            'module_type': 'ebv_episomal_module',
            'start': start, 'end': end,
            'strand': 1,
            'name': 'EBV oriP + EBNA1 episomal maintenance',
            'features': sorted(set(feat_idxs)),
            'weight': 0.97,
            'detection_method': 'rule_based',
            'rule_id': 'MR-EBV-01',
            'notes': 'EBNA1 binds oriP for stable episomal replication in human cells'
        }]

    def _detect_sv40_replication_module(self, features):
        """MR-SV-01: SV40 promoter immediately adjacent to SV40 ori (within ~50 bp)."""
        prom = self._find_all(features, ['SV40 promoter'])
        ori = self._find_all(features, ['SV40 ori'])
        if not prom or not ori:
            return []
        modules = []
        for pi, pf in prom:
            for oi, of in ori:
                dist = min(abs(pf['end'] - of['start']), abs(of['end'] - pf['start']))
                # Accept overlap or adjacent within 50 bp
                overlap = not (pf['end'] <= of['start'] or of['end'] <= pf['start'])
                if overlap or dist <= 50:
                    start = min(pf['start'], of['start'])
                    end = max(pf['end'], of['end'])
                    modules.append({
                        'module_type': 'sv40_replication_module',
                        'start': start, 'end': end,
                        'strand': pf.get('strand', 1),
                        'name': 'SV40 promoter + SV40 ori',
                        'features': [pi, oi],
                        'weight': 0.97,
                        'detection_method': 'rule_based',
                        'rule_id': 'MR-SV-01',
                        'notes': 'SV40 early promoter overlaps ori; gives T-antigen-dependent replication'
                    })
                    break
        return modules

    def _detect_gateway_dest_cassette_composite(self, features):
        """REC-GW-03: canonical Gateway destination cassette attR1 → ccdB → CmR → attR2."""
        attR1 = self._find_first(features, ['attR1'])
        attR2 = self._find_first(features, ['attR2'])
        ccdb = self._find_first(features, ['ccdB'])
        cmr = self._find_first(features, ['CmR', 'CamR', 'chloramphenicol'])
        if not (attR1[1] and attR2[1] and ccdb[1] and cmr[1]):
            return []
        feats_hit = [attR1[1], attR2[1], ccdb[1], cmr[1]]
        start = min(f['start'] for f in feats_hit)
        end = max(f['end'] for f in feats_hit)
        # Require ccdB and CmR to lie between the att sites
        r1s, r2s = attR1[1]['start'], attR2[1]['start']
        lo, hi = min(r1s, r2s), max(r1s, r2s)
        inside = lo <= ccdb[1]['start'] <= hi and lo <= cmr[1]['start'] <= hi
        if not inside:
            return []
        return [{
            'module_type': 'gateway_dest_cassette',
            'start': start, 'end': end,
            'strand': 1,
            'name': 'Gateway DEST cassette (attR1-ccdB-CmR-attR2)',
            'features': sorted({attR1[0], attR2[0], ccdb[0], cmr[0]}),
            'weight': 0.98,
            'detection_method': 'rule_based',
            'rule_id': 'REC-GW-03',
            'notes': 'Diagnostic destination vector cassette awaiting LR reaction'
        }]

    def _detect_ivt_cloning_cassette(self, features):
        """MCS-04: T7 promoter + MCS + (T3 or SP6 promoter) flanking MCS for in vitro transcription."""
        t7 = self._find_all(features, ['T7 promoter'])
        t3 = self._find_all(features, ['T3 promoter'])
        sp6 = self._find_all(features, ['SP6 promoter'])
        mcs = self._find_all(features, ['MCS', 'polylinker', 'multiple cloning'])
        if not t7 or (not t3 and not sp6):
            return []
        other = t3 if t3 else sp6
        other_name = 'T3' if t3 else 'SP6'
        modules = []
        for ti, tf in t7:
            for oi, of in other:
                dist = abs(tf['start'] - of['start'])
                if dist > 500:
                    continue
                start = min(tf['start'], of['start'])
                end = max(tf['end'], of['end'])
                mcs_idxs = [mi for mi, mf in mcs if start <= mf['start'] and mf['end'] <= end]
                modules.append({
                    'module_type': 'ivt_cloning_cassette',
                    'start': start, 'end': end,
                    'strand': 1,
                    'name': f'IVT MCS (T7 + MCS + {other_name})',
                    'features': sorted({ti, oi, *mcs_idxs}),
                    'weight': 0.95,
                    'detection_method': 'rule_based',
                    'rule_id': 'MCS-04',
                    'notes': 'Dual phage promoters flanking MCS — in vitro transcription of either strand'
                })
                break
        return modules

    def _detect_insulated_expression_block(self, features):
        """INS-01/02: paired cHS4 or β-globin insulator bracketing a payload."""
        ins = self._find_all(features, ['cHS4', 'β-globin insulator', 'beta-globin insulator',
                                        "5' β-globin", "3' β-globin", "5' beta-globin", "3' beta-globin"])
        if len(ins) < 2:
            return []
        ins_sorted = sorted(ins, key=lambda t: t[1]['start'])
        i5, f5 = ins_sorted[0]
        i3, f3 = ins_sorted[-1]
        if f3['start'] - f5['end'] < 100:  # too close together — not bracketing a payload
            return []
        return [{
            'module_type': 'insulated_expression_block',
            'start': f5['start'], 'end': f3['end'],
            'strand': 1,
            'name': 'Insulator-bracketed block (cHS4 pair)',
            'features': [idx for idx, _ in ins_sorted],
            'weight': 0.9,
            'detection_method': 'rule_based',
            'rule_id': 'INS-01+INS-02',
            'notes': f"Bracketed by {f5['name']} and {f3['name']}"
        }]

    def _detect_integrase_landing_pad(self, features):
        """REC-INT-01: phage integrase attP/attB landing pad sites (standalone)."""
        patterns = ['λ attP', 'lambda attP', 'λ attB', 'lambda attB',
                    'φC31 attP', 'phiC31 attP', 'phi C31 attP', 'phiC31 attB', 'φC31 attB',
                    'Bxb1 attP', 'Bxb1 attB', 'HK022 attP', 'phage 186 attP',
                    'φ80 attP', 'phi80 attP', 'φ21 attP', 'phi21 attP']
        hits = self._find_all(features, patterns)
        modules = []
        for idx, f in hits:
            modules.append({
                'module_type': 'integrase_landing_pad',
                'start': f['start'], 'end': f['end'],
                'strand': f.get('strand', 1),
                'name': f['name'],
                'features': [idx],
                'weight': 0.9,
                'detection_method': 'rule_based',
                'rule_id': 'REC-INT-01',
                'notes': 'Phage integrase attachment site for single-copy integration'
            })
        return modules

    def _detect_counter_selection_cds(self, features):
        """CSEL-01..05: ccdB, sacB, DTA, barnase, codA as counter-selection CDS modules."""
        table = [
            ('ccdB', 0.97, 'CSEL-01', 'CcdB DNA gyrase poison — Gateway destination counter-selection'),
            ('sacB', 0.95, 'CSEL-02', 'B. subtilis levansucrase — sucrose counter-selection'),
            ('DTA', 0.95, 'CSEL-03', 'Diphtheria toxin A chain — negative selection in mammalian cells'),
            ('barnase', 0.9, 'CSEL-04', 'Bacillus RNase — toxic without barstar'),
            ('codA', 0.9, 'CSEL-05', 'Cytosine deaminase — 5-FC prodrug negative selection'),
        ]
        modules = []
        for name_pat, weight, rule_id, note in table:
            for idx, f in self._find_all(features, [name_pat]):
                # Require CDS feature type when available
                if (f.get('type') or '').lower() not in ('cds', '', 'misc_feature', 'gene'):
                    continue
                modules.append({
                    'module_type': 'counter_selection_module',
                    'start': f['start'], 'end': f['end'],
                    'strand': f.get('strand', 1),
                    'name': f['name'],
                    'features': [idx],
                    'weight': weight,
                    'detection_method': 'rule_based',
                    'rule_id': rule_id,
                    'notes': note
                })
        return modules

    def _detect_tet_regulator_cassette(self, features):
        """REG-TET-01: tTA / rtTA / Tet-On 3G regulator feature → regulator cassette module."""
        patterns = ['tTA', 'rtTA', 'Tet-On 3G', 'tet-on 3g', 'tTA-Advanced', 'rtTA-Advanced', 'rtTA3']
        modules = []
        seen = set()
        for idx, f in self._find_all(features, patterns):
            name_lc = (f.get('name') or '').lower()
            # Avoid double-matching TetR (0.8 weight, separate rule) when looking for tTA as substring
            if 'tetr' in name_lc and 'tta' not in name_lc and 'rtta' not in name_lc:
                continue
            if idx in seen:
                continue
            seen.add(idx)
            modules.append({
                'module_type': 'tet_regulator_cassette',
                'start': f['start'], 'end': f['end'],
                'strand': f.get('strand', 1),
                'name': f['name'],
                'features': [idx],
                'weight': 0.95,
                'detection_method': 'rule_based',
                'rule_id': 'REG-TET-01',
                'notes': 'TetR-VP16 transactivator; pairs with TRE-driven payload'
            })
        return modules

    def _detect_aid_degron_system(self, features):
        """REG-AID-01: OsTIR1 + (AID | mini-AID) on the same plasmid → AID degron system."""
        tir = self._find_all(features, ['OsTIR1', 'TIR1'])
        aid = self._find_all(features, ['mini-AID', 'miniAID', 'AID tag', 'AID*', 'AID degron'])
        if not tir or not aid:
            return []
        feat_idxs = [i for grp in (tir, aid) for i, _ in grp]
        feats_hit = [features[i] for i in feat_idxs]
        start = min(f['start'] for f in feats_hit)
        end = max(f['end'] for f in feats_hit)
        return [{
            'module_type': 'aid_degron_system',
            'start': start, 'end': end,
            'strand': 1,
            'name': 'AID degron system (OsTIR1 + AID tag)',
            'features': sorted(set(feat_idxs)),
            'weight': 0.95,
            'detection_method': 'rule_based',
            'rule_id': 'REG-AID-01',
            'notes': 'Auxin-inducible degradation of AID-tagged target'
        }]

    def _detect_fkbp_frb_dimerization(self, features):
        """REG-DIM-01: FKBP + FRB on the same plasmid → rapamycin-inducible dimerization."""
        fkbp = self._find_all(features, ['FKBP'])
        frb = self._find_all(features, ['FRB'])
        if not fkbp or not frb:
            return []
        feat_idxs = [i for grp in (fkbp, frb) for i, _ in grp]
        feats_hit = [features[i] for i in feat_idxs]
        start = min(f['start'] for f in feats_hit)
        end = max(f['end'] for f in feats_hit)
        return [{
            'module_type': 'fkbp_frb_dimerization',
            'start': start, 'end': end,
            'strand': 1,
            'name': 'FKBP/FRB rapamycin dimerization',
            'features': sorted(set(feat_idxs)),
            'weight': 0.9,
            'detection_method': 'rule_based',
            'rule_id': 'REG-DIM-01',
            'notes': 'Chemically inducible dimerization with rapamycin / rapalogs'
        }]

    # ------------------------------------------------------------------
    # LAC-BW: lacZα blue/white screening modules
    # ------------------------------------------------------------------
    # The lacZα fragment supports α-complementation in lacZΔM15 host strains
    # (DH5α, XL1-Blue, JM109): a functional β-galactosidase reassembles only
    # when the host's ω fragment binds an intact α fragment, giving blue
    # colonies on X-gal/IPTG plates. An MCS engineered into lacZα disrupts
    # the reading frame when an insert is cloned in, producing white
    # colonies — the classical blue/white screen.
    #
    # Five distinct architectures appear in the Module_Library_gb plasmids:
    #   1. puc_classic     (LAC-BW-01): lac promoter → lac operator → MCS
    #                                   inside lacZα. No phage promoters.
    #                                   Examples: pUC18/19, pUC57, pUC118.
    #   2. bluescript_t7t3 (LAC-BW-02): lacZα CDS containing T7 + MCS + T3
    #                                   phage promoters, plus upstream lac
    #                                   promoter/operator. Examples:
    #                                   pBluescript SK/KS, pBC SK/KS, pBS.
    #   3. pgem_t7sp6      (LAC-BW-03): same as #2 but T7 + MCS + SP6.
    #                                   Examples: pGEM-3Z, pGEM-5Zf, pSpark.
    #   4. litmus_t7_only  (LAC-BW-04): only T7 phage promoter present
    #                                   inside lacZα (no T3 and no SP6).
    #                                   Single-strand IVT + blue/white.
    #                                   Examples: LITMUS28/29/38/39,
    #                                   pEASY-T1, pEZ BAC.
    #   5. mammalian_bw    (LAC-BW-05): dual-phage MCS+lacZα retained inside
    #                                   a non-bacterial backbone (no nearby
    #                                   E. coli lac promoter). Used as a
    #                                   bacterial blue/white shuttle layer
    #                                   on mammalian vectors. Example:
    #                                   pBK-CMV.
    #
    # Each variant emits a parent `lac_alpha_blue_white_module` with a
    # `submodule_class` field naming the variant, plus a `submodules` list
    # decomposing the module into its components (lac_promoter,
    # lac_operator, mcs, lacZ_alpha CDS, phage promoter(s), M13 primers).

    # Substring patterns matched against feature names (case-insensitive).
    # Order matters only inside _classify_lac_variant.
    # Substring patterns (matched case-insensitively) that mark a CDS as
    # a candidate α fragment. Bare 'lacz' is permissive — it is filtered
    # by a ≤500 bp length cap in _is_lac_alpha_cds to exclude full-length
    # lacZ (~3 kb).
    _LAC_ALPHA_PATTERNS = (
        'lacz-alpha', 'lacz alpha', 'lacza', 'lacz\u03b1',
        'lacz fragment', 'lacz (fragment)', 'lacz (truncated)',
        'lacz a', 'lacz',
        # pLannotate's internal BLAST DB names the α-complementation fragment
        # `BGAL_SHISS` (β-galactosidase α, Shigella/short). Treat it like a
        # bare `lacz` — trust with the length cap.
        'bgal', 'bgal shiss', 'beta-galactosidase', 'beta galactosidase',
        '\u03b2-galactosidase',
    )
    _LAC_PROMOTER_PATTERNS = ('lac promoter', 'lacuv5 promoter', 'lac uv5')
    _LAC_OPERATOR_PATTERNS = ('lac operator',)
    _CAP_BINDING_PATTERNS = ('cap binding', 'cap site', 'cap-binding')
    _MCS_PATTERNS = ('mcs', 'multiple cloning', 'polylinker')
    _T7_PROMOTER_PATTERNS = ('t7 promoter',)
    _T3_PROMOTER_PATTERNS = ('t3 promoter',)
    _SP6_PROMOTER_PATTERNS = ('sp6 promoter',)
    _M13_FWD_PATTERNS = ('m13 fwd', 'm13 forward', 'm13/puc fwd', 'm13 puc fwd', 'm13(-20)', 'm13 -20')
    _M13_REV_PATTERNS = ('m13 rev', 'm13 reverse', 'm13/puc rev', 'm13 puc rev')
    _LAC_I_PATTERNS = ('laci', 'lac i ', 'lac repressor')

    # Name patterns that unambiguously mark a feature as the α fragment
    # (as opposed to full-length lacZ). When one of these hits, trust the
    # annotation even if the feature length looks large — circular
    # plasmids where the α CDS wraps the origin are often emitted with
    # start=0, end=plasmid_length by standard parsers.
    _LAC_ALPHA_EXPLICIT_PATTERNS = (
        'lacz-alpha', 'lacz alpha', 'lacz\u03b1', 'lacza',
        'lacz fragment', 'lacz (fragment)', 'lacz (truncated)',
    )

    @staticmethod
    def _is_lac_alpha_cds(feat: Dict[str, Any], plasmid_length: Optional[int] = None) -> bool:
        """True iff the feature is an α-fragment of lacZ (not full-length lacZ).

        Full-length lacZ (~3 kb, blue colonies on X-gal directly) is
        excluded — only the truncated α fragment (~110-360 bp window
        engineered for α-complementation) supports blue/white screening.
        """
        name = (feat.get('name') or '').lower()
        if not any(p in name for p in RuleBasedModuleDetector._LAC_ALPHA_PATTERNS):
            return False
        # If the name explicitly says alpha / fragment / truncated, trust
        # it regardless of length (covers origin-wrapping CDSes where the
        # parser reports a whole-plasmid span).
        if any(p in name for p in RuleBasedModuleDetector._LAC_ALPHA_EXPLICIT_PATTERNS):
            return True
        # For ambiguous names like bare 'lacz' or 'lacz a', gate on length.
        # α fragments in blue/white vectors are typically 100-400 bp; full-
        # length lacZ is ~3 kb and never participates in α-complementation.
        try:
            span = feat.get('end') - feat.get('start')
        except (TypeError, ValueError):
            span = None
        if span is None:
            return True
        # Reject features large enough to be full-length lacZ.
        if span > 500:
            # Allow for origin-wrap artefacts where the parser reports a
            # whole-plasmid span: if we know plasmid_length and span is
            # ≥80% of it, the feature is almost certainly an α fragment
            # that wraps the circular origin (seen in pGEM-style vectors).
            if plasmid_length and span >= 0.8 * plasmid_length:
                return True
            return False
        return True

    @staticmethod
    def _spans_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
        # Fall back to a safe "no overlap" when either interval is a wrap
        # (end < start). Proximity callers short-circuit wrap features
        # before reaching this helper.
        if a_end < a_start or b_end < b_start:
            return False
        return not (a_end < b_start or b_end < a_start)

    @staticmethod
    def _feature_is_wrap(feat: Dict[str, Any], plasmid_length: Optional[int]) -> bool:
        """True if the feature straddles the circular plasmid origin.

        Two conventions are observed in practice:
          * pLannotate emits `qstart > qend` when a hit crosses the origin
            (e.g. T7 at qstart=2983, qend=2 on a 3000 bp plasmid).
          * Biopython's SeqIO collapses `join(2983..3000, 1..2)` to
            `start=0, end=3000` — a feature whose span is the entire
            plasmid. We detect this as span ≥50% of the plasmid length.
        Wrap features have no reliable positional anchor, so proximity
        logic treats them as "anywhere on the plasmid".
        """
        try:
            start, end = feat.get('start'), feat.get('end')
        except AttributeError:
            return False
        if start is None or end is None:
            return False
        if end < start:
            return True
        if plasmid_length and (end - start) >= 0.5 * plasmid_length:
            return True
        return False

    @staticmethod
    def _tightest_circular_arc(intervals, plasmid_length: Optional[int]):
        """Find the tightest arc on a circular plasmid that covers every
        interval in `intervals` (list of (start, end)). Intervals that
        wrap the plasmid origin (end < start, or span ≥50% of plasmid)
        are skipped — they carry no reliable positional anchor.

        Returns (arc_start, arc_end). If the covering arc crosses the
        origin, `arc_end < arc_start`. Returns None if no valid interval.
        """
        if not intervals:
            return None
        pts = []
        for s, e in intervals:
            if s is None or e is None:
                continue
            if plasmid_length and e < s:
                continue  # pLannotate wrap (qstart>qend)
            if plasmid_length and (e - s) >= 0.5 * plasmid_length:
                continue  # Biopython whole-plasmid wrap
            pts.append((s, e))
        if not pts:
            return None
        if not plasmid_length or plasmid_length <= 0:
            return (min(s for s, _ in pts), max(e for _, e in pts))
        # Deduplicate by (start, end) and sort by start
        pts = sorted(set(pts), key=lambda t: t[0])
        n = len(pts)
        # Find the largest gap between consecutive intervals on the circle;
        # the module arc is the complement of that gap.
        max_gap = -1
        max_gap_idx = 0
        for i in range(n):
            nxt = (i + 1) % n
            gap = pts[nxt][0] - pts[i][1]
            if nxt == 0:
                gap += plasmid_length
            if gap < 0:
                gap = 0
            if gap > max_gap:
                max_gap = gap
                max_gap_idx = i
        arc_start = pts[(max_gap_idx + 1) % n][0]
        arc_end = pts[max_gap_idx][1]
        return (arc_start, arc_end)

    @staticmethod
    def _circular_gap(a_s: int, a_e: int, b_s: int, b_e: int,
                      plasmid_length: Optional[int]) -> int:
        """Minimum bp gap between intervals [a_s,a_e] and [b_s,b_e].

        If `plasmid_length` is given, the gap is computed on a circular
        plasmid — the lesser of the forward and backward arcs. Overlapping
        intervals return 0.
        """
        # Linear overlap short-circuit
        if a_e >= b_s and b_e >= a_s:
            return 0
        if not plasmid_length:
            return (b_s - a_e) if b_s > a_e else (a_s - b_e)
        # Forward arc: a's end → b's start (going in +1 direction)
        fwd = (b_s - a_e) % plasmid_length
        # Backward arc: b's end → a's start (same direction, from b)
        bwd = (a_s - b_e) % plasmid_length
        return min(fwd, bwd)

    @staticmethod
    def _spans_within(inner_s, inner_e, outer_s, outer_e, slack: int = 30) -> bool:
        """True if `inner` lies inside `outer` allowing a small flank slack (bp)."""
        return inner_s >= outer_s - slack and inner_e <= outer_e + slack

    def _detect_lac_blue_white_module(self, features, sequence: Optional[str] = None):
        """LAC-BW-01..05: lacZα blue/white screening modules (variant-aware).

        Produces one parent module per lacZα CDS that is paired with an MCS
        or at least one nearby phage promoter. The variant is selected by
        `_classify_lac_variant`. `sequence`, when given, enables circular-
        distance math so that features straddling the plasmid origin are
        still recognised as neighbours of an α CDS on the other side.
        """
        plasmid_length = len(sequence) if sequence else None

        # All α-fragment CDSes — de-duplicated by (start, end) because
        # AIP files often stack redundant misc_feature annotations.
        alpha_hits_all = [(i, f) for i, f in enumerate(features)
                          if self._is_lac_alpha_cds(f, plasmid_length)]
        seen_spans = set()
        alpha_hits = []
        for i, f in alpha_hits_all:
            key = (f['start'], f['end'])
            if key in seen_spans:
                continue
            seen_spans.add(key)
            alpha_hits.append((i, f))
        if not alpha_hits:
            return []

        # Index helper sets once
        mcs_hits = self._find_all(features, list(self._MCS_PATTERNS))
        lac_prom = self._find_all(features, list(self._LAC_PROMOTER_PATTERNS))
        lac_op = self._find_all(features, list(self._LAC_OPERATOR_PATTERNS))
        cap_hits = self._find_all(features, list(self._CAP_BINDING_PATTERNS))
        t7_hits = self._find_all(features, list(self._T7_PROMOTER_PATTERNS))
        t3_hits = self._find_all(features, list(self._T3_PROMOTER_PATTERNS))
        sp6_hits = self._find_all(features, list(self._SP6_PROMOTER_PATTERNS))
        m13f_hits = self._find_all(features, list(self._M13_FWD_PATTERNS))
        m13r_hits = self._find_all(features, list(self._M13_REV_PATTERNS))
        laci_hits = self._find_all(features, list(self._LAC_I_PATTERNS))

        modules = []

        def _near(group, max_gap):
            """Filter `group` to entries within `max_gap` bp of the α CDS.

            Wrap-features (span >50% of plasmid) are accepted
            unconditionally — they straddle the origin and their raw
            start/end cannot be used for distance math. Overlapping
            features are always accepted.
            """
            out = []
            for idx, f in group:
                if self._feature_is_wrap(f, plasmid_length):
                    out.append((idx, f))
                    continue
                if self._spans_overlap(astart, aend, f['start'], f['end']):
                    out.append((idx, f))
                    continue
                gap = self._circular_gap(astart, aend, f['start'], f['end'], plasmid_length)
                if gap <= max_gap:
                    out.append((idx, f))
            return out

        for ai, af in alpha_hits:
            astart, aend = af['start'], af['end']

            # MCS within 120 bp of α (spanning or flanking it)
            mcs_in = _near(mcs_hits, 120)

            # Phage promoters within 250 bp of α — pBSK-style vectors
            # annotate a narrow α fragment separate from the MCS+promoter
            # block, so looser slack is needed than a strict "inside α" test.
            t7_in = _near(t7_hits, 250)
            t3_in = _near(t3_hits, 250)
            sp6_in = _near(sp6_hits, 250)

            # Need either an MCS or at least one phage promoter near the
            # α CDS, otherwise this is just a fragmentary annotation.
            if not mcs_in and not (t7_in or t3_in or sp6_in):
                continue

            # lac promoter / operator within ~600 bp of α (native lac
            # upstream region is ~120 bp; extra slack covers mammalian
            # shuttles where only the operator is retained)
            lac_prom_in = _near(lac_prom, 600)
            lac_op_in = _near(lac_op, 600)
            cap_in = _near(cap_hits, 600)

            # M13 primer sites: typically flank the α CDS within 50 bp
            m13f_in = _near(m13f_hits, 200)
            m13r_in = _near(m13r_hits, 200)

            # lacI lies ~100-400 bp downstream of the lac promoter and
            # is part of the blue/white regulatory block — include it so
            # the emitted module spans the full lac cassette.
            laci_in = _near(laci_hits, 800)

            variant, rule_id, weight = self._classify_lac_variant(
                lac_prom_in, t7_in, t3_in, sp6_in,
            )
            if variant is None:
                continue

            # Composite span: compute the tightest circular arc covering
            # every supporting feature. Wrap-features (e.g. pLannotate
            # T7 spanning the origin) are excluded from arc math because
            # they carry no reliable positional anchor; they still appear
            # in the submodule list with their raw coordinates.
            all_feat_idxs = {ai}
            spans = [(astart, aend)]
            for grp in (mcs_in, lac_prom_in, lac_op_in, cap_in,
                        t7_in, t3_in, sp6_in, m13f_in, m13r_in, laci_in):
                for idx, f in grp:
                    all_feat_idxs.add(idx)
                    spans.append((f['start'], f['end']))
            arc = self._tightest_circular_arc(spans, plasmid_length)
            if arc is None:
                # Only wrap-features were available for span math — fall
                # back to the α CDS bounds.
                mod_start, mod_end = astart, aend
            else:
                mod_start, mod_end = arc

            submodules = self._build_lac_submodules(
                features, ai, mcs_in, lac_prom_in, lac_op_in, cap_in,
                t7_in, t3_in, sp6_in, m13f_in, m13r_in, laci_in,
            )

            phage = []
            if t7_in:
                phage.append('T7')
            if t3_in:
                phage.append('T3')
            if sp6_in:
                phage.append('SP6')
            phage_str = '+'.join(phage) if phage else 'no phage promoters'

            modules.append({
                'module_type': 'lac_alpha_blue_white_module',
                'submodule_class': variant,
                'start': mod_start,
                'end': mod_end,
                'strand': af.get('strand', 1),
                'name': f'lacZα blue/white module ({variant}, {phage_str})',
                'features': sorted(all_feat_idxs),
                'submodules': submodules,
                'weight': weight,
                'detection_method': 'rule_based',
                'rule_id': rule_id,
                'notes': (
                    'α-complementation cassette: insert into MCS disrupts '
                    'lacZα reading frame → white colony on X-gal/IPTG. '
                    f'Variant: {variant}.'
                ),
            })

        return modules

    def _detect_lac_alpha_disrupted_module(self, features, sequence=None):
        """LAC-DISRUPT-01: lac-α disrupted by insert.

        Pattern: lac_operator + MCS (>=1) + (T7 | T3 | SP6 | lac_promoter)
                 are all present within ~1000 bp of each other,
                 BUT no intact lacZα CDS is detected.

        This is the canonical "successful blue/white insert" plasmid:
        the foreign insert (e.g. Cas9 in pBS-Hsp70-Cas9, a PCR product
        in pPCR-Script with an insert) replaced the α-fragment and the
        colony screens white.

        Emits one module spanning the lac-element cluster's bounding
        box with metadata.disruption_type = "insert" and a confidence
        derived from how many of the canonical pBSK elements survive.
        """
        plasmid_length = len(sequence) if sequence else None

        # If an intact lacZα CDS is detected, the regular lac-BW rule
        # handles it; skip emitting a disrupted module.
        for f in features:
            if self._is_lac_alpha_cds(f, plasmid_length):
                return []

        mcs_hits = self._find_all(features, list(self._MCS_PATTERNS))
        lac_op = self._find_all(features, list(self._LAC_OPERATOR_PATTERNS))
        lac_prom = self._find_all(features, list(self._LAC_PROMOTER_PATTERNS))
        t7_hits = self._find_all(features, list(self._T7_PROMOTER_PATTERNS))
        t3_hits = self._find_all(features, list(self._T3_PROMOTER_PATTERNS))
        sp6_hits = self._find_all(features, list(self._SP6_PROMOTER_PATTERNS))
        cap_hits = self._find_all(features, list(self._CAP_BINDING_PATTERNS))

        # Require at minimum: lac operator, an MCS, and at least one phage
        # or lac promoter — the signature of a pBluescript-derived cloning
        # vector with an insert.
        if not lac_op or not mcs_hits:
            return []
        phage_or_lac = bool(lac_prom or t7_hits or t3_hits or sp6_hits)
        if not phage_or_lac:
            return []

        # Bounding box covers every surviving lac-α anchor. No span
        # filter — an insert that disrupts lacZα can be arbitrarily
        # large, and the signature we want is exactly "lac apparatus
        # present without an intact lacZα CDS," not cluster tightness.
        all_anchors = [f for _, f in lac_op + mcs_hits + lac_prom + t7_hits + t3_hits + sp6_hits]
        starts = [int(f.get("start", 0)) for f in all_anchors]
        ends = [int(f.get("end", 0)) for f in all_anchors]
        cluster_start = min(starts)
        cluster_end = max(ends)

        # Build submodule list — every surviving pBSK element gets
        # captured so the user can see WHAT remains after the insert.
        submods = []
        def _add_submod(group, label, mtype):
            for _, f in group:
                submods.append({
                    "module_type": mtype,
                    "start": int(f.get("start", 0)),
                    "end": int(f.get("end", 0)),
                    "strand": f.get("strand", 1),
                    "name": f.get("name", label),
                })
        _add_submod(lac_op, "lac operator", "lac_operator")
        _add_submod(mcs_hits, "MCS", "multiple_cloning_site")
        _add_submod(lac_prom, "lac promoter", "lac_promoter")
        _add_submod(cap_hits, "CAP binding site", "cap_binding_site")
        _add_submod(t7_hits, "T7 promoter", "t7_promoter")
        _add_submod(t3_hits, "T3 promoter", "t3_promoter")
        _add_submod(sp6_hits, "SP6 promoter", "sp6_promoter")

        # Confidence — 5/5 elements (op, MCS, lac prom, T7, T3) = 0.95;
        # 3/5 = 0.80; 2/5 (just op+MCS) = 0.70.
        present = sum(1 for g in (lac_op, mcs_hits, lac_prom,
                                   t7_hits or t3_hits or sp6_hits,
                                   cap_hits) if g)
        weight = 0.60 + 0.08 * present  # 5 present -> 1.00; 2 -> 0.76

        return [{
            "module_type": "lac_alpha_disrupted_module",
            "start": cluster_start,
            "end": cluster_end,
            "name": "lacα disrupted (cloning insert)",
            "features": [],
            "submodules": submods,
            "weight": min(weight, 0.95),
            "detection_method": "rule_based",
            "rule_id": "LAC-DISRUPT-01",
            "metadata": {
                "disruption_type": "insert",
                "element_count": present,
                "n_mcs": len(mcs_hits),
                "n_lac_op": len(lac_op),
                "n_lac_prom": len(lac_prom),
                "phage_promoters": len(t7_hits) + len(t3_hits) + len(sp6_hits),
            },
            "notes": (
                f"lacα CDS missing; surviving pBSK elements: "
                f"{len(mcs_hits)} MCS, {len(lac_op)} lac op, "
                f"{len(lac_prom)} lac prom, "
                f"{len(t7_hits) + len(t3_hits) + len(sp6_hits)} phage prom. "
                f"Indicates a successful insert into the multiple cloning site "
                f"(blue/white screening — colony screens white)."
            ),
        }]

    @staticmethod
    def _classify_lac_variant(lac_prom_in, t7_in, t3_in, sp6_in):
        """Pick the LAC-BW-* variant rule_id from supporting features.

        Returns (variant_name, rule_id, weight) or (None, None, None) if
        no variant fits.
        """
        has_lac_prom = bool(lac_prom_in)
        has_t7 = bool(t7_in)
        has_t3 = bool(t3_in)
        has_sp6 = bool(sp6_in)

        # Prefer the most specific dual-phage classifications first.
        if has_t7 and has_t3:
            if has_lac_prom:
                return 'bluescript_t7t3', 'LAC-BW-02', 0.97
            return 'mammalian_bw', 'LAC-BW-05', 0.92
        if has_t7 and has_sp6:
            if has_lac_prom:
                return 'pgem_t7sp6', 'LAC-BW-03', 0.97
            return 'mammalian_bw', 'LAC-BW-05', 0.92
        if has_t7 and not (has_t3 or has_sp6):
            # T7-only variant (LITMUS, pEASY-T1, pEZ BAC) — a single
            # phage promoter lets you run IVT off one strand only.
            if has_lac_prom:
                return 'litmus_t7_only', 'LAC-BW-04', 0.93
            return 'mammalian_bw', 'LAC-BW-05', 0.9
        if has_sp6 and not has_t7:
            # SP6-only variant — SP6 inside α CDS without T7. No current
            # library example, but the architecture is plausible.
            if has_lac_prom:
                return 'pgem_sp6_only', 'LAC-BW-04', 0.9
            return 'mammalian_bw', 'LAC-BW-05', 0.88
        if has_lac_prom and not (has_t7 or has_t3 or has_sp6):
            return 'puc_classic', 'LAC-BW-01', 0.95
        # Fallback: lacZα + MCS but no phage promoters and no annotated
        # lac promoter. Most often a partially-annotated pUC derivative.
        return 'puc_classic', 'LAC-BW-01', 0.85

    @staticmethod
    def _build_lac_submodules(
        features, alpha_idx, mcs_in, lac_prom_in, lac_op_in, cap_in,
        t7_in, t3_in, sp6_in, m13f_in, m13r_in, laci_in,
    ):
        """Decompose the parent module into per-component submodules.

        Each submodule is a dict with keys analogous to the parent module
        (`module_type`, `start`, `end`, `strand`, `name`, `feature_index`,
        `rule_id`). Downstream consumers (Step 5 hierarchical conversion)
        can render these as nested annotations under the parent.
        """
        af = features[alpha_idx]
        subs = [{
            'module_type': 'lac_alpha_cds',
            'start': af['start'], 'end': af['end'],
            'strand': af.get('strand', 1),
            'name': af.get('name', 'lacZα'),
            'feature_index': alpha_idx,
            'rule_id': 'LAC-BW-SUB-ALPHA',
            'notes': 'α-fragment of lacZ supporting α-complementation',
        }]

        def _push(group, mod_type, rule_id, note):
            for idx, f in group:
                subs.append({
                    'module_type': mod_type,
                    'start': f['start'], 'end': f['end'],
                    'strand': f.get('strand', 1),
                    'name': f.get('name', mod_type),
                    'feature_index': idx,
                    'rule_id': rule_id,
                    'notes': note,
                })

        _push(lac_prom_in, 'lac_promoter', 'LAC-BW-SUB-PROMOTER',
              'σ70 lac promoter — repressed by LacI, induced by IPTG')
        _push(cap_in, 'cap_binding_site', 'LAC-BW-SUB-CAP',
              'CAP/cAMP binding site — catabolite activation of lac promoter')
        _push(lac_op_in, 'lac_operator', 'LAC-BW-SUB-OPERATOR',
              'LacI binding site — IPTG relieves repression')
        _push(mcs_in, 'mcs', 'LAC-BW-SUB-MCS',
              'Multiple cloning site embedded in lacZα; insertion disrupts α-complementation')
        _push(t7_in, 't7_phage_promoter', 'LAC-BW-SUB-T7',
              'T7 RNA polymerase promoter for in vitro transcription of insert')
        _push(t3_in, 't3_phage_promoter', 'LAC-BW-SUB-T3',
              'T3 RNA polymerase promoter for in vitro transcription of insert (opposite strand)')
        _push(sp6_in, 'sp6_phage_promoter', 'LAC-BW-SUB-SP6',
              'SP6 RNA polymerase promoter for in vitro transcription of insert')
        _push(m13f_in, 'm13_fwd_primer', 'LAC-BW-SUB-M13F',
              'M13 forward sequencing primer site flanking MCS')
        _push(m13r_in, 'm13_rev_primer', 'LAC-BW-SUB-M13R',
              'M13 reverse sequencing primer site flanking MCS')
        _push(laci_in, 'lac_i_gene', 'LAC-BW-SUB-LACI',
              'LacI repressor — binds lac operator in absence of IPTG, maintaining blue/white gating')

        subs.sort(key=lambda s: (s['start'], s['end']))
        return subs

    def _detect_ori_by_kb_class(self, features: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Detect replication origin modules based on KB class.

        Any feature with KB class "replication_origin" or "origin" becomes an Ori module.

        EBV oriP matches are gated on EBNA1 co-presence — when EBNA1 is absent,
        the bare "oriP" / "EBV" name hit is skipped here and handled by the
        specialised _detect_ebv_episomal (MR-EBV-01) if applicable. This prevents
        stray partial-name hits (e.g. inside WPRE regions) from being mis-labelled
        as mammalian_replication.
        """
        modules = []

        # Co-presence check: is EBNA1 annotated anywhere on the plasmid?
        has_ebna1 = any(
            ('ebna1' in (f.get('name') or '').lower().replace('-', '').replace(' ', '')
             or 'ebna 1' in (f.get('name') or '').lower()
             or 'ebna-1' in (f.get('name') or '').lower())
            for f in features
        )

        for i, feat in enumerate(features):
            kb_class = feat.get('kb_class', '').lower()
            name = feat.get('name', '').lower()

            # Check if it's a replication origin by KB class or name pattern
            is_ori_by_kb = kb_class in {'replication_origin', 'origin', 'ori'}
            is_ori_by_name = 'ori' in name and not any(x in name for x in ['origin of transfer', 'orit'])

            if is_ori_by_kb or is_ori_by_name:
                # Determine specific ori type based on name
                if '2μ' in name or '2micron' in name or '2-micron' in name or '2u' in name or '2 u' in name:
                    module_type = 'yeast_replication'
                elif 'cen' in name or 'ars' in name:
                    module_type = 'yeast_replication'
                elif 'f1' in name or 'm13' in name:
                    module_type = 'phage_replication'
                elif 'sv40' in name:
                    module_type = 'mammalian_replication'
                elif 'ebv' in name or 'orip' in name:
                    # Stricter gate: only emit EBV mammalian_replication when EBNA1
                    # is present elsewhere. Otherwise the composite MR-EBV-01 rule
                    # (which requires oriP + EBNA1 co-presence) should handle it,
                    # or the feature is a partial-name false positive.
                    if not has_ebna1:
                        continue
                    module_type = 'mammalian_replication'
                else:
                    module_type = 'bacterial_replication'

                module = {
                    'module_type': module_type,
                    'start': feat['start'],
                    'end': feat['end'],
                    'strand': feat.get('strand', 1),
                    'name': feat.get('name', 'ori'),
                    'features': [i],
                    'weight': 0.95,
                    'detection_method': 'rule_based',
                    'notes': f'Detected by KB class: {kb_class}' if is_ori_by_kb else 'Detected by name pattern'
                }
                modules.append(module)

        return modules


    def _detect_gateway_cassettes(self, features: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Detect Gateway cloning cassettes based on att site pairs.

        Gateway cloning uses paired att sites that flank a cargo region:
        - attP1/attP2: Donor vectors (pDONR) - att sites face inward
        - attB1/attB2: PCR products/expression clones - att sites face inward
        - attL1/attL2: Entry clones - att sites face inward
        - attR1/attR2: Destination vectors - att sites face inward (usually with ccdB)

        The numbered sites (1 and 2) must match for a valid cassette.
        Site 1 should be 5-prime of site 2 on the plasmid.
        """
        modules = []

        # Categorize att sites by type and number
        att_sites = {
            'attP': {'1': [], '2': [], 'other': []},
            'attB': {'1': [], '2': [], 'other': []},
            'attL': {'1': [], '2': [], 'other': []},
            'attR': {'1': [], '2': [], 'other': []},
        }

        import re
        att_pattern = re.compile(r'att([PBLR])([1-6])?', re.IGNORECASE)

        for i, feat in enumerate(features):
            name = feat.get('name', '')
            match = att_pattern.search(name)
            if match:
                att_type = 'att' + match.group(1).upper()
                att_num = match.group(2) if match.group(2) else 'other'
                if att_type in att_sites:
                    att_sites[att_type][att_num].append({
                        'index': i,
                        'feat': feat,
                        'start': feat['start'],
                        'end': feat['end'],
                        'strand': feat.get('strand', 1),
                    })

        # Detect cassettes for each att type
        cassette_types = {
            'attP': ('gateway_donor_cassette', 'Gateway Donor Cassette (attP)'),
            'attB': ('gateway_entry_cassette', 'Gateway Entry Cassette (attB)'),
            'attL': ('gateway_entry_cassette', 'Gateway Entry Clone (attL)'),
            'attR': ('gateway_destination_cassette', 'Gateway Destination Cassette (attR)'),
        }

        for att_type, (module_type, module_name) in cassette_types.items():
            sites = att_sites[att_type]

            # Match site 1 with site 2
            for site1 in sites['1']:
                for site2 in sites['2']:
                    # Site 1 should be before site 2 in sequence
                    if site1['start'] < site2['start']:
                        # Create cassette module spanning from site 1 to site 2
                        cassette_start = site1['start']
                        cassette_end = site2['end']

                        module = {
                            'module_type': module_type,
                            'start': cassette_start,
                            'end': cassette_end,
                            'strand': 1,
                            'name': module_name,
                            'features': [site1['index'], site2['index']],
                            'weight': 0.95,
                            'detection_method': 'rule_based',
                            'notes': 'Bounded by ' + site1['feat']['name'] + ' and ' + site2['feat']['name'],
                            'att_site_1': site1['feat']['name'],
                            'att_site_2': site2['feat']['name'],
                        }
                        modules.append(module)

            # Also check for higher-numbered pairs (multisite Gateway: 3-4, 5-6)
            for num1, num2 in [('3', '4'), ('5', '6')]:
                for site1 in sites.get(num1, []):
                    for site2 in sites.get(num2, []):
                        if site1['start'] < site2['start']:
                            module = {
                                'module_type': 'multisite_' + module_type,
                                'start': site1['start'],
                                'end': site2['end'],
                                'strand': 1,
                                'name': 'MultiSite ' + module_name + ' (' + num1 + '-' + num2 + ')',
                                'features': [site1['index'], site2['index']],
                                'weight': 0.90,
                                'detection_method': 'rule_based',
                                'notes': 'Bounded by ' + site1['feat']['name'] + ' and ' + site2['feat']['name'],
                            }
                            modules.append(module)

        return modules


    def _detect_floxed_regions(self, features: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Detect floxed (loxP-flanked) regions for Cre-mediated recombination.

        Common patterns:
        - LSL (loxP-STOP-loxP): Two loxP sites flanking stop cassettes, both on same strand
          with 5' ends pointing inward toward the cargo
        - loxP-flanked genes: Similar pattern for conditional gene deletion

        loxP sites should be on the same strand and facing inward:
        - Forward strand loxP on 5' side
        - Forward strand loxP on 3' side (or reverse complement facing the same direction)

        The key is that both loxP sites face toward the internal region for recombination.
        """
        modules = []

        # Find all loxP sites
        loxp_sites = []
        for i, feat in enumerate(features):
            name = feat.get('name', '').lower()
            # Match loxP, lox, lox2272, loxN, etc.
            if 'loxp' in name or 'lox2272' in name or 'lox66' in name or 'lox71' in name or 'loxn' in name:
                loxp_sites.append({
                    'index': i,
                    'feat': feat,
                    'start': feat['start'],
                    'end': feat['end'],
                    'strand': feat.get('strand', 1),
                    'name': feat.get('name', 'loxP'),
                })

        if len(loxp_sites) < 2:
            return modules

        # Sort by position
        loxp_sites.sort(key=lambda x: x['start'])

        # Find pairs of loxP sites that flank a region
        # For a valid floxed cassette:
        # - Two loxP sites on the same strand (both +1 or both -1)
        # - Or one forward and one reverse, both "pointing" toward the internal region
        used_indices = set()

        for i, site1 in enumerate(loxp_sites):
            if i in used_indices:
                continue

            for j, site2 in enumerate(loxp_sites):
                if j <= i or j in used_indices:
                    continue

                # Check if they form a valid floxed pair
                # Sites should be on the same strand for canonical LSL cassettes
                # OR in orientations that point toward each other
                strand1 = site1['strand']
                strand2 = site2['strand']

                # Calculate distance between sites
                distance = site2['start'] - site1['end']

                # Skip if too close (< 50bp) or too far (> 10kb for typical cassettes)
                if distance < 50 or distance > 10000:
                    continue

                # Valid orientations for floxed regions:
                # 1. Both forward strand (common for LSL)
                # 2. Both reverse strand
                # 3. Site1 forward, Site2 reverse (facing each other)
                is_valid_pair = (
                    (strand1 == strand2) or  # Same strand
                    (strand1 == 1 and strand2 == -1)  # Facing each other
                )

                if is_valid_pair:
                    # Create floxed module
                    cassette_start = site1['start']
                    cassette_end = site2['end']

                    # Determine cassette type based on contents
                    # Check for STOP elements between the loxP sites
                    has_stop = False
                    internal_features = []
                    for feat in features:
                        feat_start = feat['start']
                        feat_end = feat['end']
                        feat_name = feat.get('name', '').lower()

                        # Check if feature is between the loxP sites
                        if feat_start >= site1['end'] and feat_end <= site2['start']:
                            internal_features.append(feat)
                            # Check for stop/terminator elements
                            if any(term in feat_name for term in ['stop', 'terminator', 'poly(a)', 'polya', 'sv40']):
                                has_stop = True

                    if has_stop:
                        module_type = 'lsl_cassette'
                        module_name = 'LSL (loxP-STOP-loxP) Cassette'
                    else:
                        module_type = 'floxed_region'
                        module_name = 'Floxed Region (loxP-flanked)'

                    module = {
                        'module_type': module_type,
                        'start': cassette_start,
                        'end': cassette_end,
                        'strand': strand1,  # Use first site's strand
                        'name': module_name,
                        'features': [site1['index'], site2['index']],
                        'weight': 0.93,
                        'detection_method': 'rule_based',
                        'notes': f"Bounded by {site1['name']} and {site2['name']}",
                        'loxp_site_1': site1['name'],
                        'loxp_site_2': site2['name'],
                        'internal_feature_count': len(internal_features),
                    }
                    modules.append(module)

                    # Mark these sites as used
                    used_indices.add(i)
                    used_indices.add(j)
                    break  # Move to next site1

        return modules


    def _detect_frt_flanked_regions(self, features: List[Dict[str, Any]], sequence: str = None) -> List[Dict[str, Any]]:
        """
        Detect FRT-flanked regions for Flp recombination.

        FRT (Flp Recombination Target) sites are recognized by Flp recombinase.
        Two FRT sites flanking a region allow Flp-mediated excision or inversion.

        If sequence is provided and no FRT features are found, performs motif scanning.
        """
        modules = []

        # Find all FRT sites from annotations
        # Note: FRT sites are often found within 2μ ori features due to sequence similarity
        frt_sites = []
        for i, feat in enumerate(features):
            name = feat.get('name', '').lower()
            # Match FRT, FRT3, FRT5, etc.
            if 'frt' in name and 'flirt' not in name:
                frt_sites.append({
                    'index': i,
                    'feat': feat,
                    'start': feat['start'],
                    'end': feat['end'],
                    'strand': feat.get('strand', 1),
                    'name': feat.get('name', 'FRT'),
                    'from_motif_scan': False,
                    'is_2u_ori': False,
                })
            # Also check for 2μ ori features which often contain FRT sites
            elif '2u' in name.replace('μ', 'u').replace('µ', 'u') or '2-u' in name or '2micro' in name or '2 micro' in name:
                # 2μ ori features often contain FRT sequences
                frt_sites.append({
                    'index': i,
                    'feat': feat,
                    'start': feat['start'],
                    'end': feat['end'],
                    'strand': feat.get('strand', 1),
                    'name': feat.get('name', '2μ ori (FRT)'),
                    'from_motif_scan': False,
                    'is_2u_ori': True,
                })

        # If no FRT sites found from annotations but sequence is provided, scan for motifs
        if len(frt_sites) < 2 and sequence:
            motif_frt_sites = self._scan_for_frt_motifs(sequence)
            for i, site in enumerate(motif_frt_sites):
                # Avoid adding duplicate positions
                is_duplicate = any(
                    abs(site['start'] - existing['start']) < 10
                    for existing in frt_sites
                )
                if not is_duplicate:
                    frt_sites.append({
                        'index': -1 - i,  # Negative index for motif-detected sites
                        'feat': None,
                        'start': site['start'],
                        'end': site['end'],
                        'strand': site['strand'],
                        'name': f"{site['name']} (motif)",
                        'from_motif_scan': True,
                    })

        print(f"[DEBUG FRT] Found {len(frt_sites)} FRT/2μ ori sites", flush=True)
        for site in frt_sites:
            print(f"[DEBUG FRT]   - {site['name']} at {site['start']}-{site['end']}, strand={site['strand']}, is_2u_ori={site.get('is_2u_ori', False)}", flush=True)

        if len(frt_sites) < 2:
            print(f"[DEBUG FRT] Less than 2 sites, returning empty", flush=True)
            return modules

        # Sort by position
        frt_sites.sort(key=lambda x: x['start'])

        # Find pairs
        used_indices = set()
        for i, site1 in enumerate(frt_sites):
            if i in used_indices:
                continue

            for j, site2 in enumerate(frt_sites):
                if j <= i or j in used_indices:
                    continue

                distance = site2['start'] - site1['end']
                if distance < 50 or distance > 15000:
                    continue

                # FRT sites typically on same strand for excision
                strand1 = site1['strand']
                strand2 = site2['strand']

                print(f"[DEBUG FRT] Checking pair: {site1['name']} (strand={strand1}) - {site2['name']} (strand={strand2}), distance={distance}", flush=True)
                if strand1 == strand2 or (strand1 == 1 and strand2 == -1):
                    print(f"[DEBUG FRT] Creating FRT module!", flush=True)
                    module = {
                        'module_type': 'frt_flanked_region',
                        'start': site1['start'],
                        'end': site2['end'],
                        'strand': strand1,
                        'name': 'FRT-Flanked Region (Flp)',
                        'features': [site1['index'], site2['index']],
                        'weight': 0.93,
                        'detection_method': 'rule_based',
                        'notes': f"Bounded by {site1['name']} and {site2['name']}",
                        'frt_site_1': site1['name'],
                        'frt_site_2': site2['name'],
                    }
                    modules.append(module)
                    used_indices.add(i)
                    used_indices.add(j)
                    break

        return modules

    def _detect_transposon_flanked_regions(self, features: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Detect transposon ITR-flanked regions.

        Supported transposon systems:
        - Sleeping Beauty (SB): IR/DR elements, SB ITR L/R
        - PiggyBac (PB): PB ITR 5'/3', TTAA target sites
        - Tn7: Tn7L/Tn7R (attTn7)
        - Tn5: Tn5 mosaic ends (ME)
        """
        modules = []

        # Categorize transposon elements by system
        transposon_sites = {
            'sleeping_beauty': [],
            'piggybac': [],
            'tn7': [],
            'tn5': [],
        }

        for i, feat in enumerate(features):
            name = feat.get('name', '').lower()
            site_info = {
                'index': i,
                'feat': feat,
                'start': feat['start'],
                'end': feat['end'],
                'strand': feat.get('strand', 1),
                'name': feat.get('name', ''),
                'is_left': False,
                'is_right': False,
            }

            # Sleeping Beauty ITRs
            if 'sleeping beauty' in name or 'sb itr' in name or 'sb100' in name or ('ir/dr' in name and 'sb' in name):
                site_info['is_left'] = 'left' in name or '_l' in name or '5' in name or 'ir/dr(l)' in name
                site_info['is_right'] = 'right' in name or '_r' in name or '3' in name or 'ir/dr(r)' in name
                transposon_sites['sleeping_beauty'].append(site_info)

            # PiggyBac ITRs
            elif 'piggybac' in name or 'pb itr' in name or 'pb 5' in name or 'pb 3' in name or 'pb-' in name:
                site_info['is_left'] = '5' in name or 'left' in name or '_l' in name
                site_info['is_right'] = '3' in name or 'right' in name or '_r' in name
                transposon_sites['piggybac'].append(site_info)

            # Tn7 - Tn7L and Tn7R define the transposon boundaries
            # In pFastBac vectors: Tn7R is upstream, Tn7L is downstream
            elif 'tn7' in name:
                # Tn7L = left end of transposon (downstream in pFastBac)
                site_info['is_left'] = 'tn7l' in name or 'tn7 l' in name or 'left' in name
                # Tn7R = right end of transposon (upstream in pFastBac)
                site_info['is_right'] = 'tn7r' in name or 'tn7 r' in name or 'right' in name
                # Note: In pFastBac, Tn7R (position ~2510) is before Tn7L (position ~4428)
                # The cargo is between them
                transposon_sites['tn7'].append(site_info)

            # Tn5 mosaic ends
            elif 'tn5' in name or 'mosaic end' in name:
                # Tn5 has identical ends, so just track them
                transposon_sites['tn5'].append(site_info)

        # Detect flanked regions for each system
        system_names = {
            'sleeping_beauty': ('sleeping_beauty_payload', 'Sleeping Beauty Transposon Payload'),
            'piggybac': ('piggybac_payload', 'PiggyBac Transposon Payload'),
            'tn7': ('tn7_payload', 'Tn7 Transposon Payload'),
            'tn5': ('tn5_payload', 'Tn5 Transposon Payload'),
        }

        for system, sites in transposon_sites.items():
            if len(sites) < 2:
                continue

            sites.sort(key=lambda x: x['start'])
            module_type, module_name = system_names[system]

            # For systems with L/R designation, match L with R
            left_sites = [s for s in sites if s['is_left']]
            right_sites = [s for s in sites if s['is_right']]

            if left_sites and right_sites:
                # Match sites - for Tn7, the R site often appears before L in sequence
                # Try both orientations
                for site_a in left_sites + right_sites:
                    for site_b in left_sites + right_sites:
                        if site_a['index'] == site_b['index']:
                            continue
                        # Ensure different site types for Tn7 (L pairs with R)
                        if system == 'tn7':
                            if site_a['is_left'] == site_b['is_left']:
                                continue  # Both L or both R, skip
                        if site_b['start'] > site_a['end']:
                            distance = site_b['start'] - site_a['end']
                            if 50 <= distance <= 20000:
                                module = {
                                    'module_type': module_type,
                                    'start': site_a['start'],
                                    'end': site_b['end'],
                                    'strand': 1,
                                    'name': module_name,
                                    'features': [site_a['index'], site_b['index']],
                                    'weight': 0.94,
                                    'detection_method': 'rule_based',
                                    'notes': f"Bounded by {site_a['name']} and {site_b['name']}",
                                }
                                modules.append(module)
                                break
                    else:
                        continue
                    break
            else:
                # For Tn5 or unmarked ITRs, just pair first two
                if len(sites) >= 2:
                    site1, site2 = sites[0], sites[1]
                    distance = site2['start'] - site1['end']
                    if 50 <= distance <= 20000:
                        module = {
                            'module_type': module_type,
                            'start': site1['start'],
                            'end': site2['end'],
                            'strand': 1,
                            'name': module_name,
                            'features': [site1['index'], site2['index']],
                            'weight': 0.92,
                            'detection_method': 'rule_based',
                            'notes': f"Bounded by {site1['name']} and {site2['name']}",
                        }
                        modules.append(module)

        return modules

    def _detect_baculovirus_homology_regions(self, features: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Detect baculovirus homology regions for recombination-based cloning.

        Baculovirus transfer vectors use two main homology regions:
        - LEF2/ORF603 locus (left homology arm)
        - ORF1629 locus (right homology arm)

        The expression cassette (polyhedrin/p10 promoter, gene, poly(A)) sits between them.
        """
        modules = []
        print(f'[DEBUG _detect_baculovirus] Called with {len(features)} features')

        # Debug: Print all feature names
        for i, feat in enumerate(features):
            fname = feat.get('name', '').lower()
            if any(x in fname for x in ['p10', 'polyhedrin', 'polh', 'baculovirus', 'baculo', 'orf1629', 'orf603', 'lef']):
                print(f'[DEBUG _detect_baculovirus] MATCH feature {i}: {feat.get("name", "")}')

        # Find baculovirus-related features and categorize by locus
        homology_regions = []
        for i, feat in enumerate(features):
            name = feat.get('name', '').lower()

            # Match baculovirus features
            if any(x in name for x in ['polyhedrin', 'polh', 'p10', 'baculovirus', 'baculo', 'orf1629', 'orf603', 'lef2', 'lef-2', 'recombination region']):
                # Determine which locus this belongs to
                # ORF1629 locus takes precedence if mentioned
                if 'orf1629' in name:
                    locus = 'orf1629'
                elif 'lef2' in name or 'lef-2' in name or 'orf603' in name:
                    locus = 'lef2_orf603'
                elif 'polyhedrin' in name or 'polh' in name or 'p10' in name:
                    locus = 'expression'  # These are in the expression cassette, not homology arms
                else:
                    locus = 'other'

                homology_regions.append({
                    'index': i,
                    'feat': feat,
                    'start': feat['start'],
                    'end': feat['end'],
                    'strand': feat.get('strand', 1),
                    'name': feat.get('name', ''),
                    'locus': locus,
                })

        if len(homology_regions) < 2:
            return modules

        homology_regions.sort(key=lambda x: x['start'])

        # Strategy: First try to pair LEF2/ORF603 locus with ORF1629 locus
        # This creates the main baculovirus recombination cassette
        lef2_regions = [r for r in homology_regions if r['locus'] == 'lef2_orf603']
        orf1629_regions = [r for r in homology_regions if r['locus'] == 'orf1629']

        used_indices = set()

        # Pair LEF2/ORF603 with ORF1629 (the main recombination boundaries)
        if lef2_regions and orf1629_regions:
            # Use the first (leftmost) LEF2/ORF603 region and first ORF1629 region
            lef2_region = lef2_regions[0]
            orf1629_region = orf1629_regions[0]

            # Ensure ORF1629 is downstream of LEF2/ORF603
            if orf1629_region['start'] > lef2_region['end']:
                distance = orf1629_region['start'] - lef2_region['end']
                if distance <= 20000:  # Allow larger distance for full cassette
                    module = {
                        'module_type': 'baculovirus_recombination_cassette',
                        'start': lef2_region['start'],
                        'end': orf1629_region['end'],
                        'strand': 1,
                        'name': 'Baculovirus Recombination Cassette (LEF2/ORF603 - ORF1629)',
                        'features': [lef2_region['index'], orf1629_region['index']],
                        'weight': 0.92,
                        'detection_method': 'rule_based',
                        'notes': f"Bounded by {lef2_region['name']} and {orf1629_region['name']}",
                    }
                    modules.append(module)
                    # Mark all LEF2/ORF603 and ORF1629 regions as used
                    for r in lef2_regions + orf1629_regions:
                        used_indices.add(r['index'])

        # Detect expression cassettes: promoter (p10/polyhedrin) -> 3' UTR
        # These are nested within the recombination cassette
        promoter_regions = []
        utr_regions = []

        for r in homology_regions:
            name_lower = r['name'].lower()
            # Promoters
            if 'p10 promoter' in name_lower or 'polyhedrin promoter' in name_lower or 'polh promoter' in name_lower:
                promoter_regions.append(r)
            # 3' UTRs
            elif "3' utr" in name_lower or "3'utr" in name_lower:
                utr_regions.append(r)

        # Also check for polyhedrin/p10 regions that might be promoters or UTRs
        for r in homology_regions:
            name_lower = r['name'].lower()
            if r['locus'] == 'expression' and r not in promoter_regions and r not in utr_regions:
                # If it's before the cargo region (lower position), likely a promoter
                # If it's after, likely a UTR - simple heuristic based on position
                if promoter_regions:
                    if r['start'] > promoter_regions[0]['end']:
                        utr_regions.append(r)
                elif utr_regions:
                    if r['end'] < utr_regions[0]['start']:
                        promoter_regions.append(r)

        # Pair promoters with UTRs to create expression cassettes
        if promoter_regions and utr_regions:
            promoter_regions.sort(key=lambda x: x['start'])
            utr_regions.sort(key=lambda x: x['start'])

            used_promoters = set()
            used_utrs = set()

            for promoter in promoter_regions:
                if promoter['index'] in used_promoters:
                    continue

                for utr in utr_regions:
                    if utr['index'] in used_utrs:
                        continue

                    # UTR should be downstream of promoter
                    if utr['start'] > promoter['end']:
                        distance = utr['start'] - promoter['end']
                        if 50 <= distance <= 10000:
                            # Determine promoter type
                            prom_name = promoter['name'].lower()
                            if 'p10' in prom_name:
                                cassette_name = 'p10 Expression Cassette'
                            elif 'polyhedrin' in prom_name or 'polh' in prom_name:
                                cassette_name = 'Polyhedrin Expression Cassette'
                            else:
                                cassette_name = 'Baculovirus Expression Cassette'

                            module = {
                                'module_type': 'baculovirus_expression_cassette',
                                'start': promoter['start'],
                                'end': utr['end'],
                                'strand': 1,
                                'name': cassette_name,
                                'features': [promoter['index'], utr['index']],
                                'weight': 0.90,
                                'detection_method': 'rule_based',
                                'notes': f"From {promoter['name']} to {utr['name']}",
                            }
                            modules.append(module)
                            used_promoters.add(promoter['index'])
                            used_utrs.add(utr['index'])
                            break

        return modules


    def _scan_for_frt_motifs(self, sequence: str) -> List[Dict[str, Any]]:
        """
        Scan the sequence for FRT motif sites that may not be annotated.

        This is useful when pLannotate doesn't detect FRT sites due to
        overlap with other features like 2μ ori.

        Returns list of detected FRT sites with positions.
        """
        import re
        frt_sites = []

        sequence_upper = sequence.upper()

        for frt_name, frt_seq in FRT_MOTIFS.items():
            # Search forward strand
            for match in re.finditer(frt_seq, sequence_upper):
                frt_sites.append({
                    'name': frt_name,
                    'start': match.start(),
                    'end': match.end(),
                    'strand': 1,
                    'sequence': frt_seq,
                })

            # Search reverse complement
            rev_comp = self._reverse_complement(frt_seq)
            for match in re.finditer(rev_comp, sequence_upper):
                frt_sites.append({
                    'name': frt_name,
                    'start': match.start(),
                    'end': match.end(),
                    'strand': -1,
                    'sequence': rev_comp,
                })

        # Remove duplicates (longer motifs that contain shorter ones)
        # Keep the longer match if positions overlap
        frt_sites.sort(key=lambda x: (x['start'], -len(x['sequence'])))
        filtered_sites = []
        for site in frt_sites:
            overlaps = False
            for existing in filtered_sites:
                # Check for overlap
                if not (site['end'] <= existing['start'] or site['start'] >= existing['end']):
                    overlaps = True
                    break
            if not overlaps:
                filtered_sites.append(site)

        return filtered_sites

    def _reverse_complement(self, seq: str) -> str:
        """Return reverse complement of DNA sequence."""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                      'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
        return ''.join(complement.get(base, base) for base in reversed(seq))

    # Constants for the lentiviral-regulatory detector
    _LENTI_UPSTREAM_FTYPES = {'promoter', 'enhancer', 'intron'}
    _LENTI_UPSTREAM_KB_CLASSES = {'promoter', 'enhancer', 'intron', '5_utr', '5utr'}
    _LENTI_DOWNSTREAM_FTYPES = {'cds', 'polya_signal', 'terminator'}
    _LENTI_DOWNSTREAM_KB_CLASSES = {'cds', '3_utr', '3utr', 'polya_signal', 'wpre', 'terminator'}
    _LENTI_BACTERIAL_SCOPES = {'bacterial', 'prokaryotic', 'e_coli', 'bacteria'}

    def _feature_host_scope(self, feat: Dict[str, Any]) -> List[str]:
        hs = feat.get('host_scope', [])
        if isinstance(hs, str):
            hs = [hs]
        return [(s or '').lower() for s in hs]

    def _is_lenti_upstream_feature(self, feat: Dict[str, Any]) -> bool:
        """Promoter/enhancer/intron that passes animal (non-bacterial) gating."""
        ftype = (feat.get('type') or '').lower()
        kb_class = (feat.get('kb_class') or '').lower()
        subclass = (feat.get('kb_subclass') or '').lower()
        name = (feat.get('name') or '').lower()
        has_tag = (
            ftype in self._LENTI_UPSTREAM_FTYPES
            or kb_class in self._LENTI_UPSTREAM_KB_CLASSES
            or subclass in self._LENTI_UPSTREAM_KB_CLASSES
            or 'promoter' in name
            or 'enhancer' in name
            or 'intron' in name
        )
        if not has_tag:
            return False
        host = self._feature_host_scope(feat)
        # Exclude features explicitly scoped to bacterial-only host.
        if host and all(h in self._LENTI_BACTERIAL_SCOPES for h in host):
            return False
        return True

    def _is_lenti_downstream_anchor(self, feat: Dict[str, Any]) -> bool:
        """CDS, WPRE, 3'UTR, polyA signal — defines where downstream regulatory begins."""
        ftype = (feat.get('type') or '').lower()
        kb_class = (feat.get('kb_class') or '').lower()
        subclass = (feat.get('kb_subclass') or '').lower()
        name = (feat.get('name') or '').lower()
        if ftype in self._LENTI_DOWNSTREAM_FTYPES:
            return True
        if kb_class in self._LENTI_DOWNSTREAM_KB_CLASSES:
            return True
        if subclass in self._LENTI_DOWNSTREAM_KB_CLASSES:
            return True
        if any(x in name for x in ('wpre', '3\'utr', '3utr', 'polya', 'poly(a)', 'terminator')):
            return True
        return False


    # Generation classification constants
    _HIV_ACCESSORY_GENES = {'vif', 'vpr', 'vpu', 'nef', 'tat', 'rev'}
    _SIN_3PRIME_MARKERS = ('delta-u3', 'delta_u3', '\u0394u3', 'sin', 'self-inactivating')
    _LENTI_CIS_ELEMENTS = {
        # submodule_subtype : (name-keywords to match)
        'psi':            ('psi', '\u03a8', 'hiv-1 psi', 'packaging signal'),
        'rre':            ('rre', 'rev response element'),
        'cppt_cts':       ('cppt', 'cts', 'central polypurine tract', 'central termination'),
        'wpre':           ('wpre', 'woodchuck'),
        'gag_fragment':   ('gag',),  # partial gag used as encapsidation signal
        '5_ltr':          (),  # populated from feat_5 directly
        '3_ltr':          (),  # populated from feat_3 directly
    }

    def _classify_lentiviral_generation(
        self,
        features: List[Dict[str, Any]],
        feat_5: Dict[str, Any],
        feat_3: Dict[str, Any],
        payload_start: int,
        payload_end: int,
    ) -> Dict[str, Any]:
        """Classify lentiviral payload as Generation 1 / 2 / 3.

        Returns dict with keys: generation, is_sin, accessory_genes_present,
        has_external_promoter, rationale.

        Heuristic:
          - 3' LTR marked Delta-U3 / SIN → payload is SIN (Gen 3 or Gen 2).
          - HIV accessory genes (vif/vpr/vpu/nef) present inside the payload
            → Gen 1 or Gen 2.
          - External non-LTR Pol II promoter (CMV / EF-1a / RSV / SV40 /
            CAG / PGK) driving the payload → Gen 2 or Gen 3 (Gen 1 uses
            the LTR U3 as its own promoter).
          Composite decision tree:
             is_SIN + external_promoter + no_accessories => "Gen 3"
             is_SIN + external_promoter                  => "Gen 3"
             is_SIN                                      => "Gen 3"
             not is_SIN + external_promoter + few_accessories => "Gen 2"
             accessory_genes_present                     => "Gen 1"
             otherwise                                   => "Gen 2"
        """
        name_3 = (feat_3.get('name') or '').lower()
        is_sin = any(marker in name_3 for marker in self._SIN_3PRIME_MARKERS)

        accessory_hits = []
        for f in features:
            nm = (f.get('name') or '').lower()
            if f.get('start', -1) < payload_start or f.get('end', 10**9) > payload_end:
                continue
            tokens = [t.strip('-()') for t in nm.replace('/', ' ').split()]
            for gene in self._HIV_ACCESSORY_GENES:
                if gene in tokens:
                    accessory_hits.append(gene)
                    break

        accessory_genes_present = sorted(set(accessory_hits))

        _EXTERNAL_POL2_PROMOTER_NAMES = (
            'cmv', 'ef-1', 'ef1', 'ef-1-alpha', 'htlv', 'pgk', 'rsv',
            'sv40 promoter', 'cag', 'ubc', 'sffv', 'mscv',
        )
        has_external_promoter = False
        for f in features:
            nm = (f.get('name') or '').lower()
            ftype = (f.get('type') or '').lower()
            if f.get('end', -1) > feat_5.get('start', 0):
                continue  # must sit before the 5' LTR to count as external
            if ftype not in ('promoter', 'enhancer'):
                continue
            if any(k in nm for k in _EXTERNAL_POL2_PROMOTER_NAMES):
                has_external_promoter = True
                break

        # Classification
        if is_sin and has_external_promoter:
            generation = 'Gen 3'
            rationale = (
                "self-inactivating 3' LTR (ΔU3) + external Pol II promoter"
                + (f" + no accessory genes detected"
                   if not accessory_genes_present else
                   f" (residual accessory genes: {', '.join(accessory_genes_present)})")
                + "."
            )
        elif is_sin:
            generation = 'Gen 3'
            rationale = "self-inactivating 3' LTR without external-promoter evidence."
        elif accessory_genes_present and len(accessory_genes_present) >= 2:
            generation = 'Gen 1'
            rationale = (
                f"intact 3' LTR + HIV accessory genes present "
                f"({', '.join(accessory_genes_present)})."
            )
        elif has_external_promoter:
            generation = 'Gen 2'
            rationale = "intact 3' LTR + external Pol II promoter (packaging-separated)."
        else:
            generation = 'Gen 2'
            rationale = "intact 3' LTR, no SIN marker, no strong generation signal."

        return {
            'generation': generation,
            'is_sin': is_sin,
            'accessory_genes_present': accessory_genes_present,
            'has_external_promoter': has_external_promoter,
            'rationale': rationale,
        }

    def _detect_lentiviral_cis_elements(
        self,
        features: List[Dict[str, Any]],
        feat_5: Dict[str, Any],
        feat_3: Dict[str, Any],
        payload_start: int,
        payload_end: int,
    ) -> List[Dict[str, Any]]:
        """Emit lentiviral cis-regulatory elements (Ψ, RRE, cPPT/CTS, WPRE,
        partial gag, 5'/3' LTR) as peer submodules of the payload."""
        subs: List[Dict[str, Any]] = []

        # 5' and 3' LTR boundary markers
        subs.append({
            'module_type': 'lentiviral_cis_element',
            'subtype': '5_ltr',
            'start': feat_5['start'],
            'end': feat_5['end'],
            'strand': feat_5.get('strand', 1),
            'name': feat_5.get('name', "5' LTR"),
            'weight': 0.98,
            'detection_method': 'rule_based',
            'rule_id': 'LENTI-CIS-LTR-5',
            'notes': 'Lentiviral 5\' LTR — defines upstream integration boundary.',
        })
        subs.append({
            'module_type': 'lentiviral_cis_element',
            'subtype': '3_ltr',
            'start': feat_3['start'],
            'end': feat_3['end'],
            'strand': feat_3.get('strand', 1),
            'name': feat_3.get('name', "3' LTR"),
            'weight': 0.98,
            'detection_method': 'rule_based',
            'rule_id': 'LENTI-CIS-LTR-3',
            'notes': 'Lentiviral 3\' LTR — defines downstream integration boundary.',
        })

        # Scan for named cis elements inside the payload
        already = set()  # (subtype, start, end) dedupe
        for f in features:
            fs = f.get('start', -1)
            fe = f.get('end', 10**9)
            if fs < payload_start or fe > payload_end:
                continue
            nm = (f.get('name') or '').lower()
            for subtype, keywords in self._LENTI_CIS_ELEMENTS.items():
                if not keywords:
                    continue
                if any(k in nm for k in keywords):
                    key = (subtype, fs, fe)
                    if key in already:
                        continue
                    already.add(key)
                    subs.append({
                        'module_type': 'lentiviral_cis_element',
                        'subtype': subtype,
                        'start': fs,
                        'end': fe,
                        'strand': f.get('strand', 1),
                        'name': f.get('name', subtype),
                        'weight': 0.9,
                        'detection_method': 'rule_based',
                        'rule_id': f'LENTI-CIS-{subtype.upper()}',
                        'notes': f'Lentiviral cis element ({subtype}) inside payload.',
                    })
                    break
        return subs

    def _detect_ltr_payload(self, features: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Detect lentiviral payload bounded by 5' LTR and 3' LTR, and emit two
        peer regulatory modules (upstream + downstream) that interact with it.

        Emits three module types:
        - lentiviral_upstream_regulatory:
            starts at the first animal-scope promoter/enhancer feature and
            ends at the last promoter/enhancer/intron-class feature strictly
            before the 3' LTR. Can extend past the 5' LTR start (introns
            inside the LTR boundary count). Bacterial-only features excluded.
        - lentiviral_payload:
            existing LTR-to-LTR span.
        - lentiviral_downstream_regulatory:
            starts at the first CDS / WPRE / 3' UTR / polyA feature AFTER the
            3' LTR start and ends at the end of the polyA signal (or, if no
            polyA, at the end of the farthest downstream anchor).
        """
        modules = []

        # Collect any feature that looks like an LTR (type, KB class, or name)
        ltr_candidates = []
        for i, feat in enumerate(features):
            ftype = (feat.get('type') or '').upper()
            kb_class = (feat.get('kb_class') or '').upper()
            name = (feat.get('name') or '').lower()
            is_ltr = (
                ftype == 'LTR'
                or kb_class == 'LTR'
                or 'ltr' in name
                or "\u0394u3" in name
                or 'delta-u3' in name
                or 'delta_u3' in name
                or 'self-inactivating' in name
            )
            if is_ltr:
                ltr_candidates.append((i, feat))

        if len(ltr_candidates) < 2:
            return modules

        # Sort by start; the upstream (5') LTR is the leftmost, the downstream (3')
        # LTR is the rightmost. If there are more than two, we still pair first↔last
        # since that defines the integration-competent payload window.
        ltr_candidates.sort(key=lambda t: t[1].get('start', 0))
        idx_5, feat_5 = ltr_candidates[0]
        idx_3, feat_3 = ltr_candidates[-1]

        # Make the loop structure the same as before so the body below doesn't
        # need to change — wrap in a single-pair list.
        # Single canonical pair: leftmost ↔ rightmost LTR candidate
        if feat_3['start'] > feat_5['end']:
            payload_start = feat_5['start']
            payload_end = feat_3['end']

            # ── Lentiviral Payload + generation classification + cis-element submodules ──
            gen_info = self._classify_lentiviral_generation(
                features, feat_5, feat_3, payload_start, payload_end
            )
            cis_subs = self._detect_lentiviral_cis_elements(
                features, feat_5, feat_3, payload_start, payload_end
            )
            modules.append({
                'module_type': 'lentiviral_payload',
                'start': payload_start,
                'end': payload_end,
                'name': f"Lentiviral Payload ({gen_info['generation']})",
                'features': [idx_5, idx_3],
                'weight': 0.98,
                'detection_method': 'rule_based',
                'rule_id': 'LENTI-PAY-01',
                'generation': gen_info['generation'],
                'submodules': cis_subs,
                'metadata': {
                    'generation': gen_info['generation'],
                    'is_SIN': gen_info['is_sin'],
                    'accessory_genes_present': gen_info['accessory_genes_present'],
                    'has_external_promoter': gen_info['has_external_promoter'],
                },
                'notes': f"Bounded by {feat_5.get('name', '5LTR')} and {feat_3.get('name', '3LTR')}; {gen_info['rationale']}",
            })
            # Also emit each cis element as its own peer annotation
            for cis in cis_subs:
                modules.append(cis)

            # ── Lentiviral Upstream Regulatory ────────────────────────
            # Starts at the first animal-scope promoter/enhancer before
            # the 3' LTR; ends at the last upstream-class feature before
            # the 3' LTR (can sit inside the 5' LTR boundary).
            # Only consider upstream features that END at or before the 5' LTR
            # start — the upstream regulatory region must terminate at the LTR.
            upstream_matches = [
                (i, f) for i, f in enumerate(features)
                if f.get('end', 0) <= feat_5['start']
                and self._is_lenti_upstream_feature(f)
            ]
            if upstream_matches:
                ur_start = min(f['start'] for _, f in upstream_matches)
                ur_end = min(max(f['end'] for _, f in upstream_matches), feat_5['start'])
                ur_idxs = sorted({i for i, _ in upstream_matches})
                modules.append({
                    'module_type': 'lentiviral_upstream_regulatory',
                    'start': ur_start,
                    'end': ur_end,
                    'strand': feat_5.get('strand', 1),
                    'name': 'Lentiviral Upstream Regulatory',
                    'features': ur_idxs,
                    'weight': 0.96,
                    'detection_method': 'rule_based',
                    'rule_id': 'LENTI-UR-01',
                    'notes': (
                        f"{len(upstream_matches)} animal-scope promoter/enhancer/intron "
                        f"features before the 5' LTR ({feat_5['start']})."
                    ),
                })

            # ── Lentiviral Downstream Regulatory ──────────────────────
            # A polyA signal that sits *outside* (just past) the 3' LTR
            # and acts as the transcription terminator for the integrated
            # lentiviral payload. Constraints:
            #   - start  > feat_3.end   (outside the LTR, not inside it)
            #   - start <= feat_3.end + DR_MAX_GAP_BP  (close enough to
            #     count as "directly after" the payload)
            #   - feature is polyA-class (type/kb_class polyA_signal,
            #     or name contains "polya" / "poly(a)")
            # Exactly one DR is emitted per payload — the nearest polyA
            # satisfying the gap constraint.
            DR_MAX_GAP_BP = 100
            ltr_end = int(feat_3.get('end', 0))
            dr_window_end = ltr_end + DR_MAX_GAP_BP

            polya_candidates = []
            for i, f in enumerate(features):
                fs = int(f.get('start', 0))
                if fs <= ltr_end or fs > dr_window_end:
                    continue
                nm = (f.get('name') or '').lower()
                ft = (f.get('type') or '').lower()
                kbc = (f.get('kb_class') or '').lower()
                is_polya = (
                    ft == 'polya_signal'
                    or kbc == 'polya_signal'
                    or 'polya' in nm
                    or 'poly(a)' in nm
                )
                if is_polya:
                    polya_candidates.append((i, f))

            if polya_candidates:
                # Nearest polyA wins
                polya_candidates.sort(key=lambda t: t[1].get('start', 0))
                dr_idx, dr_feat = polya_candidates[0]
                modules.append({
                    'module_type': 'lentiviral_downstream_regulatory',
                    'start': int(dr_feat['start']),
                    'end': int(dr_feat['end']),
                    'strand': dr_feat.get('strand', feat_3.get('strand', 1)),
                    'name': 'Lentiviral Downstream Regulatory',
                    'features': [dr_idx],
                    'weight': 0.95,
                    'detection_method': 'rule_based',
                    'rule_id': 'LENTI-DR-01',
                    'metadata': {
                        'polya_feature': dr_feat.get('name', ''),
                        'gap_from_ltr_end_bp': int(dr_feat['start']) - ltr_end,
                    },
                    'notes': (
                        f"{dr_feat.get('name','polyA')} at "
                        f"{dr_feat.get('start')}-{dr_feat.get('end')} sits "
                        f"{int(dr_feat['start']) - ltr_end} bp past the "
                        f"3' LTR end; serves as the transcription terminator "
                        f"for the lentiviral payload."
                    ),
                })

        return modules

    def _detect_itr_payload(self, features: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Detect AAV payload bounded by ITRs.

        Looks for:
        - AAV2 ITR or AAV ITR (left/right)
        - Creates aav_payload spanning both ITRs
        """
        modules = []

        # Find ITR features
        itr_features = []

        for i, feat in enumerate(features):
            name = feat.get('name', '').lower()

            if any(pattern in name for pattern in [
                'aav', 'itr', 'inverted terminal repeat'
            ]):
                itr_features.append((i, feat))

        # Match pairs of ITRs
        if len(itr_features) >= 2:
            # Sort by position
            itr_features.sort(key=lambda x: x[1]['start'])

            # Take first and last ITR
            idx_5, feat_5 = itr_features[0]
            idx_3, feat_3 = itr_features[-1]

            if feat_3['start'] > feat_5['end']:
                module = {
                    'module_type': 'aav_payload',
                    'start': feat_5['start'],
                    'end': feat_3['end'],
                    'name': 'AAV Payload',
                    'features': [idx_5, idx_3],
                    'weight': 0.98,
                    'detection_method': 'rule_based',
                    'notes': f"Bounded by ITRs"
                }
                modules.append(module)

        return modules

    def _detect_tdna_module(self, features: List[Dict[str, Any]], sequence: str = None) -> List[Dict[str, Any]]:
        """
        Detect T-DNA module bounded by LB and RB T-DNA repeats.

        T-DNA (Transfer DNA) is the region of Agrobacterium Ti plasmid or binary vectors
        that is transferred into plant cells. The T-DNA is defined by:
        - LB (Left Border): 25bp imperfect repeat, defines the end of transfer
        - RB (Right Border): 25bp imperfect repeat, defines the start of transfer

        T-DNA transfer occurs directionally from RB to LB (5' to 3'), but the module
        is defined as the region FROM LB TO RB in the direction of the plasmid.

        For circular plasmids, handles both orientations:
        - Standard: LB position < RB position (T-DNA is between them)
        - Wrapped: RB position < LB position (T-DNA wraps around origin)
        """
        modules = []

        # Helper to get feature name - pLannotate uses 'Feature', rule-based uses 'name'
        def get_feature_name(feat):
            return str(feat.get('name', feat.get('Feature', '')))

        # Helper to get feature positions - pLannotate uses 'qstart'/'qend'
        def get_start(feat):
            val = feat.get('start', feat.get('qstart', 0))
            try:
                return int(val)
            except (ValueError, TypeError):
                return 0

        def get_end(feat):
            val = feat.get('end', feat.get('qend', 0))
            try:
                return int(val)
            except (ValueError, TypeError):
                return 0

        # Debug: Print all feature names looking for T-DNA
        print(f'[DEBUG T-DNA] Scanning {len(features)} features for T-DNA borders...')
        if features:
            print(f'[DEBUG T-DNA] Sample feature keys: {list(features[0].keys())}')

        for i, feat in enumerate(features):
            fname = get_feature_name(feat).lower()
            if 't-dna' in fname or 'border' in fname or 'tdna' in fname:
                print(f'[DEBUG T-DNA] POTENTIAL MATCH: "{get_feature_name(feat)}" at {get_start(feat)}-{get_end(feat)}')

        # Find LB and RB T-DNA repeat features
        lb_features = []
        rb_features = []

        for i, feat in enumerate(features):
            raw_name = get_feature_name(feat)
            # Normalize: replace underscores with spaces, strip parenthetical suffixes like (3)
            name = raw_name.lower().replace('_', ' ')
            # Also remove trailing parenthetical like "(3)" or "(1)"
            name = re.sub(r'\s*\(\d+\)\s*$', '', name)

            start = get_start(feat)
            end = get_end(feat)

            if start == 0 and end == 0:
                continue

            # LB T-DNA repeat patterns
            if any(pattern in name for pattern in [
                'lb t-dna', 'left border', 'lb border', 'lbt-dna',
                't-dna lb', 't-dna left', 'tdna lb', 'tdna left',
                'lb t dna'  # Handle space-separated version
            ]) or (('lb' in name or 'left' in name) and ('t-dna' in name or 't dna' in name)):
                lb_features.append({
                    'index': i,
                    'feat': feat,
                    'start': start,
                    'end': end,
                    'strand': feat.get('strand', 1),
                    'name': raw_name or 'LB T-DNA repeat',
                })
                print(f'[DEBUG T-DNA] Found LB: "{raw_name}" at {start}-{end}')

            # RB T-DNA repeat patterns
            elif any(pattern in name for pattern in [
                'rb t-dna', 'right border', 'rb border', 'rbt-dna',
                't-dna rb', 't-dna right', 'tdna rb', 'tdna right',
                'rb t dna'  # Handle space-separated version
            ]) or (('rb' in name or 'right' in name) and ('t-dna' in name or 't dna' in name)):
                rb_features.append({
                    'index': i,
                    'feat': feat,
                    'start': start,
                    'end': end,
                    'strand': feat.get('strand', 1),
                    'name': get_feature_name(feat) or 'RB T-DNA repeat',
                })
                print(f'[DEBUG T-DNA] Found RB: "{get_feature_name(feat)}" at {start}-{end}')

        print(f'[DEBUG T-DNA] Found {len(lb_features)} LB and {len(rb_features)} RB features')

        if not lb_features or not rb_features:
            return modules

        # Get sequence length for circular handling
        seq_len = len(sequence) if sequence else None

        # Sort features by position
        lb_features.sort(key=lambda x: x['start'])
        rb_features.sort(key=lambda x: x['start'])

        # Match LB with RB to create T-DNA modules
        used_lb = set()
        used_rb = set()

        for lb in lb_features:
            if lb['index'] in used_lb:
                continue

            best_rb = None
            best_distance = float('inf')

            for rb in rb_features:
                if rb['index'] in used_rb:
                    continue

                # Calculate distance - T-DNA goes from LB to RB
                if rb['start'] > lb['end']:
                    distance = rb['start'] - lb['end']
                    orientation = 'standard'
                elif seq_len and lb['start'] > rb['end']:
                    distance = (seq_len - lb['start']) + rb['end']
                    orientation = 'wrapped'
                else:
                    continue

                # T-DNA modules are typically 1-25 kb, allow up to 30kb
                if distance < 50 or distance > 30000:
                    continue

                if distance < best_distance:
                    best_distance = distance
                    best_rb = (rb, orientation)

            if best_rb:
                rb, orientation = best_rb

                if orientation == 'standard':
                    module_start = lb['start']
                    module_end = rb['end']
                else:
                    module_start = lb['start']
                    module_end = rb['end']

                module = {
                    'module_type': 'tdna_module',
                    'start': module_start,
                    'end': module_end,
                    'strand': 1,
                    'name': 'T-DNA Module',
                    'features': [lb['index'], rb['index']],
                    'weight': 0.97,
                    'detection_method': 'rule_based',
                    'notes': f"Bounded by {lb['name']} and {rb['name']}",
                    'lb_site': lb['name'],
                    'rb_site': rb['name'],
                    'orientation': orientation,
                    'is_circular_wrap': orientation == 'wrapped',
                }
                modules.append(module)
                used_lb.add(lb['index'])
                used_rb.add(rb['index'])

                print(f'[DEBUG T-DNA] Created T-DNA module: {module_start}-{module_end} ({orientation})')

        return modules

    def _detect_flank_cargo_modules(self, features: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Detect flank:cargo pattern modules.

        Patterns include:
        - Cre/loxP: loxP sites flanking cargo
        - FLP/FRT: FRT sites flanking cargo
        - Gateway: attB, attL, attR, attP site pairs
        - Tn7: Tn7L/Tn7R
        - Sleeping Beauty: SB ITR L/R
        - PiggyBac: PB ITR 5'/3'
        """
        modules = []

        # Get high-weight flank:cargo rules
        flank_rules = [r for r in self.rules
                      if r.get('weight', 0) >= 0.9
                      and 'flanks:cargo' in r.get('location_constraint', '')]

        for rule in flank_rules:
            module_type = rule['module_type']
            required_features = rule.get('features', [])

            # Find matching feature pairs
            feature_5_candidates = []
            feature_3_candidates = []

            for i, feat in enumerate(features):
                name = feat.get('name', '').lower()

                # Check if feature matches any required pattern
                for req_feat in required_features:
                    req_lower = req_feat.lower()

                    if req_lower in name:
                        # Determine if 5' or 3' based on naming
                        if any(x in name for x in ['left', '_l)', '(l)', "5'", '_5_', "loxp"]):
                            feature_5_candidates.append((i, feat))
                        elif any(x in name for x in ['right', '_r)', '(r)', "3'", '_3_']):
                            feature_3_candidates.append((i, feat))
                        else:
                            # Generic match - could be either
                            feature_5_candidates.append((i, feat))
                            feature_3_candidates.append((i, feat))

            # Match pairs
            for idx_5, feat_5 in feature_5_candidates:
                for idx_3, feat_3 in feature_3_candidates:
                    if idx_5 == idx_3:
                        continue

                    if feat_3['start'] > feat_5['end']:
                        # Check distance constraint if specified
                        distance_kb = (feat_3['end'] - feat_5['start']) / 1000

                        # Apply distance limits based on module type
                        max_distance = 1000  # Default 1 Mb
                        if 'tn7' in module_type.lower():
                            max_distance = 15

                        if distance_kb <= max_distance:
                            module = {
                                'module_type': module_type,
                                'start': feat_5['start'],
                                'end': feat_3['end'],
                                'name': rule.get('notes', module_type.replace('_', ' ').title()),
                                'features': [idx_5, idx_3],
                                'weight': rule['weight'],
                                'detection_method': 'rule_based',
                                'notes': f"Flanked by {feat_5['name']} and {feat_3['name']}"
                            }
                            modules.append(module)

        return modules

    def _detect_standalone_modules(self, features: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Detect standalone modules with single non-CDS feature and weight ≥0.9.

        Examples:
        - Origins: pUC ori, ColE1 ori, R6K γ ori, 2μ ori, CEN/ARS
        - Phage origins: f1 ori, M13 ori
        - Conjugation: oriT
        - Insulators: cHS4
        """
        modules = []

        # Get high-weight standalone rules
        standalone_rules = [r for r in self.rules
                          if r.get('weight', 0) >= 0.9
                          and 'standalone_module' in r.get('location_constraint', '')]

        def _feat_matches_token(feat_name_lower: str, token: str) -> bool:
            """Token-boundary substring match: token must be a substring of feature name.

            Never does the reverse (feature-name-in-token) because that causes
            short generic names like "ori" to match long composite tokens like
            "oripebna1" and mislabel a pUC ori as EBV oriP.
            """
            tok = token.strip().lower()
            if not tok:
                return False
            return tok in feat_name_lower

        def _rule_matches_feature_set(required_features_list, rule_type_lower, feats):
            """Return list of feature index sets that satisfy the rule.

            - presence_all: every required token (split on " + ") must match some feature.
                            Alternatives within a token use " | ".
            - presence: any one required token match fires a hit (one per feature match).
            Each required-token string may contain " + " for conjunction and
            " | " for alternatives (same conventions as the heuristics CSV).
            """
            # Flatten: the CSV often gives us a single string like 'oriP + EBNA1'
            # or 'BK virus ori | BKV ori | polyoma ori'. Merge list entries with
            # " + " so a list ["oriP + EBNA1"] and a list ["oriP", "EBNA1"]
            # both expand to the same conjuncts.
            joined = ' + '.join(str(x) for x in required_features_list)
            # Conjuncts separated by " + "; each conjunct may be a disjunction on " | ".
            conjuncts = [c.strip() for c in joined.split(' + ') if c.strip()]

            hits_per_conjunct = []  # list of list of (feat_idx, matching_token)
            for conj in conjuncts:
                alts = [a.strip() for a in conj.split(' | ') if a.strip()]
                matching = []
                for i, f in enumerate(feats):
                    nm = (f.get('name') or '').lower()
                    ftype = (f.get('type') or '').lower()
                    if ftype == 'cds':
                        # CDS filter preserved from legacy behaviour for standalone ori/regulatory rules.
                        continue
                    for alt in alts:
                        if _feat_matches_token(nm, alt):
                            matching.append((i, alt))
                            break
                hits_per_conjunct.append(matching)

            if 'presence_all' in rule_type_lower:
                # Every conjunct must have at least one matching feature.
                if not all(hits_per_conjunct):
                    return []
                # Return the feature set as the union of first matches per conjunct.
                chosen = [hits[0] for hits in hits_per_conjunct]
                return [chosen]
            else:
                # 'presence': each feature match is its own hit.
                out = []
                for hits in hits_per_conjunct:
                    for (fi, tok) in hits:
                        out.append([(fi, tok)])
                return out

        for rule in standalone_rules:
            module_type = rule['module_type']
            required_features = rule.get('features', [])
            rule_type = (rule.get('rule_type', '') or '').lower()

            if 'presence' not in rule_type:
                continue

            hit_sets = _rule_matches_feature_set(required_features, rule_type, features)
            for hit_set in hit_sets:
                if not hit_set:
                    continue
                idxs = [fi for fi, _ in hit_set]
                feats_hit = [features[i] for i in idxs]
                # Module span: from leftmost feature start to rightmost feature end.
                start = min(f['start'] for f in feats_hit)
                end = max(f['end'] for f in feats_hit)
                primary = feats_hit[0]
                module = {
                    'module_type': module_type,
                    'start': start,
                    'end': end,
                    'strand': primary.get('strand', 1),
                    'name': primary.get('name', module_type),
                    'features': idxs,
                    'weight': rule['weight'],
                    'detection_method': 'rule_based',
                    'notes': rule.get('notes', ''),
                }
                modules.append(module)

        return modules

    def filter_overlapping_modules(self, modules: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Filter overlapping modules, keeping higher weight modules.

        Exception: Nested baculovirus expression cassettes are allowed within
        baculovirus recombination cassettes.

        Args:
            modules: List of detected modules

        Returns:
            Filtered list with overlaps removed
        """
        if not modules:
            return []

        # Sort by weight (descending), then by size (larger first to establish parent modules first)
        sorted_modules = sorted(modules, key=lambda m: (-m.get('weight', 0), -(m.get('end', 0) - m.get('start', 0))))

        filtered = []
        for module in sorted_modules:
            # Check if overlaps with any already selected module
            overlaps = False
            for selected in filtered:
                if self._modules_overlap(module, selected):
                    # Exception: Allow nested baculovirus cassettes
                    if self._is_allowed_nested(module, selected):
                        continue  # Not considered an overlap
                    overlaps = True
                    break

            if not overlaps:
                filtered.append(module)
            else:
                print(f"[DEBUG FILTER] Filtered out: {module.get('name')} ({module.get('module_type')}) due to overlap", flush=True)

        return filtered

    def _is_allowed_nested(self, inner: Dict[str, Any], outer: Dict[str, Any]) -> bool:
        """
        Check if inner module is allowed to be nested within outer module.

        Allowed nesting:
        - Baculovirus expression cassettes inside baculovirus recombination cassettes
        - Any expression cassette inside a recombination cassette
        """
        inner_type = inner.get('module_type', '')
        outer_type = outer.get('module_type', '')

        # Check if inner is completely contained in outer
        inner_start = inner.get('start', 0)
        inner_end = inner.get('end', 0)
        outer_start = outer.get('start', 0)
        outer_end = outer.get('end', 0)

        is_contained = inner_start >= outer_start and inner_end <= outer_end

        if not is_contained:
            # Permit lentiviral peer modules (upstream_regulatory / payload / downstream_regulatory)
            # to partially overlap each other even without full containment —
            # upstream_regulatory is allowed to extend past the 5 LTR start and
            # thus partially overlap the payload.
            if ("lentiviral_" in inner_type and "lentiviral_" in outer_type
                    and inner_type != outer_type):
                return True
            return False

        # Allow expression cassettes inside recombination cassettes
        if 'expression' in inner_type and 'recombination' in outer_type:
            return True

        # Allow the three lentiviral peer modules to co-exist
        if "lentiviral_" in inner_type and "lentiviral_" in outer_type and inner_type != outer_type:
            return True

        # Allow nested baculovirus modules
        if 'baculovirus' in inner_type.lower() and 'baculovirus' in outer_type.lower():
            # Expression cassette inside recombination cassette
            if 'expression' in inner_type or 'Expression' in inner.get('name', ''):
                return True

        # Allow FRT sites inside FRT-flanked regions
        # FRT sites may be classified as yeast_replication when from 2u ori
        inner_name = inner.get("name", "").lower()
        outer_name = outer.get("name", "").lower()

        if outer_type == "frt_flanked_region" or "frt" in outer_name:
            # Allow yeast_replication (2u ori containing FRT) inside
            if inner_type == "yeast_replication":
                if "2u" in inner_name or "2micro" in inner_name:
                    return True
            # Allow actual FRT features
            if "frt" in inner_name:
                return True

        # Allow loxP sites inside loxP-flanked regions
        if outer_type == "loxp_flanked_region" or "lox" in outer_name:
            if "lox" in inner_name:
                return True

        return False

    def _modules_overlap(self, m1: Dict[str, Any], m2: Dict[str, Any]) -> bool:
        """Check if two modules overlap significantly (>50% of smaller module)."""
        start1, end1 = m1['start'], m1['end']
        start2, end2 = m2['start'], m2['end']

        # Calculate overlap
        overlap_start = max(start1, start2)
        overlap_end = min(end1, end2)

        if overlap_start >= overlap_end:
            return False

        overlap_len = overlap_end - overlap_start
        min_len = min(end1 - start1, end2 - start2)

        return overlap_len / min_len > 0.5
