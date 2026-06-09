"""
Mammalian Pol II Expression Cassette Detector with CDS Module Resolution (v2)
==============================================================================

Implements deterministic rule-based detection of mammalian Pol II expression cassettes
with proper CDS module boundary resolution using Kozak sequences and intron handling.

Structure:
    Upstream Regulatory Module -> CDS Module -> Downstream Regulatory Module

Key change in v2: Detect upstream and downstream regulatory modules FIRST,
then resolve CDS module boundaries within them.
"""

import re
from typing import List, Dict, Any, Optional, Tuple
from dataclasses import dataclass, field


# Kozak patterns (from heuristic_motif_detector.py)
KOZAK_STRONG = re.compile(r'GCC[AG]CCATGG', re.IGNORECASE)
KOZAK_ADEQUATE = re.compile(r'[AG]CCATGG', re.IGNORECASE)
KOZAK_WEAK = re.compile(r'[ACGT]{3}ATGG', re.IGNORECASE)

# Stop codons
STOP_CODONS = {'TAA', 'TAG', 'TGA'}


def _feat_strand(feat: dict) -> int:
    """Return the strand of a feature dict, consulting both `strand` and
    `direction` keys. `plannotate_annotations` features carry `direction`
    (1 or -1); modules built later carry `strand`. The detector accepts both
    shapes so the strand filter actually works for reverse-strand features."""
    if "strand" in feat:
        return feat["strand"]
    if "direction" in feat:
        return feat["direction"]
    return 1


def _feat_kb(feat: dict) -> dict:
    """Return the KB info dict for a feature, consulting both `kb_info` (legacy)
    and `kb_data` (plannotate_annotations) keys. Same fix as `_feat_strand`:
    without this the KB-based classifiers (feature_class lookups for polyA,
    UR, DR, host_scope filters) silently saw empty dicts and only the
    name-based fallbacks could match."""
    kb = feat.get("kb_info") or feat.get("kb_data") or {}
    return kb if isinstance(kb, dict) else {}


@dataclass
class CDSModule:
    """Resolved CDS module with proper boundaries."""
    start: int  # ATG position
    end: int    # Stop codon end or last exon end
    strand: int
    kozak_start: int
    kozak_strength: str  # 'strong', 'adequate', 'weak'
    exons: List[Tuple[int, int]] = field(default_factory=list)
    introns: List[Dict[str, Any]] = field(default_factory=list)
    stop_codon: Optional[str] = None
    stop_codon_pos: Optional[int] = None
    aa_length: int = 0
    initiated_by: str = 'kozak'  # 'kozak' or 'ires'

    def to_dict(self) -> Dict[str, Any]:
        return {
            'module_type': 'cds_module',
            'start': self.start,
            'end': self.end,
            'strand': self.strand,
            'kozak_start': self.kozak_start,
            'kozak_strength': self.kozak_strength,
            'exons': self.exons,
            'introns': self.introns,
            'stop_codon': self.stop_codon,
            'stop_codon_pos': self.stop_codon_pos,
            'aa_length': self.aa_length,
            'initiated_by': self.initiated_by,
            'detection_method': 'cds_module_resolution'
        }


@dataclass
class UpstreamRegulatoryModule:
    """Upstream regulatory module (promoter/enhancer to CDS start region)."""
    start: int
    end: int
    strand: int
    components: List[Dict[str, Any]] = field(default_factory=list)
    primary_promoter: Optional[Dict[str, Any]] = None

    def to_dict(self) -> Dict[str, Any]:
        return {
            'module_type': 'upstream_regulatory_module',
            'start': self.start,
            'end': self.end,
            'strand': self.strand,
            'components': [c.get('name', '') for c in self.components],
            'primary_promoter': self.primary_promoter.get('name') if self.primary_promoter else None,
            'detection_method': 'upstream_regulatory_detection'
        }


@dataclass
class DownstreamRegulatoryModule:
    """Downstream regulatory module (WPRE/3'UTR/polyA)."""
    start: int
    end: int
    strand: int
    components: List[Dict[str, Any]] = field(default_factory=list)
    polya: Optional[Dict[str, Any]] = None

    def to_dict(self) -> Dict[str, Any]:
        return {
            'module_type': 'downstream_regulatory_module',
            'start': self.start,
            'end': self.end,
            'strand': self.strand,
            'components': [c.get('name', '') for c in self.components],
            'polya': self.polya.get('name') if self.polya else None,
            'detection_method': 'downstream_regulatory_detection'
        }


@dataclass
class MammalianPol2Cassette:
    """Complete mammalian Pol II expression cassette."""
    start: int
    end: int
    strand: int
    upstream_regulatory: UpstreamRegulatoryModule
    cds_module: CDSModule
    downstream_regulatory: DownstreamRegulatoryModule
    weight: float = 0.95

    def to_dict(self) -> Dict[str, Any]:
        return {
            'module_type': 'mammalian_pol2_expression_cassette',
            'start': self.start,
            'end': self.end,
            'strand': self.strand,
            'upstream_regulatory': self.upstream_regulatory.to_dict(),
            'cds_module': self.cds_module.to_dict(),
            'downstream_regulatory': self.downstream_regulatory.to_dict(),
            'weight': self.weight,
            'detection_method': 'mammalian_pol2_detection'
        }


class MammalianPol2Detector:
    """
    Detect mammalian Pol II expression cassettes with proper CDS module resolution.

    v2 approach:
    1. Find promoter-polyA pairs (candidate cassettes)
    2. Detect upstream regulatory module boundaries
    3. Detect downstream regulatory module boundaries
    4. Resolve CDS module within the gap between upstream and downstream
    """

    # Feature classes for upstream regulatory elements
    # 'exon' included so non-coding 5' exons of CBh/CAG-style cassettes
    # (e.g. chicken β-actin Exon 1 in pX330) fold into the UR module
    # rather than dangling between UR and CDS.
    UPSTREAM_REGULATORY_CLASSES = {
        'promoter', 'enhancer', '5_utr', '5utr', 'utr_5', 'intron', 'exon',
        "5'utr", 'five_prime_utr'
    }

    # Bacterial / phage / synthetic-prokaryotic promoter names that must
    # NEVER be folded into a Pol II expression cassette's UR module — they
    # belong to selection cassettes, bacterial backbones, or T7/T3/SP6 IVT
    # systems and only happen to be on the same strand between a mammalian
    # promoter and the polyA that closes the Pol II cassette.
    BACTERIAL_PROMOTER_KEYWORDS = (
        "ampr promoter", "ampicillin promoter", "amp promoter",
        "kanr promoter", "kan promoter", "kanamycin promoter",
        "neor promoter", "neo promoter",
        "tetr promoter", "tet promoter",
        "lac promoter", "lacuv5", "lac uv5",
        "tac promoter", "trc promoter", "trp promoter", "tac-i", "trc-i",
        "pbad promoter", "arac promoter", "para promoter",
        "em7 promoter", "em7-promoter",
        "rsv promoter",  # used as bacterial in some shuttle backbones (kept narrow)
        "reca promoter",
        # phage / IVT promoters used as bacterial in plasmid context:
        "t7 promoter", "t3 promoter", "sp6 promoter",
        # rare but real:
        "rop promoter", "rrnb promoter",
    )

    # Feature classes for downstream regulatory elements
    DOWNSTREAM_REGULATORY_CLASSES = {
        'polya_signal', 'wpre', 'pre', '3_utr', '3utr', 'utr_3',
        "3'utr", 'three_prime_utr', 'terminator'
    }

    # Host scopes to exclude (bacterial)
    EXCLUDED_HOST_SCOPES = {'bacterial', 'prokaryotic', 'e_coli', 'bacteria'}

    # Animal host scopes to include
    ANIMAL_HOST_SCOPES = {
        'mammalian', 'insect', 'vertebrate', 'human', 'mouse', 'rat',
        'primate', 'rodent', 'animal', 'eukaryotic'
    }

    def __init__(self, features: List[Dict], sequence: str, orfs: List[Dict] = None, circular: bool = True):
        """
        Initialize detector.

        Args:
            features: List of feature dicts with kb_info
            sequence: Plasmid sequence string
            orfs: List of detected ORFs from orf_finder
            circular: Whether plasmid is circular
        """
        self.features = features
        self.sequence = sequence.upper()
        self.orfs = orfs or []
        self.circular = circular
        self.seq_len = len(sequence)

        # Debug: print feature names
        print(f"[DEBUG Pol2] Total features: {len(features)}")
        for f in features[:10]:
            print(f"  - {f.get('name', 'unnamed')}: {f.get('start')}-{f.get('end')}, kb_info={f.get('kb_info', {})}")

    def detect_mammalian_pol2_cassettes(self) -> Tuple[List['MammalianPol2Cassette'], List[Dict]]:
        """
        Detect all mammalian Pol II expression cassettes.

        v2 approach:
        1. Find animal Pol2 promoters
        2. Find polyA signals
        3. For each promoter-polyA pair:
           a. Detect upstream regulatory module (promoter -> last upstream element)
           b. Detect downstream regulatory module (first downstream element -> polyA end)
           c. Resolve CDS module in the gap between them

        Returns:
            Tuple of:
            - List of MammalianPol2Cassette objects
            - List of ORFs with overlapping CDS ORFs removed
        """
        cassettes = []

        # Step 1: Find animal Pol2 promoters
        promoters = self._find_animal_pol2_promoters()
        print(f"[DEBUG Pol2] Found {len(promoters)} animal Pol2 promoters")
        for p in promoters:
            print(f"  - {p.get('name')}: {p.get('start')}-{p.get('end')}")

        # Step 2: Find polyA signals
        polya_signals = self._find_polya_signals()
        print(f"[DEBUG Pol2] Found {len(polya_signals)} polyA signals")
        for p in polya_signals:
            print(f"  - {p.get('name')}: {p.get('start')}-{p.get('end')}")

        # Step 3: For each promoter, find compatible polyA and build cassette
        for promoter in promoters:
            strand = _feat_strand(promoter)

            # Find compatible polyA signals downstream
            compatible_polyas = self._find_compatible_polyas(promoter, polya_signals, strand)

            for polya in compatible_polyas:
                print(f"[DEBUG Pol2] Trying pair: {promoter.get('name')} -> {polya.get('name')}")

                # Step 3a: Detect upstream regulatory module
                upstream_module = self._detect_upstream_regulatory_module(promoter, polya, strand)
                if not upstream_module:
                    print(f"[DEBUG Pol2]   No upstream module found")
                    continue
                print(f"[DEBUG Pol2]   Upstream module: {upstream_module.start}-{upstream_module.end}")

                # Step 3b: Detect downstream regulatory module
                downstream_module = self._detect_downstream_regulatory_module_v2(promoter, polya, strand)
                if not downstream_module:
                    print(f"[DEBUG Pol2]   No downstream module found")
                    continue
                print(f"[DEBUG Pol2]   Downstream module: {downstream_module.start}-{downstream_module.end}")

                # Step 3c: Resolve CDS module in the gap
                cds_module = self._resolve_cds_module_in_gap(
                    upstream_module, downstream_module, strand
                )
                if not cds_module:
                    print(f"[DEBUG Pol2]   No CDS module resolved")
                    continue
                print(f"[DEBUG Pol2]   CDS module: {cds_module.start}-{cds_module.end} ({cds_module.aa_length} aa)")

                # Create complete cassette
                cassette = MammalianPol2Cassette(
                    start=upstream_module.start,
                    end=downstream_module.end,
                    strand=strand,
                    upstream_regulatory=upstream_module,
                    cds_module=cds_module,
                    downstream_regulatory=downstream_module,
                    weight=0.95
                )

                cassettes.append(cassette)
                print(f"[DEBUG Pol2]   Created cassette: {cassette.start}-{cassette.end}")

        # Deduplicate overlapping cassettes
        cassettes = self._deduplicate_cassettes(cassettes)
        print(f"[DEBUG Pol2] Final cassettes after dedup: {len(cassettes)}")

        # Get resolved CDS regions from cassettes
        resolved_cds_regions = [
            {'start': c.cds_module.start, 'end': c.cds_module.end, 'strand': c.strand}
            for c in cassettes
        ]

        # Remove overlapping ORFs
        filtered_orfs = self._remove_overlapping_orf_modules(resolved_cds_regions)

        return cassettes, filtered_orfs

    def _find_animal_pol2_promoters(self) -> List[Dict]:
        """
        Find promoters with animal host scope (not bacterial).

        Looks for:
        - Features with type 'promoter' or feature_class 'promoter'
        - Name contains 'promoter' or 'enhancer'
        - Excludes bacterial host scope
        """
        promoters = []

        for feat in self.features:
            kb = _feat_kb(feat)
            feat_type = (feat.get('type') or '').lower()
            feat_class = (kb.get('feature_class') or '').lower()
            name = (feat.get('name') or '').lower()

            # Check if it's a promoter or enhancer.
            # Beyond plain "promoter" / "enhancer", we recognize three
            # special Pol II regulatory-element shapes whose canonical
            # GenoLIB / SnapGene labels don't contain "promoter":
            #   - TRE / tetracycline response element / tetO  -> Tet-on/off
            #     (rtTA / tTA drives Pol II transcription from these)
            #   - UAS / Gal4-UAS                              -> Gal4 binary
            #   - tetO operator alone                          -> minP+tetO bundle
            is_promoter_enhancer = (
                feat_type in {'promoter', 'enhancer'} or
                feat_class in {'promoter', 'enhancer'} or
                'promoter' in name or 'enhancer' in name or
                'tre' in name.split() or  # whole-word "TRE"
                'tetracycline response element' in name or
                'teto' in name or
                ('uas' in name.split() and 'gal' in name)
            )

            if not is_promoter_enhancer:
                continue

            # Check host scope - exclude bacterial
            host_scope = kb.get('host_scope', [])
            if isinstance(host_scope, str):
                host_scope = [host_scope]
            host_scope = [h.lower() for h in host_scope if h]

            # Exclude if explicitly bacterial
            if any(h in self.EXCLUDED_HOST_SCOPES for h in host_scope):
                continue

            # Check polymerase class if available
            polymerase_class = (kb.get('polymerase_class') or '').lower()
            if polymerase_class and 'pol_iii' in polymerase_class:
                continue

            # Exclude U6, H1, 7SK (Pol III promoters by name)
            if any(x in name for x in ['u6', 'h1 ', '7sk']):
                continue

            # Include if mammalian/animal or no host scope specified
            if not host_scope or any(h in self.ANIMAL_HOST_SCOPES for h in host_scope):
                promoters.append(feat)

        return promoters

    def _find_polya_signals(self) -> List[Dict]:
        """Find polyA signal features.

        Lentiviral 3' LTRs (incl. SIN / Delta-U3 / dU3 variants) carry an
        endogenous polyA and act as the terminator for any Pol II cassette
        embedded between the LTRs (e.g. EF-1α -> Cas9 -> 3' LTR in
        lentiCRISPR v2). They get included here so the cassette can resolve.
        Separately-tagged backbone polyA signals (bGH / SV40 / hGH) match
        on name or kb_info.feature_class."""
        polya_signals = []

        for feat in self.features:
            kb = _feat_kb(feat)
            feat_class = (kb.get('feature_class') or '').lower()
            name = (feat.get('name') or '').lower()
            ftype = (feat.get('type') or '').lower()

            if feat_class == 'polya_signal':
                polya_signals.append(feat)
            elif re.search(r'poly[\s\-\(]?a|polya|bgh.*poly|sv40.*poly|hgh.*poly|hsv.*poly|tk.*poly', name):
                polya_signals.append(feat)
            # Lentiviral 3' LTR carries an endogenous polyA — accept as terminator.
            elif ftype == 'ltr' and re.search(r'3.{0,3}ltr|sin|delta.?u3|du3|Δu3', name):
                polya_signals.append(feat)

        return polya_signals

    def _find_compatible_polyas(self, promoter: Dict, polya_signals: List[Dict], strand: int) -> List[Dict]:
        """Find polyA signals that are downstream of the promoter and compatible."""
        compatible = []

        for polya in polya_signals:
            # Same strand check
            if _feat_strand(polya) != strand:
                continue

            # Order check (polyA downstream of promoter)
            if strand == 1:
                if polya['start'] <= promoter['end']:
                    continue
                distance = polya['start'] - promoter['end']
            else:
                if polya['end'] >= promoter['start']:
                    continue
                distance = promoter['start'] - polya['end']

            # Distance check (100bp to 20kb)
            if distance < 100 or distance > 20000:
                continue

            compatible.append(polya)

        # Sort by distance (closest first)
        if strand == 1:
            compatible.sort(key=lambda p: p['start'])
        else:
            compatible.sort(key=lambda p: -p['end'])

        return compatible

    def _detect_upstream_regulatory_module(
        self,
        promoter: Dict,
        polya: Dict,
        strand: int
    ) -> Optional[UpstreamRegulatoryModule]:
        """
        Detect upstream regulatory module.

        Starts with the promoter/enhancer feature.
        Ends with the last promoter, enhancer, 5'UTR, or intron feature
        before the CDS region begins.
        """
        # Detect a wrap-promoter (start > end means the promoter spans the
        # origin of the circular plasmid — e.g. chicken β-actin in pX330 at
        # 8412..186 on an 8484 bp circle).
        promoter_wraps = promoter['start'] > promoter['end']

        if strand == 1:
            region_start = promoter['start']
            region_end = polya['start']
        else:
            region_start = polya['end']
            region_end = promoter['end']

        def _in_region(feat_start, feat_end):
            if strand == 1:
                if promoter_wraps:
                    # Two segments: [promoter.start, seq_len) ∪ [0, polya.start)
                    return (feat_start >= region_start) or (feat_start < region_end)
                return region_start <= feat_start < region_end
            else:
                if promoter_wraps:
                    return (feat_end <= region_end) or (feat_end > region_start)
                return region_start < feat_end <= region_end

        # Find all upstream regulatory features in the region
        regulatory_features = []

        for feat in self.features:
            if _feat_strand(feat) != strand:
                continue
            if not _in_region(feat['start'], feat['end']):
                continue

            kb = _feat_kb(feat)
            feat_class = (kb.get('feature_class') or '').lower()
            feat_type = (feat.get('type') or '').lower()
            name = (feat.get('name') or '').lower()

            # Check if it's an upstream regulatory element
            is_upstream_reg = (
                feat_class in self.UPSTREAM_REGULATORY_CLASSES or
                feat_type in {'promoter', 'enhancer', 'intron', 'exon'} or
                any(x in name for x in ['promoter', 'enhancer', "5'utr", '5utr',
                                        'intron', 'exon', 'kozak'])
            )

            # Exclude bacterial / phage / synthetic-prokaryotic promoters by
            # name. These belong to selection cassettes or IVT systems and
            # must not pad the bounds of a Pol II expression cassette's UR
            # module.
            is_bacterial_promoter = (
                feat_type == 'promoter'
                and any(k in name for k in self.BACTERIAL_PROMOTER_KEYWORDS)
            )
            kb_host = (kb.get('host_scope') or '').lower()
            if kb_host in self.EXCLUDED_HOST_SCOPES:
                is_bacterial_promoter = True

            if is_upstream_reg and not is_bacterial_promoter:
                regulatory_features.append(feat)
            elif is_bacterial_promoter:
                print(f"  [UR skip] excluded bacterial promoter '{feat.get('name','?')}' {feat['start']}-{feat['end']}")

        if not regulatory_features:
            return None

        # Find the extent of upstream regulatory region.
        # For wrap-promoter cassettes we preserve the wrap (start > end means
        # the module spans the origin). Pre-wrap features sit at [start, seq_len),
        # post-wrap features sit at [0, region_end). module_start = promoter.start
        # (anchors the wrap); module_end = max end across post-wrap features.
        if strand == 1 and promoter_wraps:
            module_start = promoter['start']
            post_wrap = [f for f in regulatory_features if f['end'] <= region_end]
            module_end = max((f['end'] for f in post_wrap), default=promoter['end'])
        elif strand == -1 and promoter_wraps:
            module_end = promoter['end']
            pre_wrap = [f for f in regulatory_features if f['start'] >= region_start]
            module_start = min((f['start'] for f in pre_wrap), default=promoter['start'])
        else:
            module_start = min(f['start'] for f in regulatory_features)
            module_end = max(f['end'] for f in regulatory_features)

        # Sort features by position.
        # For wrap-promoter cassettes, features in [promoter.start, seq_len) come
        # before features in [0, region_end), so sort accordingly.
        def _sort_key_fwd(f):
            if promoter_wraps and f['start'] < region_start:
                # Post-origin features: place after pre-origin ones
                return self.seq_len + f['start']
            return f['start']

        if strand == 1:
            regulatory_features.sort(key=_sort_key_fwd)
        else:
            regulatory_features.sort(key=lambda f: -f['end'])

        return UpstreamRegulatoryModule(
            start=module_start,
            end=module_end,
            strand=strand,
            components=regulatory_features,
            primary_promoter=promoter
        )

    def _detect_downstream_regulatory_module_v2(
        self,
        promoter: Dict,
        polya: Dict,
        strand: int
    ) -> Optional[DownstreamRegulatoryModule]:
        """
        Detect downstream regulatory module.

        Finds the first downstream element (WPRE, 3'UTR, or polyA) and extends to polyA end.
        """
        if strand == 1:
            region_start = promoter['end']
            region_end = polya['end']
        else:
            region_start = polya['start']
            region_end = promoter['start']

        # Find all downstream regulatory features
        downstream_features = []

        for feat in self.features:
            if _feat_strand(feat) != strand:
                continue

            feat_start = feat['start']
            feat_end = feat['end']

            # Check if feature is in the cassette region
            if strand == 1:
                if not (region_start < feat_start <= region_end):
                    continue
            else:
                if not (region_start <= feat_end < region_end):
                    continue

            kb = _feat_kb(feat)
            feat_class = (kb.get('feature_class') or '').lower()
            name = (feat.get('name') or '').lower()

            # Check if it's a downstream regulatory element
            is_downstream_reg = (
                feat_class in self.DOWNSTREAM_REGULATORY_CLASSES or
                any(x in name for x in ['wpre', 'pre', "3'utr", '3utr', 'poly', 'terminator'])
            )

            if is_downstream_reg:
                downstream_features.append(feat)

        # Always include the polyA signal
        if polya not in downstream_features:
            downstream_features.append(polya)

        if not downstream_features:
            return None

        # Find the extent of downstream regulatory region
        if strand == 1:
            module_start = min(f['start'] for f in downstream_features)
            module_end = max(f['end'] for f in downstream_features)
        else:
            module_start = min(f['start'] for f in downstream_features)
            module_end = max(f['end'] for f in downstream_features)

        # Sort by position
        downstream_features.sort(key=lambda f: f['start'] if strand == 1 else -f['end'])

        return DownstreamRegulatoryModule(
            start=module_start,
            end=module_end,
            strand=strand,
            components=downstream_features,
            polya=polya
        )

    def _resolve_cds_module_in_gap(
        self,
        upstream_module: UpstreamRegulatoryModule,
        downstream_module: DownstreamRegulatoryModule,
        strand: int
    ) -> Optional[CDSModule]:
        """
        Resolve CDS module in the gap between upstream and downstream regulatory modules.

        Steps:
        1. Define the CDS region as the gap between upstream.end and downstream.start
        2. Search for Kozak sequence at the start of this gap (3bp before to 15bp after)
        3. Find stop codon before downstream module starts
        4. Handle introns if present
        """
        if strand == 1:
            gap_start = upstream_module.end
            gap_end = downstream_module.start
        else:
            gap_start = downstream_module.end
            gap_end = upstream_module.start

        print(f"[DEBUG CDS] Gap region: {gap_start}-{gap_end}")

        if gap_end <= gap_start:
            print(f"[DEBUG CDS] No gap between upstream and downstream")
            return None

        # Search for Kozak at the start of the gap.
        # Window: UR end - 10 bp tolerance upstream → UR end + 150 bp downstream.
        # The 10 bp upstream tolerance catches ATGs whose Kozak motif barely
        # straddles the last codons of a 5' UTR / exon-classified UR feature
        # (more permissive than the prior 3 bp), but the window is still tight
        # enough that we never search deep into upstream regulatory features.
        # The 150 bp downstream window catches CBh/CAG-style cassettes
        # (pX330: ~110 bp between intron-1 end and 3xFLAG-Cas9 ATG) and still
        # rejects pairings where CDS sits >150 bp past the regulatory block
        # (likely a different cassette entirely).
        if strand == 1:
            search_start = max(0, gap_start - 10)
            search_end = min(self.seq_len, gap_start + 150)
        else:
            search_start = max(0, gap_end - 150)
            search_end = min(self.seq_len, gap_end + 10)

        print(f"[DEBUG CDS] Searching for Kozak in {search_start}-{search_end}")
        print(f"[DEBUG CDS] Sequence: {self.sequence[search_start:search_end]}")

        kozak_match = self._find_kozak_in_region(search_start, search_end, strand)

        # ORF-anchored fallback: if no Kozak motif matches in the window, look
        # at self.orfs (passed from Step 2's find_orfs) for any same-strand
        # ORF whose ATG falls in the search window. BLAST-detected CDS
        # features (Cas9, GFP, etc.) often land a few codons past the real
        # ATG, and that real ATG may not have a textbook Kozak context — but
        # the orf_finder still locates it. Adopt the ORF boundaries when
        # found.
        if not kozak_match:
            for orf in self.orfs:
                if _feat_strand(orf) != strand:
                    continue
                if strand == 1:
                    if search_start <= orf['start'] < search_end:
                        print(f"[DEBUG CDS] Kozak miss; falling back to ORF at {orf['start']}-{orf['end']}")
                        kozak_match = {
                            'kozak_start': orf['start'],
                            'atg_position': orf['start'],
                            'strength': 'orf_fallback',
                            'strength_rank': 0,
                            'sequence': '(no Kozak — anchored to orf_finder ATG)',
                            'initiated_by': 'orf_fallback',
                        }
                        break
                else:
                    if search_start < orf['end'] <= search_end:
                        print(f"[DEBUG CDS] Kozak miss; falling back to ORF at {orf['start']}-{orf['end']}")
                        kozak_match = {
                            'kozak_start': orf['end'],
                            'atg_position': orf['end'],
                            'strength': 'orf_fallback',
                            'strength_rank': 0,
                            'sequence': '(no Kozak — anchored to orf_finder ATG)',
                            'initiated_by': 'orf_fallback',
                        }
                        break

        if not kozak_match:
            print(f"[DEBUG CDS] No Kozak and no ORF fallback in window")
            return None

        print(f"[DEBUG CDS] Kozak/ORF found: {kozak_match}")

        # Resolve CDS end
        cds_start = kozak_match['atg_position']
        max_cds_end = gap_end if strand == 1 else gap_start

        cds_end, exons, introns, stop_codon, stop_pos = self._resolve_cds_end(
            cds_start, max_cds_end, strand
        )

        if cds_end is None:
            print(f"[DEBUG CDS] No stop codon found")
            return None

        # Strict downstream-of-UR guard. The Kozak hunt allows a 3 bp
        # tolerance upstream of UR.end for motifs that straddle the boundary;
        # past that the CDS is upstream of the regulatory block and the
        # pairing is invalid (would mean the CDS isn't actually driven by
        # this UR).
        if strand == 1:
            if cds_start + 3 < upstream_module.end:
                print(f"[DEBUG CDS] Rejected: CDS start {cds_start} is upstream of UR end {upstream_module.end}")
                return None
        else:
            if cds_end > upstream_module.start + 3:
                print(f"[DEBUG CDS] Rejected: CDS end {cds_end} is downstream of UR start {upstream_module.start}")
                return None

        # Calculate AA length
        coding_length = sum(e[1] - e[0] for e in exons)
        aa_length = coding_length // 3

        return CDSModule(
            start=cds_start,
            end=cds_end,
            strand=strand,
            kozak_start=kozak_match['kozak_start'],
            kozak_strength=kozak_match['strength'],
            exons=exons,
            introns=introns,
            stop_codon=stop_codon,
            stop_codon_pos=stop_pos,
            aa_length=aa_length,
            initiated_by=kozak_match.get('initiated_by', 'kozak')
        )

    def _find_kozak_in_region(self, start: int, end: int, strand: int) -> Optional[Dict]:
        """Find best Kozak sequence in a region."""
        if strand == -1:
            region_seq = self._reverse_complement(self.sequence[start:end])
        else:
            region_seq = self.sequence[start:end]

        # Check patterns in priority order
        for match in KOZAK_STRONG.finditer(region_seq):
            if strand == 1:
                return {
                    'kozak_start': start + match.start(),
                    'atg_position': start + match.start() + 6,
                    'strength': 'strong',
                    'strength_rank': 3,
                    'sequence': match.group()
                }
            else:
                return {
                    'kozak_start': end - match.end(),
                    'atg_position': end - match.end() + 6,
                    'strength': 'strong',
                    'strength_rank': 3,
                    'sequence': match.group()
                }

        for match in KOZAK_ADEQUATE.finditer(region_seq):
            if strand == 1:
                return {
                    'kozak_start': start + match.start(),
                    'atg_position': start + match.start() + 4,
                    'strength': 'adequate',
                    'strength_rank': 2,
                    'sequence': match.group()
                }
            else:
                return {
                    'kozak_start': end - match.end(),
                    'atg_position': end - match.end() + 4,
                    'strength': 'adequate',
                    'strength_rank': 2,
                    'sequence': match.group()
                }

        for match in KOZAK_WEAK.finditer(region_seq):
            if strand == 1:
                return {
                    'kozak_start': start + match.start(),
                    'atg_position': start + match.start() + 3,
                    'strength': 'weak',
                    'strength_rank': 1,
                    'sequence': match.group()
                }
            else:
                return {
                    'kozak_start': end - match.end(),
                    'atg_position': end - match.end() + 3,
                    'strength': 'weak',
                    'strength_rank': 1,
                    'sequence': match.group()
                }

        return None

    def _resolve_cds_end(
        self,
        cds_start: int,
        max_end: int,
        strand: int
    ) -> Tuple[Optional[int], List[Tuple[int, int]], List[Dict], Optional[str], Optional[int]]:
        """Resolve CDS end position with intron handling."""
        exons = []
        introns = []

        # Find introns in the CDS region
        intron_features = self._find_introns_in_region(cds_start, max_end, strand)

        if not intron_features:
            # No introns - simple case
            stop_pos, stop_codon = self._find_stop_codon(cds_start, max_end, strand)
            if stop_pos:
                exons.append((cds_start, stop_pos + 3))
                return stop_pos + 3, exons, introns, stop_codon, stop_pos
            return None, [], [], None, None

        # Handle introns
        current_pos = cds_start
        frame_offset = 0

        for intron_feat in sorted(intron_features, key=lambda f: f['start']):
            intron_start = intron_feat['start']
            intron_end = intron_feat['end']

            gt_pos = self._find_splice_site(intron_start - 3, intron_start + 3, 'GT')
            if gt_pos is None:
                continue

            ag_pos = self._find_splice_site(intron_end - 3, intron_end + 3, 'AG')
            if ag_pos is None:
                continue

            stop_pos, stop_codon = self._find_stop_codon(current_pos, gt_pos, strand, frame_offset)
            if stop_pos:
                exons.append((current_pos, stop_pos + 3))
                return stop_pos + 3, exons, introns, stop_codon, stop_pos

            exon_length = gt_pos - current_pos
            exons.append((current_pos, gt_pos))
            frame_offset = (frame_offset + exon_length) % 3

            introns.append({
                'start': gt_pos,
                'end': ag_pos + 2,
                'donor': 'GT',
                'acceptor': 'AG',
                'feature': intron_feat.get('name', 'intron')
            })

            current_pos = ag_pos + 2

        stop_pos, stop_codon = self._find_stop_codon(current_pos, max_end, strand, frame_offset)
        if stop_pos:
            exons.append((current_pos, stop_pos + 3))
            return stop_pos + 3, exons, introns, stop_codon, stop_pos

        return None, [], [], None, None

    def _find_introns_in_region(self, start: int, end: int, strand: int) -> List[Dict]:
        """Find intron features in a region."""
        introns = []

        for feat in self.features:
            kb = _feat_kb(feat)
            feat_class = (kb.get('feature_class') or '').lower()
            feat_type = (feat.get('type') or '').lower()
            name = (feat.get('name') or '').lower()

            is_intron = (
                feat_class == 'intron' or
                feat_type == 'intron' or
                'intron' in name
            )

            if is_intron and start < feat['start'] < end:
                if _feat_strand(feat) == strand:
                    introns.append(feat)

        return introns

    def _find_splice_site(self, start: int, end: int, motif: str) -> Optional[int]:
        """Find splice site motif (GT or AG) in region."""
        region = self.sequence[max(0, start):min(self.seq_len, end)]

        for i in range(len(region) - 1):
            if region[i:i+2].upper() == motif:
                return start + i

        return None

    def _find_stop_codon(
        self,
        start: int,
        end: int,
        strand: int,
        frame_offset: int = 0
    ) -> Tuple[Optional[int], Optional[str]]:
        """Find first in-frame stop codon."""
        if strand == 1:
            adjusted_start = start + ((3 - frame_offset) % 3) if frame_offset else start

            for pos in range(adjusted_start, end - 2, 3):
                codon = self.sequence[pos:pos+3].upper()
                if codon in STOP_CODONS:
                    return pos, codon
        else:
            region = self._reverse_complement(self.sequence[start:end])
            adjusted_start = frame_offset

            for i in range(adjusted_start, len(region) - 2, 3):
                codon = region[i:i+3]
                if codon in STOP_CODONS:
                    pos = end - i - 3
                    return pos, codon

        return None, None

    def _remove_overlapping_orf_modules(self, resolved_cds_regions: List[Dict]) -> List[Dict]:
        """Remove CDS ORF modules that overlap with resolved CDS modules."""
        filtered_orfs = []

        for orf in self.orfs:
            overlaps = False

            for cds in resolved_cds_regions:
                if _feat_strand(orf) != cds['strand']:
                    continue

                overlap_start = max(orf['start'], cds['start'])
                overlap_end = min(orf['end'], cds['end'])

                if overlap_start < overlap_end:
                    overlap_len = overlap_end - overlap_start
                    orf_len = orf['end'] - orf['start']

                    if orf_len > 0 and overlap_len / orf_len > 0.5:
                        overlaps = True
                        break

            if not overlaps:
                filtered_orfs.append(orf)

        return filtered_orfs

    def _deduplicate_cassettes(self, cassettes: List['MammalianPol2Cassette']) -> List['MammalianPol2Cassette']:
        """Remove overlapping cassettes, keeping highest weight."""
        if not cassettes:
            return []

        sorted_cassettes = sorted(cassettes, key=lambda c: -c.weight)

        kept = []
        for cassette in sorted_cassettes:
            overlaps = False
            for kept_cassette in kept:
                if cassette.strand != kept_cassette.strand:
                    continue

                overlap_start = max(cassette.start, kept_cassette.start)
                overlap_end = min(cassette.end, kept_cassette.end)

                if overlap_end > overlap_start:
                    overlap_len = overlap_end - overlap_start
                    cassette_len = cassette.end - cassette.start

                    if cassette_len > 0 and overlap_len / cassette_len > 0.5:
                        overlaps = True
                        break

            if not overlaps:
                kept.append(cassette)

        return kept

    def _reverse_complement(self, seq: str) -> str:
        """Get reverse complement of a sequence."""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                      'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
        return ''.join(complement.get(base, base) for base in reversed(seq))


# Convenience function
def detect_mammalian_pol2_cassettes(
    features: List[Dict],
    sequence: str,
    orfs: List[Dict] = None,
    circular: bool = True
) -> Tuple[List[Dict], List[Dict]]:
    """
    Detect mammalian Pol II expression cassettes (simple interface).

    Args:
        features: List of feature dicts with kb_info
        sequence: Plasmid sequence string
        orfs: List of detected ORFs from orf_finder
        circular: Whether plasmid is circular

    Returns:
        Tuple of:
        - List of cassette dicts
        - List of filtered ORFs (overlapping CDS ORFs removed)
    """
    detector = MammalianPol2Detector(features, sequence, orfs, circular)
    cassettes, filtered_orfs = detector.detect_mammalian_pol2_cassettes()
    return [c.to_dict() for c in cassettes], filtered_orfs
