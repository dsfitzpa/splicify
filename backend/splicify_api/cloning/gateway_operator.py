"""
GatewayOperator — Gateway Cloning operator for the Splicify cloning library.

Handles BP and LR recombination reactions with:
- att site detection and orthogonality validation
- Automatic reaction type inference (BP vs LR)
- PCR primer design to add attB sites when needed
- Recombination simulation with product/byproduct sequences
- ccdB selection analysis
- Complete cost/time/risk metrics
- SeqViz visualization annotations

Gateway Biology:
- BP Reaction: attB + attP → attL (product) + attR (byproduct)
  Creates Entry clone with attL sites flanking insert
- LR Reaction: attL + attR → attB (byproduct) + attP (product)
  Creates Expression clone with insert in destination vector
"""
from __future__ import annotations

from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass, field
import logging

logger = logging.getLogger(__name__)

from .build_plan import (
    BuildStep,
    GatewayBuildPlan,
    GatewayJunctionPlan,
    OperatorMetrics,
)


@dataclass
class FragmentAnnotation:
    """Track fragment sources in Gateway product for visualization."""
    name: str
    start: int
    end: int
    source_module: str
    source_type: str  # "donor_backbone", "insert", "att_site"
    direction: int = 1
from .gateway_sites import (
    GATEWAY_ATT_SITES,
    AttSiteMatch,
    scan_att_sites,
    scan_for_ccdb,
    validate_orthogonality,
    parse_att_site_type,
    extract_core_sequence,
    get_recombination_products,
)
from .junction import Junction, build_junctions
from .lab_profile import DEFAULT_LAB_PROFILE, LabProfile

# Try to import thermodynamic calculator for primer design
try:
    from ..gibson_primers import ThermodynamicCalculator
    _THERMO_AVAILABLE = True
except Exception:
    _THERMO_AVAILABLE = False


# Gateway operator parameters
_PRIMER_ANNEAL_TM_TARGET = 60.0
_PRIMER_ANNEAL_LEN_MIN = 18
_PRIMER_ANNEAL_LEN_MAX = 28


class GatewayOperator:
    """
    Gateway Cloning operator.
    
    Usage:
        op = GatewayOperator()
        plan = op.evaluate(modules, topology="circular")
        print(plan.summary)
        print(f"Cost: ${plan.metrics.total_cost_usd:.2f}")
    """
    
    def _detect_expression_module(self, module_sequence: str, module_name: str, module_att_sites: list = None) -> dict:
        """Detect expression context by analyzing module name, promoters, and terminators.

        Checks module name for expression system hints, then scans for promoters
        and terminators to determine bacterial vs eukaryotic expression context.
        """
        from .gateway_sites import scan_att_sites

        result = {
            "is_expression_module": False,
            "orientation": "unknown",
            "promoter": None,
            "features": [],
            "is_bacterial": None
        }

        # Get att sites if not provided
        if module_att_sites is None:
            module_att_sites = scan_att_sites(module_sequence, fuzzy_threshold=0)

        # Check module name for expression system hints
        module_name_lower = module_name.lower()

        # Mammalian/eukaryotic indicators in name
        mammalian_keywords = ["mammalian", "ptracer", "ef1a", "cmv", "psv40", "pcmv", "pcdna", "plenti", "paav"]
        # Bacterial indicators in name
        bacterial_keywords = ["pet", "pgex", "pmal", "ptrc", "pbad", "plac", "bacterial", "e.coli", "ecoli"]

        name_indicates_mammalian = any(kw in module_name_lower for kw in mammalian_keywords)
        name_indicates_bacterial = any(kw in bacterial_keywords for kw in bacterial_keywords)

        # If module name gives us a clear hint and there are NO att sites, it's likely the insert
        if not module_att_sites and (name_indicates_mammalian or name_indicates_bacterial):
            result["is_expression_module"] = True
            result["orientation"] = "forward"
            result["is_bacterial"] = name_indicates_bacterial
            result["promoter"] = "source_indicated"
            result["features"].append("source_file_hint")
            logger.info(f"[Gateway] Expression context from source: {module_name} ({'bacterial' if name_indicates_bacterial else 'mammalian'})")
            return result

        # Common promoter sequences with bacterial classification
        promoters = [
            ("t7", "TAATACGACTCACTATAGGG", True),
            ("t3", "AATTAACCCTCACTAAAGGG", True),
            ("sp6", "ATTTAGGTGACACTATAG", True),
            ("lac", "TGTGTGAAATTGTTATCCGCTCACAATTCCACACAACATACGAGCCGGAAGCATAAAGTGTAAAGCCTGGGGTGCCTAATGAGTGAGCTAACTCACATTAATTGCGTTGCGCTCACTGCCCGCTTTCCAGTCGGGAAACCTGTCGTGCCAGC", True),
            ("trc", "GAGCTGTTGACAATTAATCATCGGCTCGTATAATGT", True),
            ("ara", "TCCATAAAAAATTATTCTCAACACTAAACAT", True),
            ("cmv", "GTTGACATTGATTATTGACTAGTTATTAATAGTAATCAATTACGGGGTC", False),
            ("sv40", "GCCCTTAATATAACTTCGTAT", False),
            ("ef1a", "GGGGATTGGGAAGACAATAGCAGGCATGC", False),
            ("cag", "TCGTTACATAACTTACGGTAAATGGCCCGCCTGGCTGACCGCCCAACGACCCCCGCCCATTGACGTCAATAATGACGTATGTTCCCATAGTAACGCC", False),
        ]

        # Bacterial terminators
        terminators = [
            ("t7_terminator", "CTAGCATAACCCCTTGGGGCCTCTAAACGGGTCTTGAGGGGTTTTTTG", True),
            ("rrnb_t1", "TGCCTGGCGGCAGTAGCGCGGTGGTCCCACCTGACCCCATGCC", True),
            ("rrnb_t2", "CCCGCTTTCCCTTATTATGGGTT", True),
        ]

        # Poly-A signals (eukaryotic terminators)
        polya_signals = [
            ("sv40_polya", "AACTTGTTTATTGCAGCTTATAATGGTTACAAATAAAGCAATAGCATCACAAATTTCACAAATAAAGCATTTTTTTCACTGCATTCTAGTTGTGGTTTGTCCAAACTCATCAATGTATCTTATCATGTCTGGATC", False),
            ("bgh_polya", "CTGTGCCTTCTAGTTGCCAGCCATCTGTTGTTTGCCCCTCCCCCGTGCCTTCCTTGACCCTGGAAGGTGCCACTCCCACTGTCCTTTCCTAATAAAATGAGGAAATTGCATCGCATTGTCTGAGTAGGTGTCATTCTATTCTGGGGGGTGGGGTGGGGCAGGACAGCAAGGGGGAGGATTGGGAAGACAATAGCAGGCATGCTGGGGATGCGGTGGGCTCTATGG", False),
        ]

        seq_upper = module_sequence.upper()

        # Search for promoters first
        if module_att_sites:
            first_att_pos = min(site.start for site in module_att_sites)

            for prom_name, prom_seq, is_bacterial in promoters:
                pos = seq_upper.find(prom_seq.upper())
                if pos != -1:
                    # Found promoter - check if upstream of att sites
                    if pos < first_att_pos:
                        result["is_expression_module"] = True
                        result["promoter"] = prom_name
                        result["orientation"] = "forward"
                        result["features"].append(f"promoter_{prom_name}")
                        result["is_bacterial"] = is_bacterial
                        logger.info(f"[Gateway] Detected {prom_name} promoter in {module_name} ({'bacterial' if is_bacterial else 'eukaryotic'})")
                        return result

            # If no promoter found, look for terminators to infer expression system
            for term_name, term_seq, is_bacterial in terminators + polya_signals:
                pos = seq_upper.find(term_seq.upper())
                if pos != -1:
                    # Found terminator - this indicates expression context
                    if pos > first_att_pos:
                        result["is_expression_module"] = True
                        result["promoter"] = f"terminator_{term_name}"
                        result["orientation"] = "forward"
                        result["features"].append(term_name)
                        result["is_bacterial"] = is_bacterial
                        logger.info(f"[Gateway] Detected {term_name} in {module_name} ({'bacterial' if is_bacterial else 'eukaryotic'})")
                        return result

        return result

    def _orient_insert_for_expression(self, insert_sequence: str, expression_orientation: str) -> tuple:
        """
        Orient insert sequence to match expression module orientation.

        Args:
            insert_sequence: The insert DNA sequence
            expression_orientation: 'forward' or 'reverse' from expression module

        Returns:
            tuple of (oriented_sequence, was_reversed)
        """
        if expression_orientation == "reverse":
            # Reverse complement the insert
            from Bio.Seq import Seq
            insert_seq_obj = Seq(insert_sequence)
            oriented = str(insert_seq_obj.reverse_complement())
            return oriented, True
        else:
            # Keep forward orientation
            return insert_sequence, False


    def __init__(self, lab_profile: Optional[LabProfile] = None) -> None:
        self.lab = lab_profile or DEFAULT_LAB_PROFILE
        self._thermo = ThermodynamicCalculator() if _THERMO_AVAILABLE else None
    
    def evaluate(
        self,
        modules: List[dict],
        topology: str = "circular"
    ) -> GatewayBuildPlan:
        """
        Evaluate Gateway cloning for a module list.
        
        Args:
            modules: List of module dicts with sequence, canonical_id, role
            topology: "circular" (default) or "linear"
            
        Returns:
            Complete GatewayBuildPlan with all design details
        """
        plan = GatewayBuildPlan(assembly_topology=topology)
        
        # Validate minimum modules
        if len(modules) < 2:
            plan.feasible = False
            plan.infeasibility_reasons.append(
                f"Gateway requires at least 2 modules (got {len(modules)})"
            )
            return plan
        
        plan.fragment_count = len(modules)
        
        # Step 1: Scan all modules for att sites
        module_att_sites = []
        for i, mod in enumerate(modules):
            seq = mod.get("sequence", "")
            sites = scan_att_sites(seq, fuzzy_threshold=2)
            module_att_sites.append({
                "module_index": i,
                "module_name": _module_name(mod),
                "sequence": seq,
                "att_sites": sites,
                "module": mod
            })

        # 2026-05-06: orientation-aware substrate classification per input
        # module. Surfaces design-time warnings when an input is itself a
        # self-recombining substrate (excision or inversion) — running BP/LR
        # clonase on it without a partner vector would self-delete or self-
        # invert. Stored on plan.metadata for the chat handler / viz to pick
        # up; warnings are emitted to plan.warnings.
        from .gateway_sites import classify_gateway_substrate
        plan.substrate_classifications = []
        for entry in module_att_sites:
            cls = classify_gateway_substrate(entry["att_sites"])
            for c in cls:
                rec = {
                    "module_index": entry["module_index"],
                    "module_name": entry["module_name"],
                    "kind": c["kind"],
                    "pair_label": c["pair_label"],
                    "cargo_span": [c["cargo_start"], c["cargo_end"]],
                }
                plan.substrate_classifications.append(rec)
                if c["kind"] == "excision":
                    plan.warnings.append(
                        f"[orientation] Module '{entry['module_name']}' carries a "
                        f"same-strand {c['pair_label']} pair — BP/LR clonase will "
                        f"EXCISE the cargo between them as a circular byproduct. "
                        f"Confirm a partner vector is provided or this is intended."
                    )
                elif c["kind"] == "inversion":
                    plan.warnings.append(
                        f"[orientation] Module '{entry['module_name']}' carries an "
                        f"opposite-strand outward-pointing {c['pair_label']} pair — "
                        f"BP/LR clonase will INVERT the cargo between them in place. "
                        f"This may not be the intended reaction."
                    )
                # intermolecular: no warning — standard substrate

        # Step 2: Infer reaction type (BP vs LR)
        reaction_type = self._infer_reaction_type(module_att_sites)
        plan.reaction_type = reaction_type
        
        # Step 3: Build junctions
        junctions = build_junctions(modules, topology=topology)
        
        # Step 4: Design each junction
        junction_plans = []
        for junction in junctions:
            jp = self._design_junction(
                junction,
                module_att_sites,
                reaction_type
            )
            junction_plans.append(jp)
        
        plan.junction_plans = junction_plans
        
        # Check feasibility
        for jp in junction_plans:
            if jp.warnings:
                for w in jp.warnings:
                    plan.warnings.append(
                        f"J{jp.junction_index} ({jp.left_module_name}→{jp.right_module_name}): {w}"
                    )
        
        # Step 5: Simulate recombination to get product sequence and fragment annotations
        product_seq, byproduct_seq, byproduct_desc, fragment_annotations = self._simulate_recombination(
            modules,
            junction_plans,
            reaction_type
        )

        plan.product_sequence = product_seq
        plan.byproduct_sequence = byproduct_seq
        plan.byproduct_description = byproduct_desc
        
        # Step 6: Analyze ccdB selection
        ccdb_info = self._analyze_ccdb_selection(
            product_seq,
            byproduct_seq,
            reaction_type
        )

        # CRITICAL: In BP reaction, the Entry clone (product) should NOT contain ccdB
        # If ccdB is found in product, swap product and byproduct
        if reaction_type == "BP" and ccdb_info["ccdb_in_product"] and not ccdb_info["ccdb_in_byproduct"]:
            plan.warnings.append(
                "ccdB found in product - swapping product/byproduct (Entry clone should not contain ccdB)"
            )
            # Swap sequences
            product_seq, byproduct_seq = byproduct_seq, product_seq
            plan.product_sequence = product_seq
            plan.byproduct_sequence = byproduct_seq
            plan.byproduct_description = "Linear fragment with ccdB gene (byproduct)"

            # Re-analyze after swap
            ccdb_info = self._analyze_ccdb_selection(
                product_seq,
                byproduct_seq,
                reaction_type
            )
        elif ccdb_info["ccdb_in_product"]:
            plan.warnings.append(
                "ccdB found in product - recombination may have failed or "
                "product requires DB3.1 strain for propagation"
            )

        # Step 6b: Validate origin of replication in product
        ori_validation = self._validate_origin(product_seq)
        if not ori_validation["has_ori"]:
            plan.warnings.append(
                "WARNING: No origin of replication detected in product - plasmid may not replicate in E. coli"
            )
        
        # Step 7: Build primer table (for PCR addition of att sites)
        plan.primer_table = self._build_primer_table(junction_plans, module_att_sites)
        
        # Step 8: Generate visualization annotations
        plan.product_annotations = self._build_visualization(
            product_seq,
            junction_plans,
            ccdb_info,
            fragment_annotations,
            modules,
            plan.primer_table
        )
        
        # Step 9: Build BOM, steps, and metrics
        plan.bom = self._build_bom(junction_plans, reaction_type, ccdb_info)
        plan.steps = self._build_steps(junction_plans, reaction_type)
        plan.metrics = self._compute_metrics(junction_plans, reaction_type)
        
        # Step 10: Generate summary
        plan.summary = self._generate_summary(plan)
        
        return plan
    
    def _infer_reaction_type(self, module_att_sites: List[dict]) -> str:
        """
        Infer Gateway reaction type based on detected att sites.
        
        Rules:
        - If attB and attP found: BP reaction
        - If attL and attR found: LR reaction
        - If no att sites or only attB: BP reaction (will add attB via PCR)
        
        Args:
            module_att_sites: List of module data with att sites
            
        Returns:
            "BP" or "LR"
        """
        all_site_types = set()
        
        for mod_data in module_att_sites:
            for site in mod_data["att_sites"]:
                site_type, _ = parse_att_site_type(site.site_type)
                all_site_types.add(site_type)
        
        # Check for L and R sites (LR reaction)
        if "L" in all_site_types and "R" in all_site_types:
            return "LR"
        
        # Check for B and P sites (BP reaction)
        if "B" in all_site_types or "P" in all_site_types:
            return "BP"
        
        # Default to BP reaction
        return "BP"
    
    def _design_junction(
        self,
        junction: Junction,
        module_att_sites: List[dict],
        reaction_type: str
    ) -> GatewayJunctionPlan:
        """
        Design a single Gateway junction.
        
        Finds compatible att site pairs or designs PCR primers to add them.
        
        Args:
            junction: Junction metadata
            module_att_sites: Module att site data
            reaction_type: "BP" or "LR"
            
        Returns:
            GatewayJunctionPlan for this junction
        """
        left_idx = junction.left_module_index
        right_idx = junction.right_module_index
        
        left_data = module_att_sites[left_idx]
        right_data = module_att_sites[right_idx]
        
        # Find att sites near junction boundaries
        # For left module: look for sites near the end
        # For right module: look for sites near the start
        
        left_sites = self._find_boundary_sites(
            left_data["att_sites"],
            len(left_data["sequence"]),
            "end"
        )
        
        right_sites = self._find_boundary_sites(
            right_data["att_sites"],
            len(right_data["sequence"]),
            "start"
        )
        
        # Try to find compatible pair
        compatible_pair = None
        for ls in left_sites:
            for rs in right_sites:
                if validate_orthogonality(ls, rs):
                    compatible_pair = (ls, rs)
                    break
            if compatible_pair:
                break
        
        warnings = []
        
        if compatible_pair:
            # Native sites found
            left_site, right_site = compatible_pair
            strategy = "native_sites"
            
            # Determine product sites based on reaction type
            # Use the helper function to get correct products
            from .gateway_sites import get_recombination_products
            product_left, product_right = get_recombination_products(
                left_site.site_type,
                right_site.site_type
            )
            
            return GatewayJunctionPlan(
                junction_index=junction.junction_index,
                left_module_name=left_data["module_name"],
                right_module_name=right_data["module_name"],
                left_att_site=left_site.site_type,
                right_att_site=right_site.site_type,
                strategy=strategy,
                product_left_site=product_left,
                product_right_site=product_right,
                reaction_type=reaction_type,
                warnings=warnings
            )
        
        else:
            # No native sites found - need to add attB via PCR
            # In Gateway BP: insert gets attB sites, recombines with donor vector's attP sites
            strategy = "pcr_add_sites"

            # Get ALL att sites from both modules (not just boundaries)
            left_all_sites = left_data["att_sites"]
            right_all_sites = right_data["att_sites"]

            # Determine which att sites to add based on which module has native sites
            if reaction_type == "BP":
                # Determine which module is donor (has attP) and which is insert (needs attB)
                if len(left_all_sites) > 0 and len(right_all_sites) == 0:
                    # Left has sites (donor), right needs sites (insert)
                    donor_sites = left_all_sites
                    donor_is_left = True
                elif len(right_all_sites) > 0 and len(left_all_sites) == 0:
                    # Right has sites (donor), left needs sites (insert)
                    donor_sites = right_all_sites
                    donor_is_left = False
                else:
                    # Both or neither have sites - fall back to default
                    donor_sites = []
                    donor_is_left = None

                # Extract site numbers from donor att sites
                from .gateway_sites import parse_att_site_type
                site_numbers = set()
                site_by_number = {}  # Map site number to actual site name

                for site in donor_sites:
                    site_type, site_num = parse_att_site_type(site.site_type)
                    if site_num != "?":
                        site_numbers.add(site_num)
                        site_by_number[site_num] = site.site_type

                # Assign attB sites to insert with matching site numbers
                # For circular topology: junction 0 and 1 need different site numbers
                sorted_numbers = sorted(site_numbers)

                if len(sorted_numbers) >= 2:
                    # Use site number based on junction index
                    # Junction 0: use first site number (typically 1)
                    # Junction 1: use second site number (typically 2)
                    if junction.junction_index == 0:
                        site_num = sorted_numbers[0]
                    else:
                        site_num = sorted_numbers[1]
                elif len(sorted_numbers) == 1:
                    # Only one site number available, use it for both junctions
                    site_num = sorted_numbers[0]
                else:
                    # No donor sites found, use default numbering
                    site_num = "1" if junction.junction_index == 0 else "2"

                # Build junction att site names
                if donor_is_left:
                    # Left is donor (has attP), right is insert (needs attB)
                    left_site_name = site_by_number.get(site_num, f"attP{site_num}")
                    right_site_name = f"attB{site_num}"
                elif donor_is_left == False:
                    # Right is donor (has attP), left is insert (needs attB)
                    left_site_name = f"attB{site_num}"
                    right_site_name = site_by_number.get(site_num, f"attP{site_num}")
                else:
                    # Fall back to default numbering
                    if junction.junction_index == 0:
                        left_site_name = "attB1"
                        right_site_name = "attB2"
                    else:
                        left_site_name = "attB2"
                        right_site_name = "attB1"
                    warnings.append(
                        "No donor vector with attP sites detected. Using default attB site numbering."
                    )

                # Products from BP reaction: attB + attP → attL + attR
                from .gateway_sites import get_recombination_products
                product_left, product_right = get_recombination_products(
                    left_site_name,
                    right_site_name
                )

            else:  # LR reaction
                warnings.append(
                    "LR reaction requires native attL/attR sites - cannot add via PCR"
                )
                left_site_name = "unknown"
                right_site_name = "unknown"
                product_left = "unknown"
                product_right = "unknown"
            
            return GatewayJunctionPlan(
                junction_index=junction.junction_index,
                left_module_name=left_data["module_name"],
                right_module_name=right_data["module_name"],
                left_att_site=left_site_name,
                right_att_site=right_site_name,
                strategy=strategy,
                product_left_site=product_left,
                product_right_site=product_right,
                reaction_type=reaction_type,
                warnings=warnings
            )
    
    def _find_boundary_sites(
        self,
        sites: List[AttSiteMatch],
        seq_len: int,
        boundary: str
    ) -> List[AttSiteMatch]:
        """
        Find att sites near sequence boundaries.
        
        Args:
            sites: List of att sites in sequence
            seq_len: Length of sequence
            boundary: "start" or "end"
            
        Returns:
            List of sites within 500bp of boundary
        """
        boundary_sites = []
        threshold = 500
        
        for site in sites:
            if boundary == "start":
                if site.start < threshold:
                    boundary_sites.append(site)
            else:  # end
                if site.end > seq_len - threshold:
                    boundary_sites.append(site)
        
        return boundary_sites
    
    def _simulate_recombination(
        self,
        modules: List[dict],
        junction_plans: List[GatewayJunctionPlan],
        reaction_type: str
    ) -> Tuple[str, str, str, List[FragmentAnnotation]]:
        """
        Simulate Gateway recombination to generate product and byproduct sequences.

        For BP reaction:
            Input: Insert (attB1-insert-attB2) + Donor (backbone_L + attP1-ccdB-attP2 + backbone_R)
            Product: backbone_L + attL1 + insert + attL2 + backbone_R (NO ccdB)
            Byproduct: attR1 + ccdB + attR2 (linear)

        Args:
            modules: Module list
            junction_plans: Junction design plans
            reaction_type: "BP" or "LR"

        Returns:
            (product_sequence, byproduct_sequence, byproduct_description, fragment_annotations)
        """
        if len(junction_plans) == 0:
            return "", "", "No junctions", []

        # For simple 2-module case
        if len(modules) == 2:
            jp = junction_plans[0]

            # Scan both modules for att sites to identify insert vs donor
            mod0_seq = modules[0].get("sequence", "")
            mod1_seq = modules[1].get("sequence", "")

            mod0_sites = scan_att_sites(mod0_seq, fuzzy_threshold=0, search_attB_only=False)
            mod1_sites = scan_att_sites(mod1_seq, fuzzy_threshold=0, search_attB_only=False)

            # Identify which module is donor (has attP) vs insert (has attB)
            mod0_has_attP = any("attP" in s.site_type for s in mod0_sites)
            mod1_has_attP = any("attP" in s.site_type for s in mod1_sites)

            if reaction_type == "BP":
                # BP reaction: attB (insert) + attP (donor) → attL (product) + attR (byproduct)
                if mod0_has_attP and not mod1_has_attP:
                    donor_seq = mod0_seq
                    donor_sites = mod0_sites
                    donor_module = modules[0]
                    insert_seq = mod1_seq
                    insert_module = modules[1]
                elif mod1_has_attP and not mod0_has_attP:
                    donor_seq = mod1_seq
                    donor_sites = mod1_sites
                    donor_module = modules[1]
                    insert_seq = mod0_seq
                    insert_module = modules[0]
                else:
                    # Fallback: simple concatenation if site detection fails
                    product = GATEWAY_ATT_SITES.get(jp.product_left_site, "") + mod0_seq + \
                              GATEWAY_ATT_SITES.get(jp.product_right_site, "") + mod1_seq
                    return product, "", "Byproduct (site detection failed)", []

                # Get module names for annotations
                donor_name = donor_module.get("canonical_id", "Donor")
                insert_name = insert_module.get("canonical_id", "Insert")

                # Find attP sites in donor vector to extract regions
                attP_sites = [s for s in donor_sites if "attP" in s.site_type]
                attP_sites.sort(key=lambda s: s.start)

                if len(attP_sites) >= 2:
                    # Extract donor backbone regions (exclude ccdB cassette between attP sites)
                    first_attP = attP_sites[0]
                    second_attP = attP_sites[1]

                    # Donor backbone left: from start to first attP start
                    donor_left_backbone = donor_seq[:first_attP.start]

                    # Region between attP sites (contains ccdB) - this goes to byproduct
                    ccdb_region = donor_seq[first_attP.end:second_attP.start]

                    # Donor backbone right: from second attP end to end of sequence
                    donor_right_backbone = donor_seq[second_attP.end:]

                    # Build product: donor_left + attL1 + insert + attL2 + donor_right
                    attL1_seq = GATEWAY_ATT_SITES.get(jp.product_left_site, "")
                    attL2_seq = GATEWAY_ATT_SITES.get(jp.product_right_site, "")

                    product = donor_left_backbone + attL1_seq + insert_seq + attL2_seq + donor_right_backbone

                    # Build fragment annotations
                    fragments = []
                    pos = 0

                    # Fragment 1: Donor left backbone
                    if len(donor_left_backbone) > 0:
                        fragments.append(FragmentAnnotation(
                            name=f"{donor_name} (5' backbone)",
                            start=pos,
                            end=pos + len(donor_left_backbone),
                            source_module=donor_name,
                            source_type="donor_backbone",
                            direction=1
                        ))
                        pos += len(donor_left_backbone)

                    # Fragment 2: attL1 site
                    fragments.append(FragmentAnnotation(
                        name=jp.product_left_site,
                        start=pos,
                        end=pos + len(attL1_seq),
                        source_module="Gateway Recombination",
                        source_type="att_site",
                        direction=1
                    ))
                    pos += len(attL1_seq)

                    # Fragment 3: Insert
                    fragments.append(FragmentAnnotation(
                        name=insert_name,
                        start=pos,
                        end=pos + len(insert_seq),
                        source_module=insert_name,
                        source_type="insert",
                        direction=1
                    ))
                    pos += len(insert_seq)

                    # Fragment 4: attL2 site
                    fragments.append(FragmentAnnotation(
                        name=jp.product_right_site,
                        start=pos,
                        end=pos + len(attL2_seq),
                        source_module="Gateway Recombination",
                        source_type="att_site",
                        direction=1
                    ))
                    pos += len(attL2_seq)

                    # Fragment 5: Donor right backbone
                    if len(donor_right_backbone) > 0:
                        fragments.append(FragmentAnnotation(
                            name=f"{donor_name} (3' backbone)",
                            start=pos,
                            end=pos + len(donor_right_backbone),
                            source_module=donor_name,
                            source_type="donor_backbone",
                            direction=1
                        ))

                    # Build byproduct: attR1 + ccdB_region + attR2 (linear)
                    # Get attR site names from junction plan
                    attR1_name = jp.left_att_site.replace("attB", "attR").replace("attP", "attR")
                    attR2_name = jp.right_att_site.replace("attB", "attR").replace("attP", "attR")
                    attR1_seq = GATEWAY_ATT_SITES.get(attR1_name, "")
                    attR2_seq = GATEWAY_ATT_SITES.get(attR2_name, "")

                    byproduct = attR1_seq + ccdb_region + attR2_seq
                    byproduct_desc = f"Linear fragment with {attR1_name} and {attR2_name} sites (contains ccdB)"

                else:
                    # Fallback if attP sites not found properly
                    attL_seq = GATEWAY_ATT_SITES.get(jp.product_left_site, "")
                    product = donor_seq + attL_seq + insert_seq
                    byproduct = ""
                    byproduct_desc = "Byproduct (attP site detection incomplete)"
                    fragments = []

            else:
                # LR reaction: attL + attR → attB + attP
                # Product = destination_vector with insert replacing ccdB
                left_product_seq = GATEWAY_ATT_SITES.get(jp.product_left_site, "")
                right_product_seq = GATEWAY_ATT_SITES.get(jp.product_right_site, "")
                product = left_product_seq + mod0_seq + right_product_seq + mod1_seq
                byproduct = "Entry clone backbone with attP sites"
                byproduct_desc = "Entry clone backbone with attP sites (byproduct)"
                fragments = []  # TODO: Implement LR fragment tracking

            return product, byproduct, byproduct_desc, fragments

        # For multi-module: concatenate all with att sites
        product_parts = []
        for i, mod in enumerate(modules):
            if i < len(junction_plans):
                jp = junction_plans[i]
                product_parts.append(GATEWAY_ATT_SITES.get(jp.product_left_site, ""))
            product_parts.append(mod.get("sequence", ""))

        if junction_plans:
            last_jp = junction_plans[-1]
            product_parts.append(GATEWAY_ATT_SITES.get(last_jp.product_right_site, ""))

        product = "".join(product_parts)

        return product, "", "Byproduct (linearized backbone)", []
    
    def _analyze_ccdb_selection(
        self,
        product_seq: str,
        byproduct_seq: str,
        reaction_type: str
    ) -> dict:
        """
        Analyze ccdB selection in Gateway cloning.
        
        Args:
            product_seq: Product sequence
            byproduct_seq: Byproduct sequence
            reaction_type: "BP" or "LR"
            
        Returns:
            Dict with ccdB analysis results
        """
        ccdb_in_product = len(scan_for_ccdb(product_seq)) > 0
        ccdb_in_byproduct = len(scan_for_ccdb(byproduct_seq)) > 0 if byproduct_seq else False
        
        info = {
            "ccdb_in_product": ccdb_in_product,
            "ccdb_in_byproduct": ccdb_in_byproduct,
            "selection_strain": "DH5α or other standard E. coli",
            "notes": []
        }
        
        if reaction_type == "LR":
            if ccdb_in_product:
                info["notes"].append(
                    "WARNING: ccdB found in LR product - recombination likely failed"
                )
                info["selection_strain"] = "DB3.1 (ccdB resistant) for propagation"
            else:
                info["notes"].append(
                    "Successful LR recombination - ccdB cassette replaced with insert"
                )
        elif reaction_type == "BP":
            if ccdb_in_product:
                info["notes"].append(
                    "Entry clone contains ccdB - use DB3.1 strain for propagation"
                )
                info["selection_strain"] = "DB3.1 (ccdB resistant)"
            else:
                info["notes"].append(
                    "No ccdB in product - standard E. coli strains can be used"
                )
        
        return info

    def _validate_origin(self, sequence: str) -> dict:
        """
        Validate that the sequence contains an origin of replication.

        Searches for common bacterial origins (ColE1, pMB1, p15A, pSC101, etc.)

        Args:
            sequence: DNA sequence to validate

        Returns:
            Dict with origin validation results
        """
        seq_upper = sequence.upper()

        # Common E. coli origin signatures (partial sequences)
        # These are conserved regions from common origins
        common_origins = {
            "ColE1/pMB1": [
                "TTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCTTGCAAACAAAAAAACCACCGCTAC",
                "CCTTTCGGGAGTTGTG",
                "GGCCCTCGATATACTT"
            ],
            "p15A": [
                "CTGTAAAAGCCGGC",
                "GCATGGGTGCGCTG"
            ],
            "pSC101": [
                "GTTTAAACGGTCTCCAGCTTGGCTGTTTTGGCGGATG"
            ],
            "pBR322": [
                "GAGATCCTTTTTTTCTGCGCGTAATCTGCTGCTTGCAAACAA"
            ],
            "R6K": [
                "GCGCGCCTCGTTCATTCACGTTTTTGAACCCGTGGAGGACGGGCAGAC"
            ]
        }

        detected_origins = []
        for ori_name, signatures in common_origins.items():
            for sig in signatures:
                if sig in seq_upper:
                    detected_origins.append(ori_name)
                    break

        # Also check for common ori-related keywords in annotations
        # (this is a heuristic check - actual origin may be present even if not detected)
        has_ori = len(detected_origins) > 0

        return {
            "has_ori": has_ori,
            "detected_origins": detected_origins,
            "note": "Origin detection based on conserved sequences" if has_ori else "No common origin signatures detected"
        }

    def _build_primer_table(
        self,
        junction_plans: List[GatewayJunctionPlan],
        module_att_sites: List[dict]
    ) -> List[Dict[str, Any]]:
        """
        Build primer table for PCR addition of att sites.
        
        Identifies modules that lack native att sites and designs BOTH FWD and REV
        primers for those modules. Modules with existing att sites (like pDONR vectors)
        do not need primers.
        
        Args:
            junction_plans: Junction designs
            module_att_sites: Module att site data
            
        Returns:
            List of primer dicts
        """
        primers = []
        processed_modules = set()
        
        # First pass: identify which modules need primers
        modules_needing_primers = {}

        # Detect expression modules to determine insert orientation
        expression_modules = {}
        for mod_data in module_att_sites:
            module_name = mod_data["module_name"]
            module_seq = mod_data["sequence"]
            expr_info = self._detect_expression_module(module_seq, module_name)
            if expr_info["is_expression_module"]:
                expression_modules[module_name] = expr_info
                logger.info(f"[Gateway] Detected expression module: {module_name}, "
                          f"promoter: {expr_info['promoter']}, "
                          f"orientation: {expr_info['orientation']}")

        
        for jp in junction_plans:
            if jp.strategy != "pcr_add_sites":
                continue
            
            # Find module data
            for mod_data in module_att_sites:
                if mod_data["module_name"] == jp.left_module_name:
                    if len(mod_data["att_sites"]) == 0 and jp.left_module_name not in modules_needing_primers:
                        modules_needing_primers[jp.left_module_name] = {
                            "module_data": mod_data,
                            "attB_start": None,  # Will be determined from junctions
                            "attB_end": None
                        }
                
                if mod_data["module_name"] == jp.right_module_name:
                    if len(mod_data["att_sites"]) == 0 and jp.right_module_name not in modules_needing_primers:
                        modules_needing_primers[jp.right_module_name] = {
                            "module_data": mod_data,
                            "attB_start": None,
                            "attB_end": None
                        }
        
        # Second pass: determine which attB sites to add
        for jp in junction_plans:
            if jp.strategy != "pcr_add_sites":
                continue
            
            # Track which attB sites to add to each module
            if jp.left_module_name in modules_needing_primers:
                # Left module needs attB at its END
                if jp.left_att_site.startswith("attB"):
                    modules_needing_primers[jp.left_module_name]["attB_end"] = jp.left_att_site
            
            if jp.right_module_name in modules_needing_primers:
                # Right module needs attB at its START
                if jp.right_att_site.startswith("attB"):
                    modules_needing_primers[jp.right_module_name]["attB_start"] = jp.right_att_site
        
        # Third pass: design primers for each module that needs them
        for module_name, info in modules_needing_primers.items():
            mod_data = info["module_data"]
            seq = mod_data["sequence"]
            
            # Check if we should orient this module based on expression context
            oriented_seq = seq
            was_reversed = False

            if expression_modules:
                # Check if this module is itself an expression module
                if module_name not in expression_modules:
                    # This is an insert - orient it to match expression module
                    expr_module_name = list(expression_modules.keys())[0]
                    expr_orientation = expression_modules[expr_module_name]["orientation"]

                    if expr_orientation != "unknown":
                        oriented_seq, was_reversed = self._orient_insert_for_expression(
                            seq,
                            expr_orientation
                        )
                        if was_reversed:
                            logger.info(f"[Gateway] Oriented insert {module_name} to match "
                                      f"expression module {expr_module_name} ({expr_orientation})")

                        # Update seq with oriented version
                        seq = oriented_seq
            anneal_len = 20
            
            # Design FORWARD primer (adds attB to START of sequence)
            if info["attB_start"] and len(seq) >= anneal_len:
                attB_seq = GATEWAY_ATT_SITES.get(info["attB_start"], "")
                anneal_region = seq[:anneal_len]
                # Add Ribosome Binding Site (RBS) for expression
                rbs_seq = ""  # Default: no RBS
                
                # Get the expression module info if available
                if expression_modules:
                    expr_module_name = list(expression_modules.keys())[0]
                    expr_info = expression_modules[expr_module_name]
                    promoter_name = expr_info.get("promoter", "unknown")
                    is_bacterial = expr_info.get("is_bacterial", None)
                    
                    if is_bacterial is True:
                        rbs_seq = "CAAGGAGG"  # Shine-Dalgarno
                        logger.info(f"[Gateway] Adding Shine-Dalgarno for bacterial expression ({promoter_name})")
                    elif is_bacterial is False:
                        rbs_seq = "GCCACC"  # Kozak consensus
                        logger.info(f"[Gateway] Adding Kozak for eukaryotic expression ({promoter_name})")
                    else:
                        logger.info(f"[Gateway] Unknown promoter type, no RBS added")
                else:
                    logger.info(f"[Gateway] No expression module detected, no RBS added")


                # Construct FWD primer: GGGG + RBS + attB + CC + annealing

                fwd_primer = "GGGG" + rbs_seq + attB_seq + "CC" + anneal_region
                
                primers.append({
                    "primer_name": f"{module_name}_attB_FWD",
                    "sequence": fwd_primer,
                    "length": len(fwd_primer),
                    "att_site_tail": attB_seq,
                    "annealing_region": anneal_region,
                    "purpose": f"Add {info['attB_start']} to {module_name}",
                    "tm_anneal": self._calc_tm(anneal_region) if self._thermo else 60.0
                })
            
            # Design REVERSE primer (adds attB to END of sequence via RC)
            if info["attB_end"] and len(seq) >= anneal_len:
                attB_seq_rc = _reverse_complement(GATEWAY_ATT_SITES.get(info["attB_end"], ""))
                anneal_region_rc = _reverse_complement(seq[-anneal_len:])
                rev_primer = "GGGG" + attB_seq_rc + "C" + anneal_region_rc  # Add 1 nucleotide for frame maintenance  # Add GGGG tail per Gateway manual
                
                primers.append({
                    "primer_name": f"{module_name}_attB_REV",
                    "sequence": rev_primer,
                    "length": len(rev_primer),
                    "att_site_tail": attB_seq_rc,
                    "annealing_region": anneal_region_rc,
                    "purpose": f"Add {info['attB_end']} to {module_name}",
                    "tm_anneal": self._calc_tm(anneal_region_rc) if self._thermo else 60.0
                })
        
        return primers


    def _build_visualization(
        self,
        product_seq: str,
        junction_plans: List[GatewayJunctionPlan],
        ccdb_info: dict,
        fragment_annotations: List[FragmentAnnotation],
        modules: List[dict],
        primer_table: List[Dict[str, Any]] = None
    ) -> List[Dict[str, Any]]:
        """
        Build SeqViz-compatible annotations for product visualization.

        Includes:
        1. Fragment annotations showing source of each region
        2. Preserved features from original input modules
        3. att site annotations
        4. ccdB warning annotation if present

        Args:
            product_seq: Product sequence
            junction_plans: Junction designs
            ccdb_info: ccdB analysis results
            fragment_annotations: Fragment source tracking
            modules: Original input modules

        Returns:
            List of annotation dicts
        """
        annotations = []

        # 1. Add fragment annotations (like inv_gib visualization)
        fragment_colors = {
            "donor_backbone": "#A8DADC",  # Light blue
            "insert": "#457B9D",  # Blue
            "att_site": "#F38181",  # Salmon (for attL sites)
        }

        for frag in fragment_annotations:
            annotations.append({
                "name": frag.name,
                "start": frag.start,
                "end": frag.end,
                "direction": frag.direction,
                "type": frag.source_type,
                "color": fragment_colors.get(frag.source_type, "#CCCCCC"),
                "source": frag.source_module,
                "metadata": {
                    "source_type": frag.source_type
                }
            })

        # 2. Extract and preserve features from original modules
        # Map module sequences to product positions using fragment annotations
        for frag in fragment_annotations:
            if frag.source_type in ("donor_backbone", "insert"):
                # Find the source module
                source_module = None
                for mod in modules:
                    mod_id = mod.get("canonical_id", "")
                    if mod_id in frag.source_module:
                        source_module = mod
                        break

                if source_module and "features" in source_module:
                    # Get fragment sequence from product
                    frag_seq = product_seq[frag.start:frag.end]

                    # Preserve features that fall within this fragment
                    for feature in source_module.get("features", []):
                        feat_start = feature.get("start", 0)
                        feat_end = feature.get("end", 0)
                        feat_type = feature.get("type", "misc_feature")
                        feat_name = feature.get("name", "")

                        # Check if feature is fully within the fragment boundaries
                        # (Simple mapping - assumes features are in extracted region)
                        # TODO: More sophisticated coordinate mapping for circular plasmids
                        if feat_start >= 0 and feat_end <= len(frag_seq):
                            # Map feature coordinates to product
                            product_start = frag.start + feat_start
                            product_end = frag.start + feat_end

                            if product_end <= len(product_seq):
                                annotations.append({
                                    "name": feat_name,
                                    "start": product_start,
                                    "end": product_end,
                                    "type": feat_type,
                                    "direction": feature.get("direction", 1),
                                    "color": feature.get("color", "#999999"),
                                    "source": frag.source_module,
                                    "metadata": {
                                        "preserved_from": frag.source_module,
                                        "original_feature": True
                                    }
                                })

        # 3. Scan product for att sites (for verification/labeling)
        product_sites = scan_att_sites(product_seq, fuzzy_threshold=0)

        # Color mapping for att sites
        att_colors = {
            "B": "#4ECDC4",  # Teal
            "P": "#95E1D3",  # Light teal
            "L": "#F38181",  # Salmon
            "R": "#AA96DA",  # Purple
        }

        for site in product_sites:
            site_type, site_num = parse_att_site_type(site.site_type)
            color = att_colors.get(site_type, "#CCCCCC")

            # Only add if not already covered by fragment annotation
            already_annotated = any(
                a["start"] == site.start and a["end"] == site.end
                for a in annotations
            )

            if not already_annotated:
                annotations.append({
                    "start": site.start,
                    "end": site.end,
                    "name": site.site_type,
                    "type": "recombination_site",
                    "color": color,
                    "direction": 1 if site.strand == 1 else -1,
                    "metadata": {
                        "core_sequence": site.core_sequence,
                        "match_quality": site.match_quality
                    }
                })

        # 4. Add ccdB annotations if present (WARNING - should not be in product)
        if ccdb_info["ccdb_in_product"]:
            ccdb_sites = scan_for_ccdb(product_seq)
            for start, end in ccdb_sites:
                annotations.append({
                    "start": start,
                    "end": end,
                    "name": "ccdB (BYPRODUCT MARKER)",
                    "type": "cds",
                    "color": "#8B0000",  # Dark red
                    "direction": 1,
                    "metadata": {
                        "note": "WARNING: ccdB in product indicates byproduct, not Entry clone"
                    }
                })


        # 4. Add primer binding site annotations
        if primer_table:
            # Find insert fragments in product
            insert_fragments = [frag for frag in fragment_annotations if frag.source_type == "insert"]

            for primer in primer_table:
                primer_name = primer.get("primer_name", "")
                annealing_region = primer.get("annealing_region", "")

                if not annealing_region:
                    continue

                # Determine if this is FWD or REV primer
                is_fwd = "_FWD" in primer_name or "FWD" in primer_name

                # Find where this primer binds in the product
                # For FWD primer: binds to start of insert
                # For REV primer: binds to end of insert (reverse complement)
                for frag in insert_fragments:
                    frag_seq = product_seq[frag.start:frag.end]

                    if is_fwd:
                        # FWD primer: check if annealing region matches start of insert
                        if frag_seq.startswith(annealing_region):
                            annotations.append({
                                "name": f"{primer_name} binding site",
                                "start": frag.start,
                                "end": frag.start + len(annealing_region),
                                "direction": 1,
                                "type": "primer_bind",
                                "color": "#FF6B9D",  # Pink for primers
                                "metadata": {
                                    "primer_name": primer_name,
                                    "primer_type": "forward",
                                    "annealing_length": len(annealing_region)
                                }
                            })
                    else:
                        # REV primer: check if RC of annealing region matches end of insert
                        # The annealing_region is already the RC in the primer table
                        expected_seq = _reverse_complement(annealing_region)
                        if frag_seq.endswith(expected_seq):
                            annotations.append({
                                "name": f"{primer_name} binding site",
                                "start": frag.end - len(expected_seq),
                                "end": frag.end,
                                "direction": -1,
                                "type": "primer_bind",
                                "color": "#FF6B9D",  # Pink for primers
                                "metadata": {
                                    "primer_name": primer_name,
                                    "primer_type": "reverse",
                                    "annealing_length": len(expected_seq)
                                }
                            })

        return annotations
    
    def _build_bom(
        self,
        junction_plans: List[GatewayJunctionPlan],
        reaction_type: str,
        ccdb_info: dict
    ) -> List[str]:
        """Build bill of materials."""
        bom = []
        
        # Count primers needed
        pcr_junctions = [jp for jp in junction_plans if jp.strategy == "pcr_add_sites"]
        primer_count = len(pcr_junctions) * 2
        
        if primer_count > 0:
            bom.append(f"{primer_count} custom primers for attB site addition")
            bom.append(f"PCR reagents for {len(pcr_junctions)} reactions")
            bom.append(f"Gel purification kit")
        
        # Gateway enzyme mix
        enzyme_name = "BP Clonase II" if reaction_type == "BP" else "LR Clonase II"
        bom.append(f"{enzyme_name} enzyme mix")
        
        # Vectors
        if reaction_type == "BP":
            bom.append("pDONR vector (with attP sites)")
        else:
            bom.append("Destination vector (with attR-ccdB-attR cassette)")
        
        # Competent cells
        strain = ccdb_info.get("selection_strain", "DH5α")
        bom.append(f"{strain} competent cells")
        
        # Standard reagents
        bom.append("LB agar plates with appropriate antibiotic")
        bom.append("Miniprep kit")
        bom.append("Sequencing primers (2× for verification)")
        
        return bom
    
    def _build_steps(
        self,
        junction_plans: List[GatewayJunctionPlan],
        reaction_type: str
    ) -> List[BuildStep]:
        """Build step-by-step protocol."""
        steps = []
        step_num = 1
        
        # PCR steps if needed
        pcr_junctions = [jp for jp in junction_plans if jp.strategy == "pcr_add_sites"]
        
        if pcr_junctions:
            steps.append(BuildStep(
                step_number=step_num,
                step_type="pcr",
                description=f"PCR amplify {len(pcr_junctions)} module(s) to add attB sites",
                materials=["Template DNA", "Custom primers with attB tails", "High-fidelity polymerase"],
                estimated_hours=0.5,
                estimated_days=0.25
            ))
            step_num += 1
            
            steps.append(BuildStep(
                step_number=step_num,
                step_type="gel_purification",
                description="Gel purify PCR products",
                materials=["Agarose gel", "Gel extraction kit"],
                estimated_hours=1.0,
                estimated_days=0.25
            ))
            step_num += 1
        
        # Gateway recombination
        reaction_name = "BP" if reaction_type == "BP" else "LR"
        enzyme_name = "BP Clonase II" if reaction_type == "BP" else "LR Clonase II"
        
        steps.append(BuildStep(
            step_number=step_num,
            step_type="assembly",
            description=f"Gateway {reaction_name} recombination reaction",
            materials=[
                f"{enzyme_name}",
                "Insert DNA (with attB or attL sites)",
                "Vector DNA (with attP or attR sites)",
                "Reaction buffer"
            ],
            estimated_hours=1.5,
            estimated_days=0.5
        ))
        step_num += 1
        
        # Transformation
        steps.append(BuildStep(
            step_number=step_num,
            step_type="transformation",
            description="Transform into competent cells and plate",
            materials=["Competent cells", "SOC medium", "LB agar plates with antibiotic"],
            estimated_hours=0.5,
            estimated_days=1.0
        ))
        step_num += 1
        
        # Colony screening
        steps.append(BuildStep(
            step_number=step_num,
            step_type="colony_screening",
            description="Pick colonies and miniprep",
            materials=["Sterile tips", "LB medium", "Miniprep kit"],
            estimated_hours=2.0,
            estimated_days=0.5
        ))
        step_num += 1
        
        # Sequencing
        steps.append(BuildStep(
            step_number=step_num,
            step_type="sequencing",
            description="Sanger sequencing verification of junctions",
            materials=["Sequencing primers", "Miniprep DNA"],
            estimated_hours=0.25,
            estimated_days=2.0
        ))
        
        return steps
    
    def _compute_metrics(
        self,
        junction_plans: List[GatewayJunctionPlan],
        reaction_type: str
    ) -> OperatorMetrics:
        """Compute cost, time, and risk metrics."""
        metrics = OperatorMetrics()
        
        # Count PCR reactions needed
        pcr_junctions = [jp for jp in junction_plans if jp.strategy == "pcr_add_sites"]
        pcr_count = len(pcr_junctions)
        primer_count = pcr_count * 2
        
        # Costs
        if primer_count > 0:
            metrics.primer_count = primer_count
            # Gateway attB primers: tail (~25 nt attB1/attB2) + 18-25 nt anneal
            # = ~45 mer. Use length-scaled $0.24/bp via lab.primer_cost.
            metrics.primer_cost_usd = primer_count * self.lab.primer_cost(45)
            metrics.pcr_count = pcr_count
            metrics.pcr_cost_usd = pcr_count * self.lab.pcr_rxn_cost_usd
            metrics.gel_count = pcr_count
            metrics.gel_cost_usd = pcr_count * self.lab.gel_lane_cost_usd
        
        # Gateway enzyme cost
        # BP or LR Clonase II (whichever the workflow needs); both priced
        # the same via lab catalog so the per-rxn line is unchanged.
        metrics.assembly_cost_usd = self.lab.catalog["bp_clonase_ii"].cost_per_rxn
        
        # Transformation and screening
        metrics.transformation_count = 1
        metrics.transformation_cost_usd = self.lab.transformation_cost_usd
        metrics.miniprep_count = 2  # Screen 2 colonies
        # Miniprep folded into ONT plasmid sequencing service.
        metrics.miniprep_cost_usd = 0.0
        
        # Sequencing (2 junctions)
        metrics.sequencing_count = self.lab.sequencing_reads_per_construct
        metrics.sequencing_cost_usd = metrics.sequencing_count * self.lab.sequencing_cost_usd
        
        # Total cost
        metrics.total_cost_usd = (
            metrics.primer_cost_usd +
            metrics.pcr_cost_usd +
            metrics.gel_cost_usd +
            metrics.assembly_cost_usd +
            metrics.transformation_cost_usd +
            metrics.miniprep_cost_usd +
            metrics.sequencing_cost_usd
        )
        
        # Time estimates
        if pcr_count > 0:
            metrics.total_calendar_days = 4.5  # PCR + Gateway + transform + screen + seq
        else:
            metrics.total_calendar_days = 4.0  # No PCR needed
        
        metrics.total_labor_hours = 3.0 + (pcr_count * 1.5)
        
        # Risk assessment
        if pcr_count > 0:
            metrics.pcr_risk = 0.2  # PCR addition of att sites
        else:
            metrics.pcr_risk = 0.0
        
        metrics.assembly_risk = 0.1  # Gateway is very efficient (>95%)
        
        # Overall risk
        metrics.overall_risk_score = (metrics.pcr_risk + metrics.assembly_risk) / 2
        
        # Add warning count to risk
        for jp in junction_plans:
            if jp.warnings:
                metrics.risk_flags.extend(jp.warnings)
                metrics.overall_risk_score += len(jp.warnings) * 0.05
        
        metrics.overall_risk_score = min(1.0, metrics.overall_risk_score)
        
        return metrics
    
    def _generate_summary(self, plan: GatewayBuildPlan) -> str:
        """Generate human-readable summary."""
        reaction_name = "BP" if plan.reaction_type == "BP" else "LR"
        
        native_count = sum(1 for jp in plan.junction_plans if jp.strategy == "native_sites")
        pcr_count = sum(1 for jp in plan.junction_plans if jp.strategy == "pcr_add_sites")
        
        summary_parts = [
            f"Gateway {reaction_name} Cloning ({plan.fragment_count} modules)",
            f"  • {native_count} junction(s) with native att sites",
        ]
        
        if pcr_count > 0:
            summary_parts.append(f"  • {pcr_count} junction(s) requiring PCR addition of attB sites")
        
        summary_parts.extend([
            f"  • Estimated cost: ${plan.metrics.total_cost_usd:.2f}",
            f"  • Estimated time: {plan.metrics.total_calendar_days:.1f} days",
            f"  • Risk score: {plan.metrics.overall_risk_score:.2f}"
        ])
        
        return "\n".join(summary_parts)
    
    def _calc_tm(self, seq: str) -> float:
        """Calculate melting temperature."""
        if self._thermo:
            try:
                return self._thermo.calc_tm(seq)
            except Exception:
                pass
        
        # Fallback: simple GC content-based estimate
        gc = sum(1 for b in seq.upper() if b in "GC") / max(len(seq), 1)
        return 64.0 + 41.0 * gc - 500.0 / max(len(seq), 1)


# Helper functions

def _module_name(module: dict) -> str:
    """Extract module name from module dict."""
    return module.get("canonical_id") or module.get("description") or "Unknown"


def _reverse_complement(seq: str) -> str:
    """Reverse complement a DNA sequence."""
    complement = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    try:
        from Bio.Seq import Seq
        return str(Seq(seq).reverse_complement())
    except Exception:
        # Fallback manual implementation
        return "".join(complement.get(b.upper(), "N") for b in reversed(seq))
