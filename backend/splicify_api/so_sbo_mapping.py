"""
SO / SBO role mapping for hierarchical annotation module types and pLannotate
feature types.

Maps the ad-hoc `module_type` / feature-type strings emitted by the
rule-based detector, CDS submodule parser, mammalian Pol2 detector, and
pLannotate to their Sequence Ontology (SO) role and Systems Biology
Ontology (SBO) term equivalents.

- SO URIs are of the form  http://identifiers.org/so/SO:0000167  (promoter).
- SBO URIs are of the form http://identifiers.org/biomodels.sbo/SBO:0000459
  (stimulator participant role), etc.

This lets downstream code attach canonical role URIs on every emitted
annotation and lets the SBOL3 exporter round-trip role semantics without
rebuilding the mapping.
"""

from typing import Dict, List, Optional

# --------------------------------------------------------------------------- #
# Sequence Ontology (SO) role URIs
# --------------------------------------------------------------------------- #

SO = "http://identifiers.org/so/SO:"

SO_PROMOTER                 = SO + "0000167"
SO_CDS                      = SO + "0000316"
SO_TERMINATOR               = SO + "0000141"
SO_POLYA_SIGNAL             = SO + "0000551"
SO_RBS                      = SO + "0000139"
SO_KOZAK                    = SO + "0001647"
SO_OPERATOR                 = SO + "0000057"
SO_ENHANCER                 = SO + "0000165"
SO_INSULATOR                = SO + "0000627"
SO_ORIGIN_OF_REPLICATION    = SO + "0000296"
SO_ORIGIN_OF_TRANSFER       = SO + "0000724"
SO_LTR                      = SO + "0000286"
SO_ITR                      = SO + "0000481"
SO_POLYPEPTIDE_REGION       = SO + "0000839"
SO_LOCALIZATION_SIGNAL      = SO + "0001527"
SO_NLS                      = SO + "0002001"  # nuclear_localization_signal (community extension)
SO_EPITOPE_TAG              = SO + "0001668"
SO_PROTEIN_DOMAIN           = SO + "0000417"
SO_RIBOSOMAL_SKIP           = SO + "0001891"  # 2A self-cleaving peptide
SO_FLEXIBLE_LINKER          = SO + "0001828"  # polypeptide_linker
SO_GENE                     = SO + "0000704"
SO_REGULATORY_REGION        = SO + "0005836"
SO_MCS                      = SO + "0001957"  # polylinker_site
SO_RECOMBINATION_SITE       = SO + "0000299"  # specific_recombination_site
SO_ATTB                     = SO + "0000946"
SO_ATTP                     = SO + "0000947"
SO_ATTL                     = SO + "0000948"
SO_ATTR                     = SO + "0000949"
SO_LOXP                     = SO + "0000949"  # fallback — see inversion_site
SO_FRT                      = SO + "0000350"  # FRT_site
SO_INVERTED_REPEAT          = SO + "0000294"
SO_DIRECT_REPEAT            = SO + "0000314"
SO_TRANSPOSABLE_ELEMENT     = SO + "0000101"
SO_TRANSPOSON_END           = SO + "0000197"  # uORF stand-in; specialized per transposon
SO_SGRNA                    = SO + "0001998"
SO_TRACR_RNA                = SO + "0002098"
SO_GRNA_SCAFFOLD            = SO + "0001998"  # fallback to sgRNA
SO_SHRNA                    = SO + "0001997"
SO_PRIMER_BINDING_SITE      = SO + "0005850"
SO_PROTEIN_COMPOSITE        = SO + "0000839"  # polypeptide_region
SO_REGION                   = SO + "0000001"
SO_ENGINEERED_REGION        = SO + "0000804"

# --------------------------------------------------------------------------- #
# Systems Biology Ontology (SBO) terms for interaction types / roles
# --------------------------------------------------------------------------- #

SBO = "http://identifiers.org/biomodels.sbo/SBO:"

# Interaction types
SBO_GENETIC_PRODUCTION      = SBO + "0000589"  # genetic production
SBO_TRANSCRIPTION           = SBO + "0000183"
SBO_TRANSLATION             = SBO + "0000184"
SBO_STIMULATION             = SBO + "0000170"  # used for activation
SBO_INHIBITION              = SBO + "0000169"  # used for repression
SBO_CLEAVAGE                = SBO + "0000178"
SBO_NON_COVALENT_BINDING    = SBO + "0000177"
SBO_CONTROL                 = SBO + "0000168"
SBO_CONVERSION              = SBO + "0000182"  # for recombination
SBO_DEGRADATION             = SBO + "0000179"

# Participation roles
SBO_STIMULATOR              = SBO + "0000459"  # stimulator
SBO_INHIBITOR               = SBO + "0000020"
SBO_PROMOTER_ROLE           = SBO + "0000598"  # promoter (regulator)
SBO_TEMPLATE                = SBO + "0000645"  # template
SBO_PRODUCT                 = SBO + "0000011"
SBO_REACTANT                = SBO + "0000010"
SBO_MODIFIER                = SBO + "0000019"

# --------------------------------------------------------------------------- #
# module_type → SO role
# --------------------------------------------------------------------------- #

# Canonical mapping for every module_type emitted across the codebase.
# When an entry maps to None, the exporter falls back to SO_ENGINEERED_REGION.
MODULE_TYPE_TO_SO: Dict[str, str] = {
    # --- viral payload cassettes (boundary-defined regions) ---
    "lentiviral_payload":                SO_ENGINEERED_REGION,
    "lentiviral_cis_element":            SO_REGULATORY_REGION,
    "lentiviral_upstream_regulatory":    SO_REGULATORY_REGION,
    "lentiviral_downstream_regulatory":  SO_REGULATORY_REGION,
    "aav_payload":                       SO_ENGINEERED_REGION,
    "tdna_module":                       SO_ENGINEERED_REGION,

    # --- recombination-flanked cassettes ---
    "floxed_cassette":                   SO_ENGINEERED_REGION,
    "frt_flanked_cassette":              SO_ENGINEERED_REGION,
    "gateway_entry_cassette":            SO_ENGINEERED_REGION,
    "gateway_dest_cassette":             SO_ENGINEERED_REGION,
    "gateway_recombination":             SO_RECOMBINATION_SITE,
    "cre_loxp_recombination":            SO_RECOMBINATION_SITE,
    "flp_frt_recombination":             SO_RECOMBINATION_SITE,
    "integrase_landing_pad":             SO_RECOMBINATION_SITE,

    # --- transposons ---
    "tn7_transposon":                    SO_TRANSPOSABLE_ELEMENT,
    "sleeping_beauty_transposon":        SO_TRANSPOSABLE_ELEMENT,
    "piggybac_transposon":               SO_TRANSPOSABLE_ELEMENT,

    # --- replication ---
    "bacterial_replication":             SO_ORIGIN_OF_REPLICATION,
    "yeast_replication":                 SO_ORIGIN_OF_REPLICATION,
    "phage_replication":                 SO_ORIGIN_OF_REPLICATION,
    "bac_f_replicon":                    SO_ORIGIN_OF_REPLICATION,
    "ebv_episomal_module":               SO_ORIGIN_OF_REPLICATION,
    "sv40_replication_module":           SO_ORIGIN_OF_REPLICATION,

    # --- regulation ---
    "insulator":                         SO_INSULATOR,
    "insulated_expression_block":        SO_INSULATOR,
    "tet_regulator_cassette":            SO_REGULATORY_REGION,
    "aid_degron_system":                 SO_REGULATORY_REGION,
    "fkbp_frb_dimerization":             SO_REGULATORY_REGION,
    "counter_selection_module":          SO_CDS,

    # --- MCS / screening ---
    "ivt_cloning_cassette":              SO_MCS,
    "lac_alpha_blue_white_module":       SO_MCS,

    # --- expression cassettes ---
    "pol2_cassette":                     SO_ENGINEERED_REGION,
    "pol2_expression_cassette":          SO_ENGINEERED_REGION,
    "pol3_expression_cassette":          SO_ENGINEERED_REGION,
    "guide_expression_cassette":         SO_ENGINEERED_REGION,
    "bacterial_selection_cassette":      SO_ENGINEERED_REGION,
    "mammalian_selection_cassette":      SO_ENGINEERED_REGION,
    "selection_cassette":                SO_ENGINEERED_REGION,
    "bacterial_backbone":                SO_ENGINEERED_REGION,

    # --- Pol2 cassette submodules ---
    "upstream_regulatory_module":        SO_REGULATORY_REGION,
    "downstream_regulatory_module":      SO_REGULATORY_REGION,
    "cds_module":                        SO_CDS,

    # --- CDS submodules ---
    "protein_module":                    SO_POLYPEPTIDE_REGION,
    "protein_submodule":                 SO_POLYPEPTIDE_REGION,
    "nls_module":                        SO_NLS,
    "tag_module":                        SO_EPITOPE_TAG,
    "linker_module":                     SO_FLEXIBLE_LINKER,
    "gap_module":                        SO_REGION,

    # --- lacZα submodules ---
    "lac_alpha_cds":                     SO_CDS,
    "lac_promoter":                      SO_PROMOTER,
    "lac_operator":                      SO_OPERATOR,
    "cap_binding_site":                  SO_REGULATORY_REGION,
    "mcs":                               SO_MCS,
    "t7_phage_promoter":                 SO_PROMOTER,
    "t3_phage_promoter":                 SO_PROMOTER,
    "sp6_phage_promoter":                SO_PROMOTER,
    "m13_fwd_primer":                    SO_PRIMER_BINDING_SITE,
    "m13_rev_primer":                    SO_PRIMER_BINDING_SITE,
}

# --------------------------------------------------------------------------- #
# pLannotate feature type → SO role
# (pLannotate uses a controlled vocabulary that mostly matches SO already)
# --------------------------------------------------------------------------- #

FEATURE_TYPE_TO_SO: Dict[str, str] = {
    "promoter":                 SO_PROMOTER,
    "CDS":                      SO_CDS,
    "gene":                     SO_GENE,
    "terminator":               SO_TERMINATOR,
    "polyA_signal":             SO_POLYA_SIGNAL,
    "RBS":                      SO_RBS,
    "kozak":                    SO_KOZAK,
    "operator":                 SO_OPERATOR,
    "enhancer":                 SO_ENHANCER,
    "insulator":                SO_INSULATOR,
    "rep_origin":               SO_ORIGIN_OF_REPLICATION,
    "ori":                      SO_ORIGIN_OF_REPLICATION,
    "oriT":                     SO_ORIGIN_OF_TRANSFER,
    "LTR":                      SO_LTR,
    "ITR":                      SO_ITR,
    "protein_bind":             SO_REGULATORY_REGION,
    "misc_feature":             SO_REGION,
    "misc_recomb":              SO_RECOMBINATION_SITE,
    "regulatory":               SO_REGULATORY_REGION,
    "primer_bind":              SO_PRIMER_BINDING_SITE,
    "polylinker":               SO_MCS,
    "MCS":                      SO_MCS,
    "sig_peptide":              SO_LOCALIZATION_SIGNAL,
    "mat_peptide":              SO_POLYPEPTIDE_REGION,
    "domain":                   SO_PROTEIN_DOMAIN,
    "mobile_element":           SO_TRANSPOSABLE_ELEMENT,
    "repeat_region":            SO_DIRECT_REPEAT,
    "inverted_repeat":          SO_INVERTED_REPEAT,
    "LTR_retrotransposon":      SO_TRANSPOSABLE_ELEMENT,
    "protein_generator":        SO_CDS,
}


# --------------------------------------------------------------------------- #
# SO role → SBO interaction participant role
# --------------------------------------------------------------------------- #
# Used by the Interaction builder: given a part's SO role, what participant
# role does it most naturally play in a transcription/translation reaction?
SO_TO_SBO_PARTICIPATION: Dict[str, str] = {
    SO_PROMOTER:            SBO_STIMULATOR,
    SO_ENHANCER:            SBO_STIMULATOR,
    SO_OPERATOR:            SBO_INHIBITOR,      # if bound by a repressor
    SO_RBS:                 SBO_STIMULATOR,     # translation initiation
    SO_KOZAK:               SBO_STIMULATOR,
    SO_CDS:                 SBO_TEMPLATE,
    SO_TERMINATOR:          SBO_INHIBITOR,      # terminates transcription
    SO_POLYA_SIGNAL:        SBO_MODIFIER,
    SO_INSULATOR:           SBO_MODIFIER,
}


def so_role_for_module_type(module_type: Optional[str]) -> str:
    """Return the canonical SO role URI for a module_type string.

    Falls back to SO_ENGINEERED_REGION for unknown types so every emitted
    annotation carries a valid role.
    """
    if not module_type:
        return SO_ENGINEERED_REGION
    return MODULE_TYPE_TO_SO.get(module_type, SO_ENGINEERED_REGION)


def so_role_for_feature_type(feature_type: Optional[str]) -> str:
    """Return the canonical SO role URI for a pLannotate feature type."""
    if not feature_type:
        return SO_REGION
    return FEATURE_TYPE_TO_SO.get(feature_type, SO_REGION)


def sbo_participation_for_so(so_role: Optional[str]) -> str:
    """Return the most-likely SBO participation role for a given SO role.

    Defaults to SBO_MODIFIER when unknown.
    """
    if not so_role:
        return SBO_MODIFIER
    return SO_TO_SBO_PARTICIPATION.get(so_role, SBO_MODIFIER)


def enrich_annotation_with_roles(annotation: Dict) -> Dict:
    """Populate sbo_role / so_role fields on a hierarchical_annotation dict.

    - Uses the annotation's `module_type` first; if absent, falls back to the
      `type` field (pLannotate feature type).
    - Idempotent: if `sbo_role` is already set, returns the annotation
      unchanged.
    - Mutates in place AND returns the annotation for call-site ergonomics.
    """
    if annotation.get("sbo_role"):
        return annotation

    module_type = annotation.get("module_type")
    if module_type:
        so_role = so_role_for_module_type(module_type)
    else:
        so_role = so_role_for_feature_type(annotation.get("type"))

    annotation["sbo_role"] = so_role           # SO URI (kept on key "sbo_role" for
                                               # legacy SBOL2 compat — see sbol_io.py)
    annotation["so_role"] = so_role
    annotation["sbo_participation"] = sbo_participation_for_so(so_role)
    return annotation


__all__ = [
    "MODULE_TYPE_TO_SO",
    "FEATURE_TYPE_TO_SO",
    "SO_TO_SBO_PARTICIPATION",
    "so_role_for_module_type",
    "so_role_for_feature_type",
    "sbo_participation_for_so",
    "enrich_annotation_with_roles",
    # SO constants re-exported for callers
    "SO_PROMOTER", "SO_CDS", "SO_TERMINATOR", "SO_POLYA_SIGNAL", "SO_RBS",
    "SO_KOZAK", "SO_OPERATOR", "SO_ENHANCER", "SO_INSULATOR",
    "SO_ORIGIN_OF_REPLICATION", "SO_ORIGIN_OF_TRANSFER",
    "SO_LTR", "SO_ITR", "SO_MCS", "SO_RECOMBINATION_SITE",
    "SO_ATTB", "SO_ATTP", "SO_ATTL", "SO_ATTR",
    "SO_NLS", "SO_EPITOPE_TAG", "SO_FLEXIBLE_LINKER", "SO_RIBOSOMAL_SKIP",
    "SO_POLYPEPTIDE_REGION", "SO_PROTEIN_DOMAIN",
    "SO_ENGINEERED_REGION", "SO_REGION",
    # SBO constants re-exported
    "SBO_GENETIC_PRODUCTION", "SBO_TRANSCRIPTION", "SBO_TRANSLATION",
    "SBO_STIMULATION", "SBO_INHIBITION", "SBO_CLEAVAGE",
    "SBO_NON_COVALENT_BINDING", "SBO_CONTROL", "SBO_CONVERSION",
    "SBO_STIMULATOR", "SBO_INHIBITOR", "SBO_TEMPLATE",
    "SBO_PRODUCT", "SBO_REACTANT", "SBO_MODIFIER",
    "SBO_PROMOTER_ROLE",
]
