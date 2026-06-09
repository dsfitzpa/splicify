"""SO term -> feature_type vocabulary mapping.

GenoLIB SBOL records use SO Resource URIs (lowercase "so:0000316",
uppercase "SO:0000296", or mixed); we normalise to lowercase prefix
and look up here. Vocabulary matches pLannotate's emitted Type values
so the rest of the pipeline (rule_based_module_detector,
hierarchical_annotator, interaction_builder) doesn't need to learn
any new strings.

Sequence Ontology accessions: http://www.sequenceontology.org/
"""
from __future__ import annotations

# Map SO accession (lowercase, "so:NNNNNNN" form) -> Type string.
SO_TO_TYPE: dict[str, str] = {
    "so:0000316": "CDS",                     # CDS
    "so:0000167": "promoter",                # promoter
    "so:0000141": "terminator",              # terminator
    "so:0000296": "rep_origin",              # origin_of_replication
    "so:0000165": "enhancer",                # enhancer
    "so:0000551": "polyA_signal",            # polyA_signal_sequence
    "so:0000553": "polyA_site",              # polyA_site
    "so:0000188": "intron",                  # intron
    "so:0000147": "exon",                    # exon
    "so:0000286": "LTR",                     # long_terminal_repeat
    "so:0005850": "primer_bind",             # primer_binding_site (filtered out)
    "so:0000552": "RBS",                     # Shine_Dalgarno_sequence
    "so:0000139": "RBS",                     # ribosome_entry_site
    "so:0000298": "misc_recomb",             # recombination_feature
    "so:0000704": "gene",                    # gene
    "so:0001837": "mobile_element",          # mobile_genetic_element
    "so:0005836": "regulatory",              # regulatory_region
    "so:0000724": "oriT",                    # origin_of_transfer
    "so:0000657": "repeat_region",           # repeat_region
    "so:0000174": "promoter",                # TATA_box -> promoter family
    "so:0000204": "5'UTR",                   # five_prime_UTR
    "so:0000205": "3'UTR",                   # three_prime_UTR
    "so:0000253": "tRNA",                    # tRNA
    "so:0000233": "mRNA",                    # mature_transcript
    "so:0000178": "operon",                  # operon
    "so:0000410": "mRNA",                    # mRNA
    "so:0000252": "rRNA",                    # rRNA
    "so:0000655": "ncRNA",                   # ncRNA
    "so:0001877": "ncRNA",                   # lnc_RNA
    "so:0000673": "transcript",              # transcript
    "so:0000275": "snoRNA",                  # snoRNA
    "so:0000234": "mRNA",                    # mRNA_with_5'utr_3'utr
    # SO:0000001 is the generic "region" — too vague; treated as misc_feature
    # at emit time when name doesn't refine it.
    "so:0000001": "misc_feature",
}


def so_to_type(so_resource: str, default: str = "misc_feature") -> str:
    """Map a SO Resource URI fragment ("so:0000316", "SO:0000316", or
    full URI) to a Type string. Returns ``default`` for unknown terms.
    """
    if not so_resource:
        return default
    key = so_resource.strip().lower()
    # Strip prefixes if a full URI was passed
    if "/" in key:
        key = key.rsplit("/", 1)[-1]
    if ":" not in key:
        key = f"so:{key}"
    return SO_TO_TYPE.get(key, default)
