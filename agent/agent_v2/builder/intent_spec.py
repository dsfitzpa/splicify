"""IntentSpec: the parsed user requirements that the builder
satisfies and the verifier checks against. Built either by a small
LLM call or from a deterministic role taxonomy when the user names
parts explicitly."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


# Canonical roles the builder + verifier understand. Maps roughly to
# the rule_based_module_detector's submodule_types and feature_class
# values, so the verifier can compare apples-to-apples.
KNOWN_ROLES = {
    "promoter",
    "cds",
    "polya",
    "selection_marker",
    "origin",
    "ltr",
    "lentiviral_cis",
    "wpre",
    "scaffold",          # gRNA scaffold etc.
    "spacer",            # sgRNA spacer
    "stuffer",
    "tag",               # FLAG, NLS, etc.
    "kozak",
    "polylinker",
    "att_site",
    "loxP",
    "frt",
    "itr",
    "enhancer",
}


# Module types the verifier knows how to score (subset of the
# rule_based_module_detector's emitted module_types).
KNOWN_MODULES = {
    "mammalian_pol2_expression_cassette",
    "guide_expression_cassette",
    "bacterial_selection_cassette",
    "mammalian_selection_cassette",
    "lentiviral_payload",
    "lac_alpha_blue_white_module",
    "gateway_destination_cassette",
    "tet_inducible_expression_cassette",
}


# Function tags surface higher-level intent that doesn't map cleanly
# to a single module — verifier maps each to a checklist of expected
# modules + interactions.
KNOWN_FUNCTIONS = {
    "expression",                # promoter → CDS → polyA
    "cloning_vector",            # backbone + selection + MCS or GG-acceptor
    "lentiviral_vector",
    "crispr_knockout",           # sgRNA cassette + Cas9 cassette
    "prime_editing",             # pegRNA cassette + prime editor
    "sgrna_cloning",             # POL3-GG-01 backbone
}


@dataclass
class IntentSpec:
    """What the user wants. The builder treats every required entry
    as a hard constraint; preferred entries are nice-to-have."""
    function: str                                       # one of KNOWN_FUNCTIONS, or "custom"
    host_scope: str = "mammalian"                       # "bacterial" | "mammalian" | "yeast" | …
    required_modules: list[str] = field(default_factory=list)
    required_interactions: list[str] = field(default_factory=list)
    required_roles: list[str] = field(default_factory=list)  # e.g. ["promoter", "cds", "polya"]
    preferred_features: list[str] = field(default_factory=list)  # e.g. ["CAG promoter", "Cas9"]
    forbidden_features: list[str] = field(default_factory=list)
    topology: str = "circular"

    @classmethod
    def for_expression_cassette(cls, host_scope: str = "mammalian",
                                  preferred: list[str] | None = None) -> "IntentSpec":
        """Convenience: build me a plasmid that expresses one CDS in
        the host_scope. Maps to promoter + CDS + polyA + selection +
        origin."""
        return cls(
            function="expression",
            host_scope=host_scope,
            required_modules=["mammalian_pol2_expression_cassette" if host_scope == "mammalian"
                              else "bacterial_expression_cassette"],
            required_interactions=["expression_cassette"],
            required_roles=["promoter", "cds", "polya", "selection_marker", "origin"],
            preferred_features=preferred or [],
        )

    @classmethod
    def for_sgrna_cloning_backbone(cls) -> "IntentSpec":
        return cls(
            function="sgrna_cloning",
            host_scope="mammalian",
            required_modules=["guide_expression_cassette"],
            required_roles=["promoter", "scaffold", "stuffer", "selection_marker", "origin"],
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "function": self.function,
            "host_scope": self.host_scope,
            "required_modules": list(self.required_modules),
            "required_interactions": list(self.required_interactions),
            "required_roles": list(self.required_roles),
            "preferred_features": list(self.preferred_features),
            "forbidden_features": list(self.forbidden_features),
            "topology": self.topology,
        }
