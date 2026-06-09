"use client";

import React, { useState, useMemo, useCallback, useEffect, useRef } from "react";
import dynamic from "next/dynamic";
import { GenBankDownloadButton } from "./components_GenBankDownloadButton";
import CircularPlasmidViewer, { PlasmidInteraction } from "./CircularPlasmidViewer";
import LinearSequenceViewer from "./LinearSequenceViewer";
import FlatPlasmidViewer from "./FlatPlasmidViewer";

const SeqVizAny = dynamic(async () => {
  const mod: any = await import("seqviz");
  return (mod.default ?? mod.SeqViz ?? mod) as React.ComponentType<any>;
}, { ssr: false });

export type SeqVizAnnotation = {
  name: string;
  start: number; // 0-indexed internally
  end: number;   // 0-indexed internally
  direction?: 1 | -1 | 0;
  strand?: 1 | -1 | 0;
  length?: number;
  color?: string;
  type?: string;        // Feature type (CDS, promoter, etc.)
  description?: string; // Feature description/note
  sseqid?: string;      // Subject sequence ID from BLAST
  db?: string;          // Database source: "swissprot", "snapgene", etc.
  kb_data?: {           // Knowledge base metadata for gene card
    // Common fields
    source_type?: "swissprot" | "feature_kb";
    
    // SwissProt fields
    protein_name?: string;
    gene_name?: string;
    organism?: string;
    taxonomy_id?: string;
    protein_existence?: string;
    entry_name?: string;
    
    // Feature KB fields
    feature_id?: string;
    feature_name?: string;
    feature_type?: string;
    feature_class?: string;      // promoter, enhancer, cds_payload, terminator, etc.
    subclass?: string;
    host_scope?: string[];       // ["mammalian", "bacterial", etc.]
    delivery_scope?: string[];   // ["plasmid", "viral", etc.]
    descriptions?: string[];     // Array of description strings
    annotation_source?: string;  // "SnapGene", "FPbase", etc.
    polymerase_class?: string;   // "pol_ii", "pol_iii", etc.
    orientation_requirements?: string;
    frame_semantics?: string;
  } | null;
  layer?: "feature" | "module" | "motif" | "gap" | "cloning_feature" | "translation";  // Annotation layer type
  added_by_design?: boolean;  // True when the feature was injected by the CRISPR/cloning design pipeline (e.g. pegRNAs, primers).
  module_type?: string;      // Module type (for module layer)
  module_family?: string;    // Module family
  motif_type?: string;       // Motif type (start_codon, stop_codon, etc.)
  payload_id?: string;      // Payload ID (for module layer)
  source?: string;           // Annotation source (motif_detector, etc.)
  sequence?: string;         // Sequence at this annotation position
  metadata?: Record<string, any>;  // Module metadata for heuristic information
  // Cloning-feature-layer metadata (populated when layer === "cloning_feature")
  feature_family?: "restriction_site_II" | "restriction_site_IIs" | "gateway_att" | "primer_design_warning";
  subtype?: string;
  cut_profile?: {
    cut_top: number;
    cut_bottom: number;
    overhang_seq: string;
    overhang_type: "5prime" | "3prime" | "blunt";
    overhang_len: number;
  } | null;
};

// Cloning-feature response payload from the backend scanner.
export type CloningFeaturesPayload = {
  features: any[];
  cut_count_per_enzyme: Record<string, number>;
  non_cutters: string[];
  enabled_sets: string[];
};

export type CutterCountFilter = "none" | "unique" | "2" | "3" | "all";

type AnnotationFormData = {
  name: string;
  start: string; // 1-indexed for display
  end: string;   // 1-indexed for display
  direction: "forward" | "reverse" | "none";
  type?: string;        // Feature type (read-only)
  description?: string; // Feature description (read-only)
};

type AnnotationFormErrors = {
  name?: string;
  start?: string;
  end?: string;
};

type Selection = {
  start: number;
  end: number;
  clockwise: boolean;
  ref?: string;
};

export type RestrictionSiteAnnotation = {
  name: string;
  start: number;
  end: number;
  direction: 1 | -1;
  re_type: "type2_re" | "type2s_re" | "gateway" | "cre_lox" | "flp_frt";
  recognition_seq: string;
  color: string;
};

type ReSiteFilter = "none" | "type2s_re" | "type2_re" | "gateway" | "cre_lox" | "all";

type Props = {
  initialSequence: string;
  initialAnnotations: SeqVizAnnotation[];
  title?: string;
  plasmidName?: string;  // Original uploaded filename (without extension)
  viewerMode?: "both" | "linear" | "circular";
  /** Initial state of the left pane when viewerMode="both". Defaults to "circular"; pass "flat" for linear genomic slices.
   *  Caller responsibility: set this once on first mount; subsequent updates do NOT re-sync the toggle. */
  initialLeftPaneView?: "circular" | "flat";
  height?: number;
  circular?: boolean;
  plannotateEndpoint?: string;
  restrictionSites?: RestrictionSiteAnnotation[];
  analyzeIntentEndpoint?: string;
  onAnalysisComplete?: (analysis: string, moduleGraph: any) => void;  // Callback for chat integration
  autoAnnotateOnMount?: boolean;  // Fire handleLLMAnnotate once when sequence is loaded
  // Pre-loaded annotation extras from an out-of-band annotate call
  // (e.g. /agent_v2/annotate-on-upload). When supplied, the viewer seeds
  // its interactions + cloning-features state from these props, so the
  // "Show interaction chords" toggle is enabled without re-running
  // plannotate inside the component.
  initialInteractions?: PlasmidInteraction[];
  initialCloningFeatures?: CloningFeaturesPayload | null;
  // Region the parent wants the viewer to focus on (e.g. coordinates
  // pulled from interpreter-agent citations). The viewer treats the
  // midpoint as the centre-on position for the linear viewer and as
  // the focused-arc target for the circular viewer.
  focusRegion?: { start: number; end: number } | null;
};

// CSS to apply Source Sans Pro only to labels, sequence text uses Hack font
const viewerStyles = `
  @import url('https://cdn.jsdelivr.net/npm/hack-font@3/build/web/hack.css');

  .seqviz-viewer-container {
    font-family: "Source Sans Pro", sans-serif;
  }
  /* Sequence text uses Hack font with specified color */
  .seqviz-viewer-container .la-vz-linear-scroller text tspan,
  .seqviz-viewer-container [class*="sequence"] text,
  .seqviz-viewer-container [class*="Sequence"] text,
  .seqviz-viewer-container .la-vz-seqs text,
  .seqviz-viewer-container .la-vz-seq text,
  .seqviz-viewer-container text[class*="seq"],
  .seqviz-viewer-container g[class*="seq"] text,
  .seqviz-viewer-container tspan {
    font-family: "Hack", monospace !important;
    fill: #72777e !important;
  }

  /* ---------- Plasmid-viewer toolbar (nav-bar style: flat, no bubbles) ---------- */
  .splicify-toolbar {
    background: transparent;
    border: none;
    border-bottom: 1px solid rgba(219, 239, 231, 0.10);
    border-radius: 0;
    box-shadow: none;
    padding: 4px 0 14px;
    display: flex;
    align-items: center;
    justify-content: center;
    flex-wrap: wrap;
    gap: 8px;
    font-family: var(--font-body);
  }
  .splicify-tb-group {
    display: flex;
    align-items: center;
    gap: 18px;
    flex-wrap: wrap;
  }
  .splicify-tb-sep {
    width: 1px;
    align-self: stretch;
    background-color: rgba(219, 239, 231, 0.12);
    margin: 4px 10px;
  }
  .splicify-tb-btn {
    background: transparent;
    color: rgba(219, 239, 231, 0.70);
    border: none;
    border-bottom: 1.5px solid transparent;
    border-radius: 0;
    padding: 4px 2px 6px;
    font-size: 13.5px;
    font-weight: 400;
    line-height: 1.2;
    cursor: pointer;
    transition: color 140ms ease, border-color 140ms ease;
    white-space: nowrap;
    font-family: var(--font-body);
  }
  .splicify-tb-btn:hover:not(:disabled) {
    color: var(--brass-400);
  }
  .splicify-tb-btn:active:not(:disabled) {
    color: var(--brass-500);
  }
  .splicify-tb-btn:disabled {
    opacity: 0.4;
    cursor: not-allowed;
  }
  .splicify-tb-btn.is-purple { color: #d8a3e8; }
  .splicify-tb-btn.is-purple:hover:not(:disabled) { color: #f3e7f7; }
  .splicify-tb-btn.is-magenta { color: #f5a8c8; }
  .splicify-tb-btn.is-magenta:hover:not(:disabled) { color: #ffd0e3; }
  .splicify-tb-btn.is-outline {
    color: var(--brass-400);
  }
  .splicify-tb-btn.is-light {
    color: var(--mint-200);
    font-weight: 600;
  }
  .splicify-tb-btn.is-active {
    color: var(--brass-400);
    font-weight: 600;
    border-bottom-color: var(--brass-400);
  }
  .splicify-tb-btn.is-used {
    color: rgba(219, 239, 231, 0.35);
  }
  .splicify-tb-btn.is-used:hover:not(:disabled) {
    color: rgba(219, 239, 231, 0.55);
  }
  .splicify-tb-toggle {
    display: inline-flex;
    align-items: stretch;
    background: transparent;
    border-radius: 0;
    overflow: visible;
    border: none;
    gap: 14px;
    margin-left: 4px;
  }
  .splicify-tb-toggle button {
    background: transparent;
    color: rgba(219, 239, 231, 0.55);
    border: none;
    border-bottom: 1.5px solid transparent;
    cursor: pointer;
    padding: 4px 2px 6px;
    font-size: 13.5px;
    font-weight: 400;
    font-family: var(--font-body);
    transition: color 140ms ease, border-color 140ms ease;
  }
  .splicify-tb-toggle button:hover {
    color: var(--brass-400);
  }
  .splicify-tb-toggle button.is-on {
    background: transparent;
    color: var(--brass-400);
    font-weight: 600;
    border-bottom-color: var(--brass-400);
  }
`;

export default function InteractiveSequenceViewer({
  initialSequence,
  initialAnnotations,
  title = "Visualization preview",
  plasmidName,
  viewerMode = "both",
  initialLeftPaneView = "circular",
  height = 500,
  circular = true,
  plannotateEndpoint,
  restrictionSites,
  analyzeIntentEndpoint,
  onAnalysisComplete,
  autoAnnotateOnMount = false,
  initialInteractions,
  initialCloningFeatures,
  focusRegion,
}: Props) {
  // Core state
  const [sequence, setSequence] = useState(initialSequence);
  const [annotations, setAnnotations] = useState<SeqVizAnnotation[]>(initialAnnotations);

  // Undo / redo history. Each entry is a snapshot of {seq, anns}; the active
  // state is intentionally NOT in the stack — it is what setState currently
  // holds. lastSnapshotRef tracks the previous committed state so the watcher
  // effect knows what to push onto the past stack on the next change.
  const historyPastRef = useRef<{ seq: string; anns: SeqVizAnnotation[] }[]>([]);
  const historyFutureRef = useRef<{ seq: string; anns: SeqVizAnnotation[] }[]>([]);
  const lastSnapshotRef = useRef<{ seq: string; anns: SeqVizAnnotation[] } | null>(null);
  const skipNextSnapshotRef = useRef(false);

  // Find / Cmd+F panel state
  const [findOpen, setFindOpen] = useState(false);
  const [findQuery, setFindQuery] = useState("");
  const [findIdx, setFindIdx] = useState(0);
  const findInputRef = useRef<HTMLInputElement>(null);

  // Functional interactions (from annotate_sequence_llm / rule-based detector)
  const [interactions, setInteractions] = useState<PlasmidInteraction[]>(initialInteractions || []);
  const [showInteractionChords, setShowInteractionChords] = useState(true);
  const [highlightedInteractionId, setHighlightedInteractionId] = useState<string | null>(null);
  const [interactionDescription, setInteractionDescription] = useState<string>("");
  const [interactionDescLoading, setInteractionDescLoading] = useState(false);


  // Selection state
  const [selection, setSelection] = useState<Selection | null>(null);
  const selectionRef = useRef<Selection | null>(null); // Sync ref for double-click handler
  const sequenceRef = useRef<string>("");
  const [centerOnPosition, setCenterOnPosition] = useState<number | null>(null);
  const [topOnPosition, setTopOnPosition] = useState<number | null>(null);

  // Parent-driven focus: when `focusRegion` changes, sync the selection
  // (so the highlighted band appears on both viewers) and centre the
  // linear strip on the midpoint. The circular viewer uses `selection`
  // for its highlight arc, so setting it gives both views a visual cue.
  useEffect(() => {
    if (!focusRegion) return;
    const { start, end } = focusRegion;
    if (typeof start !== "number" || typeof end !== "number") return;
    const lo = Math.min(start, end);
    const hi = Math.max(start, end);
    if (hi <= lo) return;
    setSelection({ start: lo, end: hi, clockwise: true });
    setCenterOnPosition(Math.round((lo + hi) / 2));
  }, [focusRegion?.start, focusRegion?.end]);

  // Modal state
  const [modalOpen, setModalOpen] = useState(false);
  const [modalMode, setModalMode] = useState<"add" | "edit">("add");
  const [editingIndex, setEditingIndex] = useState<number | null>(null);
  const [activeTab, setActiveTab] = useState<"info" | "edit">("info");

  // Translation-annotation creation: small dialog asks the user for the
  // reading-frame direction before computing the AA strip on the selected
  // range. State is held here so the dialog can read selection at the
  // moment the user opens it.
  const [translationModalOpen, setTranslationModalOpen] = useState(false);
  const [translationDirection, setTranslationDirection] = useState<"forward" | "reverse">("forward");

  // Sequence editor state
  const [sequenceEditorOpen, setSequenceEditorOpen] = useState(false);
  const [editedSequence, setEditedSequence] = useState("");
  const [sequenceError, setSequenceError] = useState("");

  // Add/Change sequence modal state
  const [addChangeModalOpen, setAddChangeModalOpen] = useState(false);
  const [newSequenceText, setNewSequenceText] = useState("");
  const [addChangeError, setAddChangeError] = useState("");

  // Import-Annotations modal state
  const [importModalOpen, setImportModalOpen] = useState(false);
  const [importLoading, setImportLoading] = useState(false);
  const [importError, setImportError] = useState("");
  const [importMaxMismatches, setImportMaxMismatches] = useState(0);
  const [importFileName, setImportFileName] = useState<string>("");
  const [importResult, setImportResult] = useState<{
    summary: { n_input: number; n_matched: number; n_annotations: number; n_unmatched: number };
    unmatched: Array<{ name: string; sequence: string; reason: string }>;
  } | null>(null);
  const importFileInputRef = useRef<HTMLInputElement>(null);

  // Design-Primers modal state
  const [primerModalOpen, setPrimerModalOpen] = useState(false);
  const [primerLoading, setPrimerLoading] = useState(false);
  const [primerError, setPrimerError] = useState("");
  // Guide-design modal state
  const [guideModalOpen, setGuideModalOpen] = useState(false);
  const [guideLoading, setGuideLoading] = useState(false);
  const [guideError, setGuideError] = useState("");
  const [guideForm, setGuideForm] = useState({
    region_start: "",
    region_end: "",
    pam: "NGG",
    guide_length: "20",
    pam_position: "3prime" as "3prime" | "5prime",
    max_guides: "20",
    min_score: "0",
    score_method: "doench2014" as "doench2014" | "heuristic",
  });
  const [guideResult, setGuideResult] = useState<null | {
    summary: { n_candidates: number; n_returned: number; score_method: string;
               score_method_requested: string; doench2014_eligible: boolean };
    n_added: number;
    guides: any[];
    design_region_1based?: string;
  }>(null);

  // --- Design pegRNA (prime editing) ---
  const [pegrnaModalOpen, setPegrnaModalOpen] = useState(false);
  const [pegrnaLoading, setPegrnaLoading] = useState(false);
  const [pegrnaError, setPegrnaError] = useState("");
  const [pegrnaForm, setPegrnaForm] = useState({
    edit_start: "",
    edit_end: "",
    alt: "",
    edit_type: "substitution" as "substitution" | "insertion" | "deletion",
    n_results: "3",
    use_pe3: true,
  });
  const [pegrnaResult, setPegrnaResult] = useState<null | {
    summary: { n_sgrnas_scanned: number; n_valid_pegRNAs: number; n_candidates: number;
               n_returned: number; use_pe3: boolean; edit_type: string;
               edit_ref: string; edit_alt: string; model: string };
    n_added: number;
  }>(null);

  // Left-pane view toggle (Circular <-> Flat). Right pane is always Linear.
  // Initial value comes from the prop so callers can default genomic slices to flat
  // without the user having to click the toggle.
  const [leftPaneView, setLeftPaneView] = useState<"circular" | "flat">(initialLeftPaneView);

  // Viewer layout: split shows both panes (default); circular shows only the
  // left pane (circular/flat per leftPaneView); linear shows only the linear
  // pane. Full-screen modes expand the visible pane to square proportions while
  // preserving the container's horizontal length.
  const [viewerLayout, setViewerLayout] = useState<"split" | "circular" | "linear">("split");
  const viewerOuterRef = useRef<HTMLDivElement>(null);
  const [viewerOuterHeight, setViewerOuterHeight] = useState<number>(height);
  useEffect(() => {
    const el = viewerOuterRef.current;
    if (!el) return;
    const update = () => setViewerOuterHeight(el.offsetHeight);
    update();
    const ro = new ResizeObserver(update);
    ro.observe(el);
    return () => ro.disconnect();
  }, [viewerLayout, height]);

  const [primerForm, setPrimerForm] = useState({
    application: "fragment" as "fragment" | "sanger" | "illumina",
    region_start: "",     // 1-indexed for display
    region_end: "",
    excluded_start: "",
    excluded_end: "",
    product_size_min: "100",
    product_size_max: "300",
    primer_min_tm: "",
    primer_opt_tm: "",
    primer_max_tm: "",
    primer_min_size: "",
    primer_opt_size: "",
    primer_max_size: "",
    num_return: "5",
    pair_label: "",
    primer_fwd_name: "", // overrides default per-app suffix if set
    primer_rev_name: "", // overrides default per-app suffix if set
  });

  // Default forward / reverse name suffixes per application. Used to seed
  // the auto-generated primer names; the user can still override via the
  // primer_fwd_name / primer_rev_name fields.
  const PRIMER_NAME_SUFFIX: Record<string, { fwd: string; rev: string }> = {
    fragment: { fwd: "_F",  rev: "_R"  },
    sanger:   { fwd: "_SF", rev: "_SR" },
    illumina: { fwd: "_iF", rev: "_iR" },
  };
  const ILLUMINA_MAX_AMPLICON_BP = 600;
  const [primerResult, setPrimerResult] = useState<any>(null);

  // Annotation form state
  const [formData, setFormData] = useState<AnnotationFormData>({
    name: "",
    start: "",
    end: "",
    direction: "forward",
  });
  const [formErrors, setFormErrors] = useState<AnnotationFormErrors>({});

  // Annotation state
  const [plannotateLoading, setPlannotateLoading] = useState(false);
  const [plannotateError, setPlannotateError] = useState("");
  // Toggle for new annotation pipeline
  const [heuristicLoading, setHeuristicLoading] = useState(false);
  const [llmLoading, setLlmLoading] = useState(false);
  const [cloningScanLoading, setCloningScanLoading] = useState(false);

  // Analyze intent state
  const [analyzeIntentLoading, setAnalyzeIntentLoading] = useState(false);
  const [analyzeIntentError, setAnalyzeIntentError] = useState("");
  const [intentAnalysis, setIntentAnalysis] = useState<string | null>(null);
  const [moduleGraphData, setModuleGraphData] = useState<any>(null);  // For download
  const [uploadedFileName, setUploadedFileName] = useState<string | null>(null);  // Track filename from direct upload

  // Effective plasmid name: uploaded filename takes precedence over prop
  const effectivePlasmidName = uploadedFileName || plasmidName;

  // RE site toggle state (legacy, for externally-passed restrictionSites)
  const [reSiteFilter, setReSiteFilter] = useState<ReSiteFilter>("none");
  const [showReSiteMenu, setShowReSiteMenu] = useState(false);

  // Cloning-feature toggle state (backend-scanner driven)
  const [cloningFeatures, setCloningFeatures] = useState<CloningFeaturesPayload | null>(initialCloningFeatures ?? null);
  // Default OFF: cloning features (restriction sites, Gateway att, PCR
  // warnings) stay hidden until the user clicks "Scan Cloning Features".
  // Restriction-site density can dominate the viewer on large plasmids;
  // user opts in when they actually want to plan a clone.
  const [showCloningFeatures, setShowCloningFeatures] = useState(false);
  const [showCloningMenu, setShowCloningMenu] = useState(false);
  const [showCloningReII, setShowCloningReII] = useState(true);
  const [showCloningReIIs, setShowCloningReIIs] = useState(true);
  const [showCloningGateway, setShowCloningGateway] = useState(true);
  const [showCloningPcr, setShowCloningPcr] = useState(true);
  const [cutterFilter, setCutterFilter] = useState<CutterCountFilter>("all");
  const [showNonCutters, setShowNonCutters] = useState(true);

  // Module + interaction visibility (controlled from the Annotations dropdown
  // in the toolbar). Modules are filtered out of the rendered annotation set
  // when off; interaction chords reuse the existing showInteractionChords flag.
  const [showModules, setShowModules] = useState(true);
  const [showAnnotationsMenu, setShowAnnotationsMenu] = useState(false);

  // Modular circular viewer is the only viewer now (SeqViz fallback removed).
  const useCustomViewer = true;

  // Container ref for focus management
  const containerRef = useRef<HTMLDivElement>(null);
  const viewerRef = useRef<HTMLDivElement>(null);
  const fileInputRef = useRef<HTMLInputElement>(null);

  // Parse GenBank file content
  const parseGenBank = useCallback((content: string): { sequence: string; annotations: SeqVizAnnotation[] } | null => {
    try {
      const lines = content.split(/\r?\n/);
      let inOrigin = false;
      let inFeatures = false;
      let sequenceParts: string[] = [];
      const annotations: SeqVizAnnotation[] = [];
      let currentFeature: { type: string; location: string; qualifiers: Record<string, string> } | null = null;
      let currentQualifierKey = "";
      let currentQualifierValue = "";

      for (const line of lines) {
        // Check for ORIGIN section (sequence data)
        if (line.startsWith("ORIGIN")) {
          inOrigin = true;
          inFeatures = false;
          // Save last feature if any
          if (currentFeature) {
            const ann = parseFeatureToAnnotation(currentFeature, 0); // seqLen unknown yet
            if (ann) annotations.push(ann);
            currentFeature = null;
          }
          continue;
        }

        // Check for end of record
        if (line.startsWith("//")) {
          inOrigin = false;
          inFeatures = false;
          break;
        }

        // Check for FEATURES section
        if (line.startsWith("FEATURES")) {
          inFeatures = true;
          continue;
        }

        // Parse sequence in ORIGIN section
        if (inOrigin) {
          // Sequence lines start with a number, followed by sequence data
          const seqMatch = line.match(/^\s*\d+\s+(.+)$/);
          if (seqMatch) {
            // Remove spaces from sequence
            sequenceParts.push(seqMatch[1].replace(/\s/g, ""));
          }
          continue;
        }

        // Parse features
        if (inFeatures) {
          // New feature starts with 5 spaces + feature key + location
          const featureMatch = line.match(/^     (\S+)\s+(.+)$/);
          if (featureMatch && !line.startsWith("                     /")) {
            // Save previous feature
            if (currentFeature) {
              if (currentQualifierKey && currentQualifierValue) {
                currentFeature.qualifiers[currentQualifierKey] = currentQualifierValue.replace(/^"|"$/g, "");
              }
              const ann = parseFeatureToAnnotation(currentFeature, 0);
              if (ann) annotations.push(ann);
            }
            currentFeature = {
              type: featureMatch[1],
              location: featureMatch[2],
              qualifiers: {},
            };
            currentQualifierKey = "";
            currentQualifierValue = "";
            continue;
          }

          // Qualifier line starts with 21 spaces + /
          const qualifierMatch = line.match(/^                     \/(\w+)(?:=(.*))?$/);
          if (qualifierMatch && currentFeature) {
            // Save previous qualifier
            if (currentQualifierKey && currentQualifierValue) {
              currentFeature.qualifiers[currentQualifierKey] = currentQualifierValue.replace(/^"|"$/g, "");
            }
            currentQualifierKey = qualifierMatch[1];
            currentQualifierValue = qualifierMatch[2] || "";
            continue;
          }

          // Continuation of qualifier value (21 spaces, no /)
          if (line.startsWith("                     ") && !line.startsWith("                     /") && currentQualifierKey) {
            currentQualifierValue += line.trim();
            continue;
          }

          // Location continuation
          if (currentFeature && !currentQualifierKey && line.match(/^\s+[^\/]/)) {
            currentFeature.location += line.trim();
          }
        }
      }

      // Save last feature
      if (currentFeature) {
        if (currentQualifierKey && currentQualifierValue) {
          currentFeature.qualifiers[currentQualifierKey] = currentQualifierValue.replace(/^"|"$/g, "");
        }
        const ann = parseFeatureToAnnotation(currentFeature, 0);
        if (ann) annotations.push(ann);
      }

      const sequence = sequenceParts.join("").toUpperCase();

      // Now update annotations with correct seqLen
      const finalAnnotations = annotations.map(ann => ({
        ...ann,
        length: ann.start <= ann.end ? ann.end - ann.start : (sequence.length - ann.start) + ann.end,
      }));

      // Filter out source features
      const filteredAnnotations = finalAnnotations.filter(a =>
        a.name !== "synthetic construct" &&
        a.name !== "source" &&
        !a.name.toLowerCase().includes("mol_type")
      );

      return { sequence, annotations: filteredAnnotations };
    } catch (error) {
      console.error("Error parsing GenBank file:", error);
      return null;
    }
  }, []);

  // Helper to parse a feature into an annotation
  const parseFeatureToAnnotation = (
    feature: { type: string; location: string; qualifiers: Record<string, string> },
    seqLen: number
  ): SeqVizAnnotation | null => {
    if (feature.type === "source") return null;

    // Parse location - handle complement and join
    let location = feature.location;
    let strand: 1 | -1 | 0 = 1;

    // Check for complement
    if (location.startsWith("complement(")) {
      strand = -1;
      location = location.replace(/^complement\(/, "").replace(/\)$/, "");
    }

    // Parse simple range or join
    let start = 0;
    let end = 0;

    if (location.startsWith("join(")) {
      // Handle join - take first and last positions for circular spanning
      const joinContent = location.replace(/^join\(/, "").replace(/\)$/, "");
      const parts = joinContent.split(",");
      if (parts.length >= 2) {
        const firstPart = parts[0].match(/(\d+)\.\.(\d+)/);
        const lastPart = parts[parts.length - 1].match(/(\d+)\.\.(\d+)/);
        if (firstPart && lastPart) {
          start = parseInt(firstPart[1], 10) - 1; // Convert to 0-indexed
          end = parseInt(lastPart[2], 10);
        }
      }
    } else {
      // Simple range
      const rangeMatch = location.match(/(\d+)\.\.(\d+)/);
      if (rangeMatch) {
        start = parseInt(rangeMatch[1], 10) - 1; // Convert to 0-indexed
        end = parseInt(rangeMatch[2], 10);
      } else {
        // Single position
        const singleMatch = location.match(/(\d+)/);
        if (singleMatch) {
          start = parseInt(singleMatch[1], 10) - 1;
          end = start + 1;
        }
      }
    }

    // Get name from label or note qualifier
    const name = feature.qualifiers.label ||
                 feature.qualifiers.gene ||
                 feature.qualifiers.product ||
                 feature.qualifiers.note ||
                 feature.type;

    // Get description - prefer note, but use product or other descriptive qualifiers
    const description = feature.qualifiers.note ||
                        feature.qualifiers.product ||
                        feature.qualifiers.function ||
                        "";

    const addedByDesign =
      (feature.qualifiers.added_by || "").toLowerCase().includes("crispr") ||
      (feature.qualifiers.added_by || "").toLowerCase().includes("design");

    return {
      name,
      start,
      end,
      direction: strand,
      strand: strand,
      type: feature.type,
      description: description,
      added_by_design: addedByDesign || undefined,
    };
  };

  // Handle file upload
  const handleFileUpload = useCallback((event: React.ChangeEvent<HTMLInputElement>) => {
    const file = event.target.files?.[0];
    if (!file) return;

    // Extract filename without extension for use as plasmid name
    const fileName = file.name;
    const nameWithoutExt = fileName.replace(/\.(gb|gbk|genbank)$/i, "");
    setUploadedFileName(nameWithoutExt);

    const reader = new FileReader();
    reader.onload = (e) => {
      const content = e.target?.result as string;
      if (!content) return;

      const result = parseGenBank(content);
      if (result) {
        setSequence(result.sequence);
        setAnnotations(result.annotations);
        setSelection(null);
        // Plasmid-specific state from the previous upload — interaction
        // chords, module graph, cached cloning-feature scan, etc. — must
        // be cleared so they do not bleed onto the new plasmid.
        setInteractions([]);
        setHighlightedInteractionId(null);
        setInteractionDescription("");
        setModuleGraphData(null);
        setIntentAnalysis(null);
        setCloningFeatures(null);
      } else {
        alert("Failed to parse GenBank file. Please check the file format.");
      }
    };
    reader.readAsText(file);

    // Reset input so same file can be uploaded again
    event.target.value = "";
  }, [parseGenBank]);

  // Trigger file upload dialog
  const handleUploadClick = useCallback(() => {
    fileInputRef.current?.click();
  }, []);

  // Process annotations for display (remove arrows from overlaps, rename overlap format)
  const processedAnnotations = useMemo(() => {
    return annotations.map((a, idx) => {
      let name = a.name ?? "";
      const nameLower = name.toLowerCase();
      const ann = { ...a, id: `ann-${idx}` };

      // Convert old overlap format "X→Y overlap" to new format "overlap X-Y"
      const arrowOverlapMatch = name.match(/^(.+?)→(.+?)\s+overlap\s*(.*)$/i);
      if (arrowOverlapMatch) {
        const [, fromName, toName, suffix] = arrowOverlapMatch;
        name = `overlap ${fromName}-${toName}${suffix ? ' ' + suffix : ''}`;
        ann.name = name;
      }

      if (nameLower.includes("overlap") || nameLower.includes("homology")) {
        return { ...ann, direction: 0 as const, strand: 0 as const };
      }
      return ann;
    });
  }, [annotations]);

  // Filtered RE site annotations based on active filter
  const visibleReSites = useMemo(() => {
    if (!restrictionSites?.length || reSiteFilter === "none") return [];
    if (reSiteFilter === "all") return restrictionSites;
    return restrictionSites.filter((s) => s.re_type === reSiteFilter);
  }, [restrictionSites, reSiteFilter]);

  // Cutter-count filter for backend-scanner restriction features.
  // "unique" = enzymes with exactly 1 site in the template.
  // "2" / "3" = enzymes with exactly 2 / 3 sites.
  // "all" = every recognition hit regardless of site count.
  // "none" = show only the non-cutters list (no recognition hits rendered).
  const enzymeMatchesCutterFilter = useCallback((enzymeName: string): boolean => {
    if (!cloningFeatures) return false;
    const count = cloningFeatures.cut_count_per_enzyme[enzymeName] ?? 0;
    switch (cutterFilter) {
      case "unique": return count === 1;
      case "2": return count === 2;
      case "3": return count === 3;
      case "all": return count >= 1;
      case "none": return false;
      default: return false;
    }
  }, [cloningFeatures, cutterFilter]);

  // Partition processedAnnotations into cloning-feature and non-cloning-feature sets,
  // then apply per-category and cutter-count filters. Modules are dropped from
  // the non-cloning bucket when the Annotations dropdown has them disabled.
  const { nonCloningAnns, cloningFeatureAnns } = useMemo(() => {
    const nonCloning: typeof processedAnnotations = [];
    const cloning: typeof processedAnnotations = [];
    for (const a of processedAnnotations) {
      if ((a as any).layer === "cloning_feature") cloning.push(a);
      else if (!showModules && (a as any).layer === "module") continue;
      else nonCloning.push(a);
    }
    return { nonCloningAnns: nonCloning, cloningFeatureAnns: cloning };
  }, [processedAnnotations, showModules]);

  const visibleCloningFeatures = useMemo(() => {
    if (!showCloningFeatures) return [];
    return cloningFeatureAnns.filter((a) => {
      const fam = (a as any).feature_family;
      if (fam === "restriction_site_II") {
        if (!showCloningReII) return false;
        return enzymeMatchesCutterFilter(a.name);
      }
      if (fam === "restriction_site_IIs") {
        if (!showCloningReIIs) return false;
        return enzymeMatchesCutterFilter(a.name);
      }
      if (fam === "gateway_att") return showCloningGateway;
      if (fam === "primer_design_warning") return showCloningPcr;
      return false;
    });
  }, [
    showCloningFeatures, cloningFeatureAnns,
    showCloningReII, showCloningReIIs, showCloningGateway, showCloningPcr,
    enzymeMatchesCutterFilter,
  ]);

  // Merged annotation list: non-cloning + visible cloning + legacy RE sites
  const allAnnotations = useMemo(() => [
    ...nonCloningAnns,
    ...visibleCloningFeatures,
    ...visibleReSites.map((s, i) => ({
      id: `re-${i}`,
      name: s.name,
      start: s.start,
      end: s.end,
      direction: s.direction as 1 | -1 | 0,
      strand: s.direction as 1 | -1 | 0,
      color: s.color,
    })),
  ], [nonCloningAnns, visibleCloningFeatures, visibleReSites]);

  // Non-cutters for the given Type II enabled set (from backend scan).
  const nonCuttersList = useMemo(() => {
    return cloningFeatures?.non_cutters ?? [];
  }, [cloningFeatures]);

  // Handle selection from SeqViz / CircularPlasmidViewer.
  // Linear viewer scrolls so the START of the selection (the start of a
  // feature, or the start of a drag-select) lands at the top of the
  // viewport — so the user sees the beginning of what they selected.
  const handleSelection = useCallback((sel: Selection | null) => {
    selectionRef.current = sel; // Update ref synchronously for double-click handler
    setSelection(sel);
    if (sel) {
      const seqLen = sequenceRef.current?.length ?? 0;
      if (seqLen > 0) {
        const s = ((sel.start % seqLen) + seqLen) % seqLen;
        setTopOnPosition(s);
      }
    }
  }, []);

  useEffect(() => { sequenceRef.current = sequence; }, [sequence]);

  // History tracker — records every sequence/annotation change as an undo step
  // unless skipNextSnapshotRef was set (i.e., the change came from undo/redo).
  useEffect(() => {
    if (skipNextSnapshotRef.current) {
      skipNextSnapshotRef.current = false;
      lastSnapshotRef.current = { seq: sequence, anns: annotations };
      return;
    }
    const prev = lastSnapshotRef.current;
    if (prev && (prev.seq !== sequence || prev.anns !== annotations)) {
      historyPastRef.current.push(prev);
      if (historyPastRef.current.length > 200) historyPastRef.current.shift();
      historyFutureRef.current = []; // any new edit invalidates the redo stack
    }
    lastSnapshotRef.current = { seq: sequence, anns: annotations };
  }, [sequence, annotations]);

  const performUndo = useCallback(() => {
    if (historyPastRef.current.length === 0) return;
    const prev = historyPastRef.current.pop()!;
    historyFutureRef.current.push({ seq: sequence, anns: annotations });
    skipNextSnapshotRef.current = true;
    setSequence(prev.seq);
    setAnnotations(prev.anns);
    setSelection(null);
    selectionRef.current = null;
  }, [sequence, annotations]);

  const performRedo = useCallback(() => {
    if (historyFutureRef.current.length === 0) return;
    const next = historyFutureRef.current.pop()!;
    historyPastRef.current.push({ seq: sequence, anns: annotations });
    skipNextSnapshotRef.current = true;
    setSequence(next.seq);
    setAnnotations(next.anns);
    setSelection(null);
    selectionRef.current = null;
  }, [sequence, annotations]);

  // Insert/replace text at a position; shifts annotations to follow the splice.
  const spliceSequence = useCallback(
    (removeStart: number, removeLen: number, insertText: string) => {
      const newSeq =
        sequence.slice(0, removeStart) +
        insertText +
        sequence.slice(removeStart + removeLen);
      const delta = insertText.length - removeLen;
      const removeEnd = removeStart + removeLen;
      setSequence(newSeq);
      setAnnotations((prev) =>
        prev
          .map((a) => {
            if (a.end <= removeStart) return a;
            if (a.start >= removeEnd) {
              return { ...a, start: a.start + delta, end: a.end + delta };
            }
            const ns = a.start < removeStart ? a.start : removeStart;
            const ne = a.end > removeEnd ? a.end + delta : removeStart;
            if (ne <= ns) return null;
            return { ...a, start: ns, end: ne };
          })
          .filter((a): a is SeqVizAnnotation => a !== null)
      );
      setSelection({
        start: removeStart,
        end: removeStart + insertText.length,
        clockwise: true,
      });
    },
    [sequence]
  );

  // Copy the current selection (handles origin-wrap) onto the clipboard.
  // If the selection exactly matches an annotation with direction = -1
  // (reverse strand), copy the reverse complement so the user gets the
  // gene-sense (5'->3' of the reverse strand) sequence instead of the
  // forward-strand slice.
  const copySelection = useCallback(() => {
    if (!selection || selection.start === selection.end) return;
    const a = selection.start;
    const b = selection.end;
    const seq = sequenceRef.current ?? sequence;
    let sub = a < b ? seq.slice(a, b) : seq.slice(a) + seq.slice(0, b);
    const lo = Math.min(a, b);
    const hi = Math.max(a, b);
    const matchedReverse = annotations.find(
      (ann) => ann.start === lo && ann.end === hi && ann.direction === -1
    );
    if (matchedReverse) {
      const COMP: Record<string, string> = {
        A: "T", T: "A", G: "C", C: "G", N: "N",
        a: "t", t: "a", g: "c", c: "g", n: "n",
      };
      sub = sub.split("").reverse().map((b2) => COMP[b2] ?? b2).join("");
    }
    if (typeof navigator !== "undefined" && navigator.clipboard?.writeText) {
      navigator.clipboard.writeText(sub).catch(() => {});
    }
  }, [selection, sequence, annotations]);

  const pasteFromClipboard = useCallback(async () => {
    if (typeof navigator === "undefined" || !navigator.clipboard?.readText) return;
    let txt = "";
    try {
      txt = await navigator.clipboard.readText();
    } catch {
      return;
    }
    const cleaned = txt.replace(/\s/g, "").toUpperCase();
    if (!cleaned || !/^[ACGTNRYSWKMBDHV]+$/.test(cleaned)) return;
    let removeStart = sequence.length;
    let removeLen = 0;
    if (selection) {
      const a = Math.min(selection.start, selection.end);
      const b = Math.max(selection.start, selection.end);
      removeStart = a;
      removeLen = b - a;
    }
    spliceSequence(removeStart, removeLen, cleaned);
  }, [selection, sequence, spliceSequence]);

  // Find: build matches across forward + reverse complement strand each time
  // the query changes. Returns absolute bp coordinates on the forward strand.
  const findMatches = useMemo(() => {
    const q = findQuery.replace(/\s/g, "").toUpperCase();
    if (!q || q.length < 2 || sequence.length === 0) return [] as { start: number; end: number; reverse: boolean }[];
    const matches: { start: number; end: number; reverse: boolean }[] = [];
    const seq = sequence.toUpperCase();
    // Forward matches
    let fromIdx = 0;
    while (true) {
      const i = seq.indexOf(q, fromIdx);
      if (i < 0) break;
      matches.push({ start: i, end: i + q.length, reverse: false });
      fromIdx = i + 1;
    }
    // Reverse-complement matches: search the rev-comp of q on the forward strand.
    const COMP: Record<string, string> = { A: "T", T: "A", G: "C", C: "G", N: "N", R: "Y", Y: "R", S: "S", W: "W", K: "M", M: "K", B: "V", V: "B", D: "H", H: "D" };
    const rc = q.split("").reverse().map((b) => COMP[b] ?? b).join("");
    if (rc !== q) {
      fromIdx = 0;
      while (true) {
        const i = seq.indexOf(rc, fromIdx);
        if (i < 0) break;
        matches.push({ start: i, end: i + rc.length, reverse: true });
        fromIdx = i + 1;
      }
    }
    matches.sort((a, b) => a.start - b.start);
    return matches;
  }, [findQuery, sequence]);

  // When the active match changes, set selection so the viewers highlight
  // it AND scroll the linear viewer so the top of the match is the first
  // visible row. The circular viewer reads `selection` directly.
  useEffect(() => {
    if (!findOpen || findMatches.length === 0) return;
    const safeIdx = ((findIdx % findMatches.length) + findMatches.length) % findMatches.length;
    const m = findMatches[safeIdx];
    setSelection({ start: m.start, end: m.end, clockwise: !m.reverse });
    selectionRef.current = { start: m.start, end: m.end, clockwise: !m.reverse };
    setTopOnPosition(m.start);
  }, [findOpen, findMatches, findIdx]);

  // Selection from linear viewer: update selection but do NOT re-center it.
  const handleSelectionNoRecenter = useCallback((sel: Selection | null) => {
    selectionRef.current = sel;
    setSelection(sel);
  }, []);

  // Open edit annotation modal
  const openEditAnnotationModal = useCallback((index: number) => {
    const ann = annotations[index];
    // Convert to 1-indexed for display
    setFormData({
      name: ann.name,
      start: String(ann.start + 1),
      end: String(ann.end),
      direction: ann.direction === -1 ? "reverse" : ann.direction === 0 ? "none" : "forward",
      type: ann.type || "",
      description: ann.description || "",
    });
    setFormErrors({});
    setActiveTab("info");  // Default to info tab when opening
    setModalMode("edit");
    setEditingIndex(index);
    setModalOpen(true);
  }, [annotations]);

  // Helper function to find annotation by label text
  // Prioritizes exact matches and longer/more specific matches
  const findAnnotationByLabel = useCallback((labelText: string): number => {
    if (!labelText) return -1;

    // Exact match - highest priority
    const exactIdx = annotations.findIndex(a => a.name === labelText);
    if (exactIdx !== -1) return exactIdx;

    // Collect all potential matches with scores
    const matches: { index: number; score: number }[] = [];

    for (let i = 0; i < annotations.length; i++) {
      const name = annotations[i].name;
      let score = 0;

      // Annotation name equals label (case insensitive)
      if (name.toLowerCase() === labelText.toLowerCase()) {
        score = 1000;
      }
      // Exact prefix match - label is truncated version of annotation name
      else if (name.startsWith(labelText)) {
        // Prefer shorter annotation names (more specific match)
        // e.g., "Fragment_1" should match "Fragment_1" not "Fragment_1 FWD primer"
        score = 500 - name.length;
      }
      // Annotation name starts with label (label is prefix)
      else if (labelText.startsWith(name)) {
        // Prefer longer annotation names (more of the label matched)
        score = 400 + name.length;
      }
      // Label contains full annotation name
      else if (labelText.includes(name)) {
        score = 300 + name.length;
      }
      // Annotation name contains full label
      else if (name.includes(labelText)) {
        score = 200 - (name.length - labelText.length);
      }

      if (score > 0) {
        matches.push({ index: i, score });
      }
    }

    // Sort by score (highest first) and return best match
    if (matches.length > 0) {
      matches.sort((a, b) => b.score - a.score);
      return matches[0].index;
    }

    return -1;
  }, [annotations]);

  // Find annotation by position (start/end match)
  const findAnnotationByPosition = useCallback((start: number, end: number): number => {
    // Find exact match first
    const exactIdx = annotations.findIndex(a => a.start === start && a.end === end);
    if (exactIdx !== -1) return exactIdx;

    // Find annotation that contains this position range
    for (let i = 0; i < annotations.length; i++) {
      const a = annotations[i];
      // Check if the selection is within this annotation
      if (a.start <= start && a.end >= end) {
        return i;
      }
      // Check for origin-crossing annotations
      if (a.start > a.end) {
        if (start >= a.start || end <= a.end) {
          return i;
        }
      }
    }

    // Find closest annotation by start position
    let closestIdx = -1;
    let closestDist = Infinity;
    for (let i = 0; i < annotations.length; i++) {
      const dist = Math.abs(annotations[i].start - start);
      if (dist < closestDist) {
        closestDist = dist;
        closestIdx = i;
      }
    }
    return closestDist <= 10 ? closestIdx : -1;
  }, [annotations]);

  // Handle annotation double-click via DOM event listener
  useEffect(() => {
    const viewerElement = viewerRef.current;
    if (!viewerElement) return;

    const handleDoubleClick = (e: MouseEvent) => {
      const target = e.target as HTMLElement;
      const tagName = target.tagName.toLowerCase();

      // Helper to find label in immediate parent group only (more restrictive)
      const findLabelInImmediateContext = (element: Element): string | null => {
        // Check the element itself
        if (element.tagName.toLowerCase() === 'text' || element.tagName.toLowerCase() === 'tspan') {
          return element.textContent?.trim() || null;
        }

        // Only check immediate parent group (g element)
        const parent = element.parentElement;
        if (parent && parent.tagName.toLowerCase() === 'g') {
          // Look for text element that's a direct child of the same group
          const texts = parent.querySelectorAll(':scope > text');
          for (const text of texts) {
            const content = text.textContent?.trim();
            if (content && content.length > 0) {
              return content;
            }
          }
        }

        return null;
      };

      // Method 1: Check if clicked on text element (annotation label)
      if (tagName === 'text' || tagName === 'tspan') {
        const textContent = target.textContent?.trim() || '';
        const idx = findAnnotationByLabel(textContent);
        if (idx !== -1) {
          e.preventDefault();
          e.stopPropagation();
          openEditAnnotationModal(idx);
          return;
        }
      }

      // Method 2: Use current selection if available (works well for circular view)
      // Use ref instead of state because state update is async and may not be ready yet
      const currentSelection = selectionRef.current;
      if (currentSelection && currentSelection.start !== currentSelection.end) {
        const idx = findAnnotationByPosition(currentSelection.start, currentSelection.end);
        if (idx !== -1) {
          e.preventDefault();
          e.stopPropagation();
          openEditAnnotationModal(idx);
          return;
        }
      }

      // Method 3: Check if clicked on a colored shape (path, rect, polygon)
      if (tagName === 'path' || tagName === 'rect' || tagName === 'polygon') {
        const fill = target.getAttribute('fill');

        // Skip non-annotation elements (white, black, none, or gray backgrounds)
        if (!fill || fill === 'none' || fill === '#ffffff' || fill === '#000000' ||
            fill === 'white' || fill === 'black' || fill.startsWith('rgb(255') ||
            fill.startsWith('rgb(0,') || fill === '#f5f5f5' || fill === '#fafafa') {
          return;
        }

        // Find the label in immediate context only (don't search siblings)
        const label = findLabelInImmediateContext(target);
        if (label) {
          const idx = findAnnotationByLabel(label);
          if (idx !== -1) {
            e.preventDefault();
            e.stopPropagation();
            openEditAnnotationModal(idx);
            return;
          }
        }

        // Try parent's parent (for nested group structures)
        const grandparent = target.parentElement?.parentElement;
        if (grandparent && grandparent.tagName.toLowerCase() === 'g') {
          const texts = grandparent.querySelectorAll(':scope > text, :scope > g > text');
          for (const text of texts) {
            const content = text.textContent?.trim();
            if (content) {
              const idx = findAnnotationByLabel(content);
              if (idx !== -1) {
                e.preventDefault();
                e.stopPropagation();
                openEditAnnotationModal(idx);
                return;
              }
            }
          }
        }
      }

      // Method 4: Traverse up DOM tree looking for annotation group with label
      let clickedElement: Element | null = target.parentElement;
      let depth = 0;
      while (clickedElement && clickedElement !== viewerElement && depth < 4) {
        if (clickedElement.tagName.toLowerCase() === 'g') {
          const texts = clickedElement.querySelectorAll(':scope > text');
          for (const text of texts) {
            const content = text.textContent?.trim();
            if (content) {
              const idx = findAnnotationByLabel(content);
              if (idx !== -1) {
                e.preventDefault();
                e.stopPropagation();
                openEditAnnotationModal(idx);
                return;
              }
            }
          }
        }
        clickedElement = clickedElement.parentElement;
        depth++;
      }
    };

    viewerElement.addEventListener('dblclick', handleDoubleClick);
    return () => viewerElement.removeEventListener('dblclick', handleDoubleClick);
  }, [annotations, openEditAnnotationModal, findAnnotationByLabel, findAnnotationByPosition]);

  // Keyboard handler — Delete/Backspace, Cmd/Ctrl+Z (undo), +Shift+Z or +Y (redo),
  // Cmd/Ctrl+C (copy selection), Cmd/Ctrl+V (paste), Cmd/Ctrl+F (find), Esc (close find).
  useEffect(() => {
    const handleKeyDown = (e: KeyboardEvent) => {
      // Skip when typing in any modal/input — but always honour Cmd+F regardless,
      // and let the find panel close on Esc even if its own input is focused.
      const inEditableField =
        e.target instanceof HTMLInputElement || e.target instanceof HTMLTextAreaElement;
      const inModal = modalOpen || sequenceEditorOpen || addChangeModalOpen;

      const isMeta = e.metaKey || e.ctrlKey;

      // Find — open / focus
      if (isMeta && (e.key === "f" || e.key === "F")) {
        if (inModal) return; // dont steal find inside modal flow
        e.preventDefault();
        setFindOpen(true);
        // Focus on next tick (after render).
        setTimeout(() => findInputRef.current?.focus(), 0);
        return;
      }

      // Esc — close find panel (works whether or not the find input is focused).
      if (e.key === "Escape" && findOpen) {
        e.preventDefault();
        setFindOpen(false);
        return;
      }

      if (inEditableField || inModal) return;

      // Undo / redo
      if (isMeta && (e.key === "z" || e.key === "Z")) {
        e.preventDefault();
        if (e.shiftKey) performRedo();
        else performUndo();
        return;
      }
      if (isMeta && (e.key === "y" || e.key === "Y")) {
        e.preventDefault();
        performRedo();
        return;
      }

      // Copy selected sequence
      if (isMeta && (e.key === "c" || e.key === "C")) {
        if (selection && selection.start !== selection.end) {
          e.preventDefault();
          copySelection();
        }
        return;
      }

      // Paste at the current selection (replacing it) or at end of sequence.
      if (isMeta && (e.key === "v" || e.key === "V")) {
        e.preventDefault();
        pasteFromClipboard();
        return;
      }

      // Delete key with selection - delete selected sequence
      if ((e.key === "Delete" || e.key === "Backspace") && selection && selection.start !== selection.end) {
        e.preventDefault();
        const start = Math.min(selection.start, selection.end);
        const end = Math.max(selection.start, selection.end);

        // Delete the selected region
        const newSeq = sequence.slice(0, start) + sequence.slice(end);
        setSequence(newSeq);

        // Adjust annotations
        const deleteLen = end - start;
        setAnnotations((prev) =>
          prev
            .map((a) => {
              if (a.end <= start) return a;
              if (a.start >= end) {
                return { ...a, start: a.start - deleteLen, end: a.end - deleteLen };
              }
              let newStart = a.start < start ? a.start : start;
              let newEnd = a.end > end ? a.end - deleteLen : start;
              if (newEnd <= newStart) return null;
              return { ...a, start: newStart, end: newEnd };
            })
            .filter((a): a is SeqVizAnnotation => a !== null)
        );

        setSelection(null);
      }
    };

    window.addEventListener("keydown", handleKeyDown);
    return () => window.removeEventListener("keydown", handleKeyDown);
  }, [
    selection,
    sequence,
    modalOpen,
    sequenceEditorOpen,
    addChangeModalOpen,
    findOpen,
    performUndo,
    performRedo,
    copySelection,
    pasteFromClipboard,
  ]);

  // Validate annotation form (1-indexed input)
  const validateForm = useCallback((): boolean => {
    const errors: AnnotationFormErrors = {};

    if (!formData.name.trim()) {
      errors.name = "Name is required";
    }

    const start = parseInt(formData.start, 10);
    const end = parseInt(formData.end, 10);

    if (isNaN(start) || start < 1) {
      errors.start = "Start must be at least 1";
    } else if (start > sequence.length) {
      errors.start = `Start must not exceed sequence length (${sequence.length})`;
    }

    if (isNaN(end) || end < 1) {
      errors.end = "End must be at least 1";
    } else if (end > sequence.length) {
      errors.end = `End must not exceed sequence length (${sequence.length})`;
    }

    setFormErrors(errors);
    return Object.keys(errors).length === 0;
  }, [formData, sequence.length]);

  // Validate DNA sequence
  const validateDnaSequence = useCallback((seq: string): string => {
    const cleaned = seq.replace(/\s/g, "").toUpperCase();
    const invalidChars = cleaned.replace(/[ATGCN]/g, "");
    if (invalidChars.length > 0) {
      const unique = [...new Set(invalidChars)].join(", ");
      return `Invalid characters: ${unique}. Only A, T, G, C, N allowed.`;
    }
    return "";
  }, []);

  // Open the direction-prompt dialog for creating a translation annotation
  // from the current selection. Default direction follows the selection's
  // existing clockwise flag (forward when clockwise, reverse otherwise).
  const handleOpenTranslationModal = useCallback(() => {
    if (!selection || selection.start === selection.end) return;
    setTranslationDirection(selection.clockwise === false ? "reverse" : "forward");
    setTranslationModalOpen(true);
  }, [selection]);

  const handleCreateTranslation = useCallback(() => {
    if (!selection || selection.start === selection.end) return;
    const start = Math.min(selection.start, selection.end);
    const end = Math.max(selection.start, selection.end);
    const dir: 1 | -1 = translationDirection === "reverse" ? -1 : 1;
    const region = sequence.slice(start, end);
    const CODON_TABLE: Record<string, string> = {
      TTT:"F",TTC:"F",TTA:"L",TTG:"L",CTT:"L",CTC:"L",CTA:"L",CTG:"L",
      ATT:"I",ATC:"I",ATA:"I",ATG:"M",GTT:"V",GTC:"V",GTA:"V",GTG:"V",
      TCT:"S",TCC:"S",TCA:"S",TCG:"S",CCT:"P",CCC:"P",CCA:"P",CCG:"P",
      ACT:"T",ACC:"T",ACA:"T",ACG:"T",GCT:"A",GCC:"A",GCA:"A",GCG:"A",
      TAT:"Y",TAC:"Y",TAA:"*",TAG:"*",CAT:"H",CAC:"H",CAA:"Q",CAG:"Q",
      AAT:"N",AAC:"N",AAA:"K",AAG:"K",GAT:"D",GAC:"D",GAA:"E",GAG:"E",
      TGT:"C",TGC:"C",TGA:"*",TGG:"W",CGT:"R",CGC:"R",CGA:"R",CGG:"R",
      AGT:"S",AGC:"S",AGA:"R",AGG:"R",GGT:"G",GGC:"G",GGA:"G",GGG:"G",
    };
    const COMP: Record<string, string> = { A:"T",T:"A",G:"C",C:"G",N:"N" };
    const dnaForTranslate = dir === 1
      ? region.toUpperCase()
      : region.toUpperCase().split("").reverse().map((b) => COMP[b] || "N").join("");
    let aa = "";
    for (let i = 0; i + 2 < dnaForTranslate.length; i += 3) {
      aa += CODON_TABLE[dnaForTranslate.slice(i, i + 3)] || "X";
    }
    if (aa.endsWith("*")) aa = aa.slice(0, -1);

    const newAnn: SeqVizAnnotation = {
      name: `Translation (${aa.length} aa)`,
      start,
      end,
      direction: dir,
      strand: dir,
      length: end - start,
      color: "#673AB7",
      layer: "translation",
      module_type: "translation",
      source: "user_selection",
      metadata: {
        aa_length: aa.length,
        aa_sequence: aa,
        feature_regions: [],
        orf_detected: false,
        user_created: true,
      },
    };
    setAnnotations((prev) => [...prev, newAnn]);
    setTranslationModalOpen(false);
    setSelection(null);
  }, [selection, translationDirection, sequence]);

  // Open add annotation modal
  const handleAddAnnotation = useCallback(() => {
    const defaultStart = selection ? String(Math.min(selection.start, selection.end) + 1) : "";
    const defaultEnd = selection ? String(Math.max(selection.start, selection.end)) : "";
    // Find/search sets selection.clockwise=false on reverse-strand matches —
    // surface that as direction="reverse" so the user does not have to flip
    // it manually when adding an annotation for a reverse-strand match.
    const defaultDirection: "forward" | "reverse" | "none" =
      selection && selection.clockwise === false ? "reverse" : "forward";

    setFormData({
      name: "",
      start: defaultStart,
      end: defaultEnd,
      direction: defaultDirection,
    });
    setFormErrors({});
    setModalMode("add");
    setEditingIndex(null);
    setModalOpen(true);
  }, [selection]);

  // Save annotation (add or edit)
  const handleSaveAnnotation = useCallback(() => {
    if (!validateForm()) return;

    const start = parseInt(formData.start, 10) - 1;
    const end = parseInt(formData.end, 10);
    const direction = formData.direction === "reverse" ? -1 : formData.direction === "none" ? 0 : 1;

    const newAnnotation: SeqVizAnnotation = {
      name: formData.name.trim(),
      start,
      end,
      direction: direction as 1 | -1 | 0,
      strand: direction as 1 | -1 | 0,
      length: start <= end ? end - start : (sequence.length - start) + end,
    };

    if (modalMode === "add") {
      setAnnotations((prev) => [...prev, newAnnotation]);
    } else if (editingIndex !== null) {
      setAnnotations((prev) => {
        const updated = [...prev];
        updated[editingIndex] = newAnnotation;
        return updated;
      });
    }

    setModalOpen(false);
    setEditingIndex(null);
    setSelection(null);
  }, [formData, modalMode, editingIndex, sequence.length, validateForm]);

  // Delete annotation
  const handleDeleteAnnotation = useCallback((index: number) => {
    setAnnotations((prev) => prev.filter((_, i) => i !== index));
    setModalOpen(false);
    setEditingIndex(null);
  }, []);

  // Open sequence editor
  const handleEditSequence = useCallback(() => {
    setEditedSequence(sequence);
    setSequenceError("");
    setSequenceEditorOpen(true);
  }, [sequence]);

  // Save sequence from full editor
  const handleSaveSequence = useCallback(() => {
    const cleaned = editedSequence.replace(/\s/g, "").toUpperCase();
    const error = validateDnaSequence(cleaned);
    if (error) {
      setSequenceError(error);
      return;
    }

    const newLength = cleaned.length;
    setSequence(cleaned);
    setSequenceEditorOpen(false);

    if (newLength !== sequence.length) {
      setAnnotations((prev) =>
        prev
          .map((a) => ({
            ...a,
            start: Math.min(a.start, newLength - 1),
            end: Math.min(a.end, newLength),
          }))
          .filter((a) => a.end > a.start)
      );
    }
  }, [editedSequence, sequence.length, validateDnaSequence]);

  // Open Add/Change Sequence modal
  const handleOpenAddChangeModal = useCallback(() => {
    if (selection && selection.start !== selection.end) {
      // Range selected - show current sequence
      const start = Math.min(selection.start, selection.end);
      const end = Math.max(selection.start, selection.end);
      setNewSequenceText(sequence.slice(start, end));
    } else {
      // No range - empty for insertion
      setNewSequenceText("");
    }
    setAddChangeError("");
    setAddChangeModalOpen(true);
  }, [selection, sequence]);

  // Apply Add/Change Sequence
  const handleApplyAddChange = useCallback(() => {
    const cleaned = newSequenceText.replace(/\s/g, "").toUpperCase();
    const error = validateDnaSequence(cleaned);
    if (error && cleaned.length > 0) {
      setAddChangeError(error);
      return;
    }

    if (selection && selection.start !== selection.end) {
      // Replace selected range
      const start = Math.min(selection.start, selection.end);
      const end = Math.max(selection.start, selection.end);
      const oldLen = end - start;
      const newLen = cleaned.length;
      const diff = newLen - oldLen;

      const newSeq = sequence.slice(0, start) + cleaned + sequence.slice(end);
      setSequence(newSeq);

      // Adjust annotations
      setAnnotations((prev) =>
        prev
          .map((a) => {
            if (a.end <= start) return a;
            if (a.start >= end) {
              return { ...a, start: a.start + diff, end: a.end + diff };
            }
            // Overlapping annotation
            let newStart = a.start;
            let newEnd = a.end;
            if (a.start >= start) newStart = start;
            if (a.end <= end) {
              newEnd = start + newLen;
            } else {
              newEnd = a.end + diff;
            }
            if (newEnd <= newStart) return null;
            return { ...a, start: newStart, end: newEnd };
          })
          .filter((a): a is SeqVizAnnotation => a !== null)
      );
    } else {
      // Insert at cursor position (default to end if no selection)
      const insertPos = selection ? selection.start : sequence.length;
      const insertLen = cleaned.length;

      if (insertLen === 0) {
        setAddChangeModalOpen(false);
        return;
      }

      const newSeq = sequence.slice(0, insertPos) + cleaned + sequence.slice(insertPos);
      setSequence(newSeq);

      // Adjust annotations
      setAnnotations((prev) =>
        prev.map((a) => {
          if (a.end <= insertPos) return a;
          if (a.start >= insertPos) {
            return { ...a, start: a.start + insertLen, end: a.end + insertLen };
          }
          return { ...a, end: a.end + insertLen };
        })
      );
    }

    setAddChangeModalOpen(false);
    setSelection(null);
  }, [newSequenceText, selection, sequence, validateDnaSequence]);

  // Delete selected sequence (from modal)
  const handleDeleteSelected = useCallback(() => {
    if (!selection || selection.start === selection.end) return;

    const start = Math.min(selection.start, selection.end);
    const end = Math.max(selection.start, selection.end);
    const deleteLen = end - start;

    const newSeq = sequence.slice(0, start) + sequence.slice(end);
    setSequence(newSeq);

    setAnnotations((prev) =>
      prev
        .map((a) => {
          if (a.end <= start) return a;
          if (a.start >= end) {
            return { ...a, start: a.start - deleteLen, end: a.end - deleteLen };
          }
          let newStart = a.start < start ? a.start : start;
          let newEnd = a.end > end ? a.end - deleteLen : start;
          if (newEnd <= newStart) return null;
          return { ...a, start: newStart, end: newEnd };
        })
        .filter((a): a is SeqVizAnnotation => a !== null)
    );

    setAddChangeModalOpen(false);
    setSelection(null);
  }, [selection, sequence]);

  // Close modals
  const handleCloseModal = useCallback(() => {
    setModalOpen(false);
    setEditingIndex(null);
    setFormErrors({});
  }, []);

  const handleCloseSequenceEditor = useCallback(() => {
    setSequenceEditorOpen(false);
    setSequenceError("");
  }, []);

  const handleCloseAddChangeModal = useCallback(() => {
    setAddChangeModalOpen(false);
    setAddChangeError("");
  }, []);

  // Import Annotations from CSV
  const handleOpenImportModal = useCallback(() => {
    setImportError("");
    setImportResult(null);
    setImportMaxMismatches(0);
    setImportFileName("");
    if (importFileInputRef.current) importFileInputRef.current.value = "";
    setImportModalOpen(true);
  }, []);

  const handleClickImportFile = useCallback(() => {
    importFileInputRef.current?.click();
  }, []);

  const handleImportFileSelected = useCallback((e: React.ChangeEvent<HTMLInputElement>) => {
    const f = e.target.files?.[0];
    setImportFileName(f ? f.name : "");
  }, []);

  const handleCloseImportModal = useCallback(() => {
    setImportModalOpen(false);
    setImportError("");
  }, []);

  const handleDownloadImportTemplate = useCallback(() => {
    const csv = [
      "name,sequence,type,location,length,description",
      "example_guide,GCTCAAGATCGTCCTTCCAA,primer,,,Example sgRNA-style oligo",
      "example_misc,ATGGCCGGCAAACTG,misc_feature,,,Example longer feature",
      "",
    ].join("\n");
    const blob = new Blob([csv], { type: "text/csv;charset=utf-8" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = "annotation_import_template.csv";
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  }, []);

  // Minimal CSV parser: splits on newlines and commas. Trims whitespace and
  // strips wrapping double-quotes from each cell. Sequence/name fields are
  // not expected to contain embedded commas or newlines.
  const parseImportCsv = useCallback((text: string): Array<Record<string, string>> => {
    const rows: Array<Record<string, string>> = [];
    const lines = text.split(/\r?\n/).filter((l) => l.trim().length > 0);
    if (lines.length === 0) return rows;
    const splitLine = (line: string) =>
      line.split(",").map((c) => c.trim().replace(/^"(.*)"$/, "$1"));
    const header = splitLine(lines[0]).map((h) => h.toLowerCase());
    for (let i = 1; i < lines.length; i++) {
      const cells = splitLine(lines[i]);
      const row: Record<string, string> = {};
      header.forEach((h, j) => { row[h] = cells[j] ?? ""; });
      rows.push(row);
    }
    return rows;
  }, []);

  const handleSubmitImport = useCallback(async () => {
    if (!sequence) {
      setImportError("Load a plasmid sequence first.");
      return;
    }
    const file = importFileInputRef.current?.files?.[0];
    if (!file) {
      setImportError("Choose a CSV file to import.");
      return;
    }
    setImportLoading(true);
    setImportError("");
    setImportResult(null);
    try {
      const text = await file.text();
      const rows = parseImportCsv(text);
      if (rows.length === 0) {
        setImportError("CSV is empty or has only a header row.");
        setImportLoading(false);
        return;
      }
      const entries = rows
        .map((r) => ({
          name: r["name"] || "",
          sequence: r["sequence"] || "",
          type: r["type"] || undefined,
          location: r["location"] || undefined,
          length: r["length"] || undefined,
          description: r["description"] || undefined,
        }))
        .filter((e) => e.name || e.sequence);
      if (entries.length === 0) {
        setImportError("No rows with a name or sequence found. Check headers (name, sequence required).");
        setImportLoading(false);
        return;
      }

      const res = await fetch("/api/plannotate/import-annotations", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          sequence,
          circular,
          entries,
          max_mismatches: importMaxMismatches,
        }),
      });
      const data = await res.json();
      if (!data.ok) {
        setImportError(data.error || "Import failed");
        setImportLoading(false);
        return;
      }

      const newAnns: SeqVizAnnotation[] = (data.annotations || []).map((a: any) => ({
        name: a.name,
        start: a.start,
        end: a.end,
        direction: a.direction === -1 ? -1 : 1,
        strand: a.direction === -1 ? -1 : 1,
        color: a.color || "#84B0DC",
        type: a.type || "misc_feature",
        layer: "feature" as const,
        description: a.description || "",
        kb_data: a.kb_data || null,
        source: a.source || "csv_import",
      }));

      // Idempotent merge: dedupe the FULL combined list by
      // (name, start, end, direction, source). This is robust against
      // React Strict Mode invoking the updater twice, against double-clicks,
      // and against stale state from a prior import session — the resulting
      // list will never contain two annotations with the same key.
      setAnnotations((prev) => {
        const key = (a: SeqVizAnnotation) =>
          `${a.name}\u0001${a.start}\u0001${a.end}\u0001${a.direction ?? 0}\u0001${a.source ?? ""}`;
        const seen = new Set<string>();
        const out: SeqVizAnnotation[] = [];
        for (const a of [...prev, ...newAnns]) {
          const k = key(a);
          if (seen.has(k)) continue;
          seen.add(k);
          out.push(a);
        }
        return out;
      });
      setImportResult({ summary: data.summary, unmatched: data.unmatched || [] });
    } catch (e: any) {
      setImportError(e.message || "Network error");
    } finally {
      setImportLoading(false);
    }
  }, [sequence, circular, importMaxMismatches, parseImportCsv]);

  // Application-driven default size + Tm + exclusion windows. The
  // anchor feature is the annotation whose name (and span) seeds the
  // pair label and the excluded region; if no feature is supplied we
  // try to derive one from the current selection.
  const PRIMER_APP_DEFAULTS: Record<string, { sizeMin: string; sizeMax: string; pad: number }> = {
    fragment: { sizeMin: "100", sizeMax: "300", pad: 0 },
    sanger:   { sizeMin: "250", sizeMax: "500", pad: 75 },
    illumina: { sizeMin: "150", sizeMax: "290", pad: 30 },
  };
  const computePrimerDefaults = useCallback(
    (
      app: "fragment" | "sanger" | "illumina",
      anchor?: SeqVizAnnotation | null,
    ) => {
      const seqLen = sequence ? sequence.length : 0;
      const cfg = PRIMER_APP_DEFAULTS[app] || PRIMER_APP_DEFAULTS.fragment;

      // Resolve the anchor feature: explicit > derived from selection.
      let feat: SeqVizAnnotation | null = anchor || null;
      if (!feat && selection && selection.start !== selection.end) {
        const idx = findAnnotationByPosition(
          Math.min(selection.start, selection.end),
          Math.max(selection.start, selection.end),
        );
        if (idx >= 0) feat = annotations[idx];
      }

      // Region: prefer selection (1-indexed), else feature span ± a
      // generous halo, else the full sequence.
      let regionStart1 = 1;
      let regionEnd1 = seqLen;
      if (selection && selection.start !== selection.end) {
        regionStart1 = Math.min(selection.start, selection.end) + 1;
        regionEnd1 = Math.max(selection.start, selection.end);
      } else if (feat) {
        regionStart1 = Math.max(1, feat.start + 1 - 250);
        regionEnd1 = Math.min(seqLen, feat.end + 250);
      }

      // Excluded region around the anchor feature, padded per app.
      let exStart1 = "";
      let exEnd1 = "";
      if (feat && cfg.pad > 0) {
        const featStart1 = feat.start + 1;
        const featEnd1 = feat.end;
        const ex1 = Math.max(regionStart1, featStart1 - cfg.pad);
        const ex2 = Math.min(regionEnd1, featEnd1 + cfg.pad);
        if (ex2 > ex1) {
          exStart1 = String(ex1);
          exEnd1 = String(ex2);
          // Pull the design region out beyond the exclusion so primer3
          // has flanking room to place primers in.
          regionStart1 = Math.max(1, ex1 - 200);
          regionEnd1 = Math.min(seqLen, ex2 + 200);
        }
      }

      return {
        region_start: String(regionStart1),
        region_end: String(regionEnd1),
        excluded_start: exStart1,
        excluded_end: exEnd1,
        product_size_min: cfg.sizeMin,
        product_size_max: cfg.sizeMax,
        primer_opt_tm: "60",
        pair_label: feat?.name || "",
      };
    },
    [sequence, selection, annotations, findAnnotationByPosition],
  );

  // Primer design: open / close
  const handleOpenPrimerModal = useCallback(() => {
    setPrimerError("");
    setPrimerResult(null);
    setPrimerForm((f) => {
      const d = computePrimerDefaults(f.application);
      return {
        ...f,
        region_start: d.region_start,
        region_end: d.region_end,
        excluded_start: d.excluded_start,
        excluded_end: d.excluded_end,
        product_size_min: d.product_size_min,
        product_size_max: d.product_size_max,
        primer_opt_tm: d.primer_opt_tm,
        pair_label: d.pair_label,
        primer_fwd_name: "",
        primer_rev_name: "",
      };
    });
    setPrimerModalOpen(true);
  }, [computePrimerDefaults]);

  // Shortcut: design flanking primers around an sgRNA-type annotation.
  // Pre-fills the primer form with: region = sgRNA ± 250 bp, excluded
  // region = sgRNA ± 75 bp, Tm opt = 60 °C, pair label = "<name>_amplicon",
  // and explicit forward / reverse primer names of "<name>_Fwd" / "<name>_Rev".
  const handleOpenPrimerModalForSgrna = useCallback((ann: SeqVizAnnotation) => {
    if (!sequence) return;
    setPrimerError("");
    setPrimerResult(null);
    setPrimerForm((f) => {
      const d = computePrimerDefaults(f.application, ann);
      return {
        ...f,
        region_start: d.region_start,
        region_end: d.region_end,
        excluded_start: d.excluded_start,
        excluded_end: d.excluded_end,
        product_size_min: d.product_size_min,
        product_size_max: d.product_size_max,
        primer_opt_tm: d.primer_opt_tm,
        pair_label: ann.name,
        primer_fwd_name: "",
        primer_rev_name: "",
      };
    });
    // Close the annotation info modal (if open) before opening the design modal
    setModalOpen(false);
    setEditingIndex(null);
    setPrimerModalOpen(true);
  }, [sequence, computePrimerDefaults]);

  const handleClosePrimerModal = useCallback(() => {
    setPrimerModalOpen(false);
    setPrimerError("");
  }, []);

  const handlePrimerFormChange = useCallback(
    (field: keyof typeof primerForm) =>
      (e: React.ChangeEvent<HTMLInputElement | HTMLSelectElement>) => {
      const value = e.target.value;
      if (field === "application") {
        const app = (value as "fragment" | "sanger" | "illumina");
        const d = computePrimerDefaults(app);
        setPrimerForm((f) => ({
          ...f,
          application: app,
          // Re-seed only the size / Tm / exclusion bands (and the pair
          // label if the user has not explicitly typed one). We do NOT
          // overwrite explicit fwd/rev primer names.
          excluded_start: d.excluded_start,
          excluded_end: d.excluded_end,
          product_size_min: d.product_size_min,
          product_size_max: d.product_size_max,
          primer_opt_tm: d.primer_opt_tm,
          pair_label: f.pair_label || d.pair_label,
        }));
        return;
      }
      setPrimerForm((f) => ({ ...f, [field]: value }));
    },
    [computePrimerDefaults]
  );

  const handleSubmitPrimerDesign = useCallback(async () => {
    if (!sequence) {
      setPrimerError("Load a plasmid sequence first.");
      return;
    }
    const seqLen = sequence.length;
    const r1 = parseInt(primerForm.region_start, 10);
    const r2 = parseInt(primerForm.region_end, 10);
    if (!Number.isFinite(r1) || !Number.isFinite(r2) || r1 < 1 || r2 > seqLen || r1 >= r2) {
      setPrimerError(`Region must be 1-indexed within 1..${seqLen} with start < end.`);
      return;
    }
    // Slice the chosen region as the primer3 template; result coordinates
    // get translated back to plasmid space below.
    const offset = r1 - 1;
    const template = sequence.slice(offset, r2);

    let excludedPayload: { excluded_start: number; excluded_length: number } | undefined;
    if (primerForm.excluded_start.trim() && primerForm.excluded_end.trim()) {
      const e1 = parseInt(primerForm.excluded_start, 10);
      const e2 = parseInt(primerForm.excluded_end, 10);
      if (!Number.isFinite(e1) || !Number.isFinite(e2) || e1 < r1 || e2 > r2 || e1 >= e2) {
        setPrimerError(`Excluded region must lie inside the design region (${r1}..${r2}).`);
        return;
      }
      excludedPayload = {
        excluded_start: e1 - r1,           // offset relative to template
        excluded_length: e2 - e1 + 1,
      };
    }

    setPrimerLoading(true);
    setPrimerError("");
    setPrimerResult(null);

    const numericOrUndef = (s: string): number | undefined => {
      const n = parseFloat(s);
      return Number.isFinite(n) ? n : undefined;
    };
    const intOrUndef = (s: string): number | undefined => {
      const n = parseInt(s, 10);
      return Number.isFinite(n) ? n : undefined;
    };

    const productMin = intOrUndef(primerForm.product_size_min) ?? 100;
    const productMax = intOrUndef(primerForm.product_size_max) ?? 300;
    if (primerForm.application === "illumina" && productMax > ILLUMINA_MAX_AMPLICON_BP) {
      setPrimerError(
        "Illumina primer design is limited to amplicons \u2264 " +
        ILLUMINA_MAX_AMPLICON_BP +
        " bp (current max: " + productMax +
        " bp). Reduce product_size_max or pick a different application."
      );
      setPrimerLoading(false);
      return;
    }
    const body: any = {
      fragments_in: template,
      application: primerForm.application,
      product_size_min: productMin,
      product_size_max: productMax,
      num_return: intOrUndef(primerForm.num_return) ?? 5,
    };
    if (excludedPayload) Object.assign(body, excludedPayload);
    const tmin = numericOrUndef(primerForm.primer_min_tm);
    const topt = numericOrUndef(primerForm.primer_opt_tm);
    const tmax = numericOrUndef(primerForm.primer_max_tm);
    if (tmin !== undefined) body.primer_min_tm = tmin;
    if (topt !== undefined) body.primer_opt_tm = topt;
    if (tmax !== undefined) body.primer_max_tm = tmax;
    const smin = intOrUndef(primerForm.primer_min_size);
    const sopt = intOrUndef(primerForm.primer_opt_size);
    const smax = intOrUndef(primerForm.primer_max_size);
    if (smin !== undefined) body.primer_min_size = smin;
    if (sopt !== undefined) body.primer_opt_size = sopt;
    if (smax !== undefined) body.primer_max_size = smax;

    try {
      const res = await fetch("/api/pcr/design-primers", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(body),
      });
      const data = await res.json();
      if (!res.ok || data.ok === false) {
        setPrimerError(data.error || data.detail || "Primer design failed");
        setPrimerLoading(false);
        return;
      }

      // Translate primer positions back to plasmid coordinates.
      // primer3 returns: PRIMER_LEFT  -> [start, length]   (5' on forward strand)
      //                  PRIMER_RIGHT -> [3'-end, length]  (3' on forward strand for the rev primer)
      const leftStartTpl = data.left_pos?.start;
      const leftLen = data.left_pos?.len;
      const rightThreeTpl = data.right_pos?.start_3prime;
      const rightLen = data.right_pos?.len;
      if (
        typeof leftStartTpl !== "number" || typeof leftLen !== "number" ||
        typeof rightThreeTpl !== "number" || typeof rightLen !== "number"
      ) {
        setPrimerError("primer3 returned malformed position data");
        setPrimerLoading(false);
        return;
      }
      const leftStartPlasmid = leftStartTpl + offset;
      const leftEndPlasmid = leftStartPlasmid + leftLen;            // half-open
      const rightStartPlasmid = (rightThreeTpl - rightLen + 1) + offset;
      const rightEndPlasmid = rightThreeTpl + 1 + offset;            // half-open
      const ampliconStart = leftStartPlasmid;
      const ampliconEnd = rightEndPlasmid;
      const amplicon = sequence.slice(ampliconStart, ampliconEnd);

      const labelPrefix = primerForm.pair_label.trim() || `primer_${Date.now().toString(36).slice(-4)}`;
      const sfx = PRIMER_NAME_SUFFIX[primerForm.application] || PRIMER_NAME_SUFFIX.fragment;
      const fwdName = primerForm.primer_fwd_name.trim() || `${labelPrefix}${sfx.fwd}`;
      const revName = primerForm.primer_rev_name.trim() || `${labelPrefix}${sfx.rev}`;
      const sharedMeta = {
        primer_design: true,
        pair_index: data.pair_index,
        pair_penalty: data.pair_penalty,
        product_size: data.product_size,
        amplicon_start_1based: ampliconStart + 1,
        amplicon_end_1based: ampliconEnd,
        amplicon_length: amplicon.length,
        amplicon: amplicon,
        excluded_region: data.excluded_region || null,
        design_region_1based: `${r1}..${r2}`,
      };
      const leftAnn: SeqVizAnnotation = {
        name: fwdName,
        start: leftStartPlasmid,
        end: leftEndPlasmid,
        direction: 1,
        strand: 1,
        color: "#84B0DC",
        type: "primer",
        layer: "feature",
        description: `Forward primer (Primer3)`,
        source: "primer3_design",
        metadata: {
          ...sharedMeta,
          primer_role: "forward",
          primer_sequence: data.left_primer,
          primer_tm: data.left_tm,
          mispriming_sites: data.left_mispriming_sites || [],
          thermo_scores: data.left_scores || {},
        },
      };
      const rightAnn: SeqVizAnnotation = {
        name: revName,
        start: rightStartPlasmid,
        end: rightEndPlasmid,
        direction: -1,
        strand: -1,
        color: "#F58A5E",
        type: "primer",
        layer: "feature",
        description: `Reverse primer (Primer3)`,
        source: "primer3_design",
        metadata: {
          ...sharedMeta,
          primer_role: "reverse",
          primer_sequence: data.right_primer,
          primer_tm: data.right_tm,
          mispriming_sites: data.right_mispriming_sites || [],
          thermo_scores: data.right_scores || {},
        },
      };

      setAnnotations((prev) => {
        const key = (a: SeqVizAnnotation) =>
          `${a.name}\u0001${a.start}\u0001${a.end}\u0001${a.direction ?? 0}\u0001${a.source ?? ""}`;
        const seen = new Set<string>();
        const out: SeqVizAnnotation[] = [];
        for (const a of [...prev, leftAnn, rightAnn]) {
          const k = key(a);
          if (seen.has(k)) continue;
          seen.add(k);
          out.push(a);
        }
        return out;
      });
      setPrimerResult({
        application: data.application || primerForm.application,
        selection_method: data.selection_method,
        selection_rationale: data.selection_rationale,
        candidate_scores: data.candidate_scores || [],
        num_candidates_considered: data.num_candidates_considered,
        left_primer: data.left_primer,        // full ordered primer (adapter+anneal for Illumina)
        right_primer: data.right_primer,
        left_annealing: data.left_annealing || data.left_primer,
        right_annealing: data.right_annealing || data.right_primer,
        left_adapter: data.left_adapter || "",
        right_adapter: data.right_adapter || "",
        left_tm: data.left_tm,                // annealing-only Tm
        right_tm: data.right_tm,
        product_size: data.product_size,
        amplicon,
        amplicon_start_1based: ampliconStart + 1,
        amplicon_end_1based: ampliconEnd,
        pair_penalty: data.pair_penalty,
        sanger_scores: data.sanger_scores || [],
        amplicon_name: labelPrefix,
        fwd_name: fwdName,
        rev_name: revName,
        left_mispriming_sites: data.left_mispriming_sites || [],
        right_mispriming_sites: data.right_mispriming_sites || [],
        left_scores: data.left_scores || {},
        right_scores: data.right_scores || {},
        design_region_1based: `${r1}..${r2}`,
      });
    } catch (e: any) {
      setPrimerError(e.message || "Network error");
    } finally {
      setPrimerLoading(false);
    }
  }, [sequence, primerForm]);

  // ---- Guide design (sgRNA) ----
  const handleOpenGuideModal = useCallback(() => {
    setGuideError("");
    setGuideResult(null);
    const sel = selection;
    const selStart = sel ? Math.min(sel.start, sel.end) : 0;
    const selEnd = sel ? Math.max(sel.start, sel.end) : (sequence ? sequence.length : 0);
    setGuideForm((f) => ({
      ...f,
      region_start: String(selStart + 1),
      region_end: String(selEnd),
    }));
    setGuideModalOpen(true);
  }, [selection, sequence]);

  const handleCloseGuideModal = useCallback(() => {
    setGuideModalOpen(false);
    setGuideError("");
  }, []);

  const handleGuideFormChange = useCallback(
    (field: keyof typeof guideForm) =>
      (e: React.ChangeEvent<HTMLInputElement | HTMLSelectElement>) => {
        setGuideForm((f) => ({ ...f, [field]: e.target.value as any }));
      },
    []
  );

  const handleSubmitGuideDesign = useCallback(async () => {
    if (!sequence) {
      setGuideError("Load a plasmid sequence first.");
      return;
    }
    const seqLen = sequence.length;
    const r1 = parseInt(guideForm.region_start, 10);
    const r2 = parseInt(guideForm.region_end, 10);
    if (!Number.isFinite(r1) || !Number.isFinite(r2) || r1 < 1 || r2 > seqLen || r1 >= r2) {
      setGuideError(`Region must be 1-indexed within 1..${seqLen} with start < end.`);
      return;
    }
    const guide_length = parseInt(guideForm.guide_length, 10);
    if (!Number.isFinite(guide_length) || guide_length < 16 || guide_length > 25) {
      setGuideError("Guide length must be between 16 and 25 nt.");
      return;
    }
    const max_guides = parseInt(guideForm.max_guides, 10) || 20;
    const min_score = parseFloat(guideForm.min_score) || 0;
    const pam = guideForm.pam.trim().toUpperCase();
    if (!/^[ACGTRYSWKMBDHVN]+$/.test(pam)) {
      setGuideError("PAM must use IUPAC bases (A C G T R Y S W K M B D H V N).");
      return;
    }
    setGuideLoading(true);
    setGuideError("");
    setGuideResult(null);
    try {
      const res = await fetch("/api/plannotate/design-guides", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          sequence,
          region_start: r1,
          region_end: r2,
          pam,
          guide_length,
          pam_position: guideForm.pam_position,
          max_guides,
          min_score,
          score_method: guideForm.score_method,
        }),
      });
      const data = await res.json();
      if (!res.ok || data.ok === false) {
        setGuideError(data.error || data.detail || "Guide design failed");
        setGuideLoading(false);
        return;
      }
      const guides = (data.guides || []) as any[];
      const colorForScore = (s: number) => {
        const t = Math.max(0, Math.min(1, s / 100));
        const r = Math.round(220 + (132 - 220) * t);
        const g = Math.round(110 + (199 - 110) * t);
        const b = Math.round(120 + (140 - 120) * t);
        return `rgb(${r},${g},${b})`;
      };
      const newAnns: SeqVizAnnotation[] = guides.map((g: any) => ({
        name: g.name,
        start: g.start,
        end: g.end,
        direction: (g.direction === -1 ? -1 : 1) as 1 | -1,
        strand: (g.direction === -1 ? -1 : 1) as 1 | -1,
        color: colorForScore(g.score),
        type: "sgRNA",
        layer: "feature" as const,
        description: `sgRNA · ${g.score_method || data.summary?.score_method || "score"} ${g.score} · GC ${(g.gc_fraction * 100).toFixed(0)}% · ${g.n_offtargets} off-target${g.n_offtargets === 1 ? "" : "s"}`,
        source: "guide_design",
        metadata: {
          guide_design: true,
          spacer: g.spacer,
          pam: g.pam,
          score: g.score,
          score_method: g.score_method || data.summary?.score_method,
          score_components: g.score_components,
          context_30mer: g.context_30mer,
          gc_fraction: g.gc_fraction,
          max_homopolymer: g.max_homopolymer,
          n_offtargets: g.n_offtargets,
          design_region_1based: data.summary?.region_1based,
          pam_setting: data.summary?.pam,
          pam_position: data.summary?.pam_position,
        },
      }));
      // Compute the additions count from the *current* annotations state
      // (closure capture). React batches setState calls so reading nAdded
      // inside the updater races with setGuideResult; doing the dedupe up
      // front gives the popup a correct count regardless of when the
      // updater fires.
      const key = (a: SeqVizAnnotation) =>
        `${a.name}\u0001${a.start}\u0001${a.end}\u0001${a.direction ?? 0}\u0001${a.source ?? ""}`;
      const existingKeys = new Set(annotations.map(key));
      const additions = newAnns.filter((a) => !existingKeys.has(key(a)));
      const nAdded = additions.length;
      setAnnotations((prev) => {
        const seen = new Set<string>();
        const out: SeqVizAnnotation[] = [];
        for (const a of [...prev, ...newAnns]) {
          const k = key(a);
          if (seen.has(k)) continue;
          seen.add(k);
          out.push(a);
        }
        return out;
      });
      setGuideResult({
        summary: data.summary,
        n_added: nAdded,
        guides,
        design_region_1based: data.summary?.region_1based || `${r1}..${r2}`,
      });
    } catch (e: any) {
      setGuideError(e.message || "Network error");
    } finally {
      setGuideLoading(false);
    }
  }, [sequence, guideForm, annotations]);

  // ---- pegRNA design (prime editing — easy_prime port) ----
  const handleOpenPegrnaModal = useCallback(() => {
    setPegrnaError("");
    setPegrnaResult(null);
    const sel = selection;
    const s = sel ? Math.min(sel.start, sel.end) : 0;
    const e = sel ? Math.max(sel.start, sel.end) : (sequence ? sequence.length : 0);
    const refSlice = sequence ? sequence.slice(s, Math.max(s + 1, e)) : "";
    setPegrnaForm((f) => ({
      ...f,
      edit_start: String(s + 1),
      edit_end: String(Math.max(s + 1, e)),
      alt: f.alt || refSlice,
    }));
    setPegrnaModalOpen(true);
  }, [selection, sequence]);

  const handleClosePegrnaModal = useCallback(() => {
    setPegrnaModalOpen(false);
    setPegrnaError("");
  }, []);

  const handlePegrnaFormChange = useCallback(
    (field: keyof typeof pegrnaForm) =>
      (e: React.ChangeEvent<HTMLInputElement | HTMLSelectElement>) => {
        const v = e.target.type === "checkbox"
          ? (e.target as HTMLInputElement).checked
          : (e.target.value as any);
        setPegrnaForm((f) => ({ ...f, [field]: v }));
      },
    []
  );

  const handleSubmitPegrnaDesign = useCallback(async () => {
    if (!sequence) {
      setPegrnaError("Load a plasmid sequence first.");
      return;
    }
    const seqLen = sequence.length;
    const s = parseInt(pegrnaForm.edit_start, 10);
    const e = parseInt(pegrnaForm.edit_end, 10);
    if (!Number.isFinite(s) || !Number.isFinite(e) || s < 1 || e > seqLen || s > e) {
      setPegrnaError(`Edit range must be 1-indexed within 1..${seqLen} with start <= end.`);
      return;
    }
    const n_results = parseInt(pegrnaForm.n_results, 10) || 3;
    const alt = (pegrnaForm.alt || "").trim().toUpperCase();
    if (pegrnaForm.edit_type === "substitution" && alt.length !== (e - s + 1)) {
      setPegrnaError(`Substitution requires alt length (${alt.length}) to match selected range (${e - s + 1}).`);
      return;
    }
    if (pegrnaForm.edit_type === "insertion" && alt.length === 0) {
      setPegrnaError("Insertion requires a non-empty alt sequence.");
      return;
    }
    setPegrnaLoading(true);
    setPegrnaError("");
    setPegrnaResult(null);
    try {
      const res = await fetch("/api/plannotate/design-pegrnas", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          sequence,
          edit_start: s,
          edit_end: e,
          alt,
          edit_type: pegrnaForm.edit_type,
          n_results,
          use_pe3: pegrnaForm.use_pe3,
        }),
      });
      const data = await res.json();
      if (!res.ok || data.ok === false) {
        setPegrnaError(data.error || data.detail || "pegRNA design failed");
        setPegrnaLoading(false);
        return;
      }
      const pegs = (data.pegrnas || []) as any[];
      const colorForRank = (rank: number) => {
        if (rank === 1) return "#9bd4a8";   // top — green
        if (rank === 2) return "#cfd96e";   // mid — yellow-green
        return "#d99a6c";                   // 3+ — orange
      };
      const ngrnaSeen = new Set<string>();
      const ngrnaAnns: SeqVizAnnotation[] = [];
      for (const p of pegs) {
        if (!p.ngrna) continue;
        const ngK = `${p.ngrna.spacer}|${p.ngrna.start}|${p.ngrna.end}|${p.ngrna.strand}`;
        if (ngrnaSeen.has(ngK)) continue;
        ngrnaSeen.add(ngK);
        const ngDir = (p.ngrna.strand === "-" ? -1 : 1) as 1 | -1;
        ngrnaAnns.push({
          name: `ngRNA_${p.rank}_${(p.ngrna.spacer || "").slice(0, 6)}`,
          start: p.ngrna.start,
          end: p.ngrna.end,
          direction: ngDir,
          strand: ngDir,
          color: p.is_pe3b ? "#c97abd" : "#7d9dd0",
          type: "ngRNA",
          layer: "feature" as const,
          description: `ngRNA · ${p.is_pe3b ? "PE3b" : "PE3"} · nick-to-peg ${p.ngrna.nick_to_pegRNA} bp · cas9 ${p.ngrna.cas9_score}`,
          source: "pegrna_design",
          metadata: {
            ngrna_design: true,
            paired_pegrna_name: p.name,
            paired_pegrna_rank: p.rank,
            spacer: p.ngrna.spacer,
            original_spacer: p.ngrna.original_spacer,
            pam: p.ngrna.pam,
            cas9_score: p.ngrna.cas9_score,
            nick_to_pegRNA: p.ngrna.nick_to_pegRNA,
            is_pe3b: p.is_pe3b,
            edit_type: p.edit_type,
            edit_ref: p.edit_ref,
            edit_alt: p.edit_alt,
            edit_start_1based: p.edit_start_1based,
            edit_end_1based: p.edit_end_1based,
          },
        });
      }
      const newAnns: SeqVizAnnotation[] = pegs.map((p: any) => ({
        name: p.name,
        start: p.spacer_start,
        end: p.spacer_end,
        direction: (p.direction === -1 ? -1 : 1) as 1 | -1,
        strand: (p.direction === -1 ? -1 : 1) as 1 | -1,
        color: colorForRank(p.rank),
        type: "pegRNA",
        layer: "feature" as const,
        description: `pegRNA #${p.rank} · eff ${p.predicted_efficiency} · ${p.edit_type}${p.is_dpam ? " · dPAM" : ""}${p.is_pe3b ? " · PE3b" : ""}`,
        source: "pegrna_design",
        metadata: {
          pegrna_design: true,
          rank: p.rank,
          predicted_efficiency: p.predicted_efficiency,
          spacer: p.spacer,
          pam: p.pam,
          cas9_score: p.cas9_score,
          rtt: p.rtt,
          rtt_length: p.rtt_length,
          rtt_gc: p.rtt_gc,
          pbs: p.pbs,
          pbs_length: p.pbs_length,
          pbs_gc: p.pbs_gc,
          scaffold: p.scaffold,
          full_pegrna: p.full_pegrna,
          full_pegrna_length: p.full_pegrna_length,
          is_dpam: p.is_dpam,
          is_pe3b: p.is_pe3b,
          ngrna: p.ngrna,
          edit_type: p.edit_type,
          edit_ref: p.edit_ref,
          edit_alt: p.edit_alt,
          edit_start_1based: p.edit_start_1based,
          edit_end_1based: p.edit_end_1based,
          score_components: p.score_components,
          strand: p.strand,
        },
      }));
      const key = (a: SeqVizAnnotation) =>
        `${a.name}\u0001${a.start}\u0001${a.end}\u0001${a.direction ?? 0}\u0001${a.source ?? ""}`;
      const existingKeys = new Set(annotations.map(key));
      const allNew: SeqVizAnnotation[] = [...newAnns, ...ngrnaAnns];
      const additions = allNew.filter((a) => !existingKeys.has(key(a)));
      const nAdded = additions.length;
      setAnnotations((prev) => {
        const seen = new Set<string>();
        const out: SeqVizAnnotation[] = [];
        for (const a of [...prev, ...allNew]) {
          const k = key(a);
          if (seen.has(k)) continue;
          seen.add(k);
          out.push(a);
        }
        return out;
      });
      setPegrnaResult({ summary: data.summary, n_added: nAdded });
    } catch (e: any) {
      setPegrnaError(e.message || "Network error");
    } finally {
      setPegrnaLoading(false);
    }
  }, [sequence, pegrnaForm, annotations]);

  // Copy-to-clipboard helper for primer / amplicon strings in the modal.
  const handleCopyToClipboard = useCallback(async (s: string) => {
    try {
      await navigator.clipboard.writeText(s);
    } catch {
      // Fallback: selection-based copy
      const ta = document.createElement("textarea");
      ta.value = s;
      ta.style.position = "fixed";
      ta.style.left = "-9999px";
      document.body.appendChild(ta);
      ta.select();
      try { document.execCommand("copy"); } catch {}
      document.body.removeChild(ta);
    }
  }, []);

  // In-place reverse-complement of the plasmid: flip the sequence and
  // recompute annotation coordinates (start/end mirror across L) and
  // strand. Cloning features re-derive from sequence on next annotate so
  // we leave restrictionSites alone; if they go stale the user can
  // re-run Cloning Features.
  const handleRevComp = useCallback(() => {
    if (!sequence) return;
    const L = sequence.length;
    const COMP: Record<string, string> = {
      A: "T", T: "A", G: "C", C: "G", N: "N",
      a: "t", t: "a", g: "c", c: "g", n: "n",
      R: "Y", Y: "R", S: "S", W: "W", K: "M", M: "K",
      B: "V", V: "B", D: "H", H: "D",
    };
    const revcomp = (s: string) =>
      s.split("").reverse().map((b) => COMP[b] ?? b).join("");
    const newSeq = revcomp(sequence);
    setSequence(newSeq);
    setAnnotations((prev) =>
      prev.map((a) => {
        const flipDir = (d: any) => (d === 1 ? -1 : d === -1 ? 1 : d);
        const flippedStart = L - a.end;
        const flippedEnd = L - a.start;
        return {
          ...a,
          start: flippedStart,
          end: flippedEnd,
          direction: flipDir(a.direction),
          strand: flipDir((a as any).strand),
        };
      })
    );
    setSelection(null);
  }, [sequence]);

  // pLannotate + grammar annotator: annotate the current sequence and merge results
  const handlePlannotate = useCallback(async () => {
    if (!plannotateEndpoint || !sequence) return;
    setPlannotateLoading(true);
    setPlannotateError("");
    try {
      const res = await fetch(plannotateEndpoint, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ 
          sequence, 
          circular, 
          detailed: true, 
          hierarchical: true  // Pass pipeline selection
        }),
      });
      const data = await res.json();
      if (!data.ok) {
        setPlannotateError(data.error || "pLannotate failed");
        return;
      }

      // Note: GPT analysis is NOT automatically captured from annotation
      // Use the "Analyze Plasmid" button separately for GPT-powered analysis

      // Map regular annotations from pLannotate
      const newAnns: SeqVizAnnotation[] = (data.annotations || []).map((a: any) => ({
        name: a.name || "feature",
        start: a.start,
        end: a.end,
        direction: a.direction === -1 ? -1 : (a.direction === 0 ? 0 : 1),
        strand: a.direction === -1 ? -1 : (a.direction === 0 ? 0 : 1),
        color: a.color || "#7C3AED",
        type: a.type || "misc_feature",
        description: a.description || "",
        sseqid: a.sseqid,
        db: a.db,
        kb_data: a.kb_data || null,
        layer: "feature",
      }));

      // Map hierarchical annotations (modules, motifs, gaps, and cloning features)
      const hierarchicalAnns: SeqVizAnnotation[] = (data.hierarchical_annotations || []).map((a: any) => ({
        name: a.name || "annotation",
        start: a.start,
        end: a.end,
        direction: a.direction === -1 ? -1 : (a.direction === 0 ? 0 : 1),
        strand: a.direction === -1 ? -1 : (a.direction === 0 ? 0 : 1),
        color: a.color || "#6b5b95",
        type: a.type || a.motif_type || "misc_feature",
        description: a.description || "",
        layer: a.layer || "feature",
        motif_type: a.motif_type,
        module_type: a.module_type,
        source: a.source,
        payload_id: a.payload_id,
        metadata: a.metadata,
        sequence: a.sequence,
        feature_family: a.feature_family,
        subtype: a.subtype,
        cut_profile: a.cut_profile,
      }));

      // Combine all annotations
      const allNewAnns = [...newAnns, ...hierarchicalAnns];

      // Replace annotations with new results (clear previous)
      setAnnotations(allNewAnns);
      setCloningFeatures(data.cloning_features || null);
      setInteractions(data.interactions || []);
      if ((data.interactions || []).length > 0) {
        setInteractionDescLoading(true);
        fetch("/api/plannotate/describe_interactions", {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ interactions: data.interactions, plasmid_name: title }),
        })
          .then((r) => r.json())
          .then((j) => setInteractionDescription(j.markdown || j.summary || ""))
          .catch(() => setInteractionDescription(""))
          .finally(() => setInteractionDescLoading(false));
      } else {
        setInteractionDescription("");
      }
    } catch (e: any) {
      setPlannotateError(e.message || "Network error");
    } finally {
      setPlannotateLoading(false);
    }
  }, [plannotateEndpoint, sequence, circular]);


  // Heuristic annotation: POST sequence to heuristic endpoint
  const handleHeuristicAnnotate = useCallback(async () => {
    if (!plannotateEndpoint || !sequence) return;
    setHeuristicLoading(true);
    setPlannotateError("");
    try {
      // Use the heuristic endpoint instead of the standard one
      const heuristicEndpoint = "/api/plannotate-heuristic"; // Use dedicated route instead of: plannotateEndpoint.replace(
      const res = await fetch(heuristicEndpoint, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ sequence, circular, detailed: true }),
      });
      const data = await res.json();
      if (!data.ok) {
        setPlannotateError(data.error || "Heuristic annotation failed");
        return;
      }

      // Map regular annotations from pLannotate
      const newAnns: SeqVizAnnotation[] = (data.annotations || []).map((a: any) => ({
        name: a.name || "feature",
        start: a.start,
        end: a.end,
        direction: a.direction === -1 ? -1 : (a.direction === 0 ? 0 : 1),
        strand: a.direction === -1 ? -1 : (a.direction === 0 ? 0 : 1),
        color: a.color || "#7C3AED",
        type: a.type || "misc_feature",
        description: a.description || "",
        sseqid: a.sseqid,
        db: a.db,
        kb_data: a.kb_data || null,
        layer: "feature",
      }));

      // Map hierarchical annotations with heuristic scoring info
      const hierarchicalAnns: SeqVizAnnotation[] = (data.hierarchical_annotations || []).map((a: any) => ({
        name: a.name || "annotation",
        start: a.start,
        end: a.end,
        direction: a.direction === -1 ? -1 : (a.direction === 0 ? 0 : 1),
        strand: a.direction === -1 ? -1 : (a.direction === 0 ? 0 : 1),
        color: a.color || "#6b5b95",
        type: a.type || a.motif_type || "misc_feature",
        description: a.description || "",
        layer: a.layer || "module",
        motif_type: a.motif_type,
        module_type: a.module_type,
        score: a.score,
        confidence: a.confidence,
        rules_fired: a.rules_fired,
        payload_id: a.payload_id,
        metadata: a.metadata,
        feature_family: a.feature_family,
        subtype: a.subtype,
        cut_profile: a.cut_profile,
      }));

      // Combine all annotations
      const allNewAnns = [...newAnns, ...hierarchicalAnns];

      // Replace annotations with new results (clear previous)
      setAnnotations(allNewAnns);
      setCloningFeatures(data.cloning_features || null);
      setInteractions(data.interactions || []);
      if ((data.interactions || []).length > 0) {
        setInteractionDescLoading(true);
        fetch("/api/plannotate/describe_interactions", {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ interactions: data.interactions, plasmid_name: title }),
        })
          .then((r) => r.json())
          .then((j) => setInteractionDescription(j.markdown || j.summary || ""))
          .catch(() => setInteractionDescription(""))
          .finally(() => setInteractionDescLoading(false));
      } else {
        setInteractionDescription("");
      }

      console.log("[Heuristic] Summary:", data.summary);
    } catch (e: any) {
      setPlannotateError(e.message || "Network error");
    } finally {
      setHeuristicLoading(false);
    }
  }, [plannotateEndpoint, sequence, circular]);

  const handleLLMAnnotate = useCallback(async () => {
    if (!plannotateEndpoint || !sequence) return;
    setLlmLoading(true);
    setPlannotateError("");

    try {
      const res = await fetch("/api/plannotate-llm", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ sequence, circular, detailed: true }),
      });

      const data = await res.json();
      if (!data.ok) {
        setPlannotateError(data.error || "LLM annotation failed");
        return;
      }

      // Map regular features
      const newAnns = (data.annotations || []).map((a: any) => ({
        name: a.name || "feature",
        start: a.start,
        end: a.end,
        direction: a.direction === -1 ? -1 : (a.direction === 0 ? 0 : 1),
        strand: a.direction === -1 ? -1 : (a.direction === 0 ? 0 : 1),
        color: a.color || "#7C3AED",
        type: a.type || "misc_feature",
        layer: "feature",
        description: a.description || "",
        sseqid: a.sseqid,
        db: a.db,
        kb_data: a.kb_data || null,
      }));

      // Map hierarchical modules (including cloning-feature layer)
      const moduleAnns = (data.hierarchical_annotations || []).map((a: any) => ({
        name: a.name || "module",
        start: a.start,
        end: a.end,
        direction: a.direction === -1 ? -1 : (a.direction === 0 ? 0 : 1),
        strand: a.direction === -1 ? -1 : (a.direction === 0 ? 0 : 1),
        color: a.color,
        layer: a.layer || "module",
        module_type: a.module_type,
        metadata: a.metadata,
        payload_id: a.payload_id,
        feature_family: a.feature_family,
        subtype: a.subtype,
        cut_profile: a.cut_profile,
      }));

      setAnnotations([...newAnns, ...moduleAnns]);
      setCloningFeatures(data.cloning_features || null);
      setInteractions(data.interactions || []);
      if ((data.interactions || []).length > 0) {
        setInteractionDescLoading(true);
        fetch("/api/plannotate/describe_interactions", {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ interactions: data.interactions, plasmid_name: title }),
        })
          .then((r) => r.json())
          .then((j) => setInteractionDescription(j.markdown || j.summary || ""))
          .catch(() => setInteractionDescription(""))
          .finally(() => setInteractionDescLoading(false));
      } else {
        setInteractionDescription("");
      }
      console.log("[LLM] Summary:", data.summary);

    } catch (e: any) {
      setPlannotateError(e.message || "Network error");
    } finally {
      setLlmLoading(false);
    }
  }, [plannotateEndpoint, sequence, circular]);

  // Programmatic auto-annotate: when the parent passes autoAnnotateOnMount,
  // fire the same /api/plannotate-llm flow the viewer-internal Annotate button
  // uses, exactly once per loaded sequence.
  const autoAnnotateFiredRef = useRef<string | null>(null);
  useEffect(() => {
    if (!autoAnnotateOnMount) return;
    if (!sequence) return;
    if (llmLoading) return;
    if (autoAnnotateFiredRef.current === sequence) return;
    autoAnnotateFiredRef.current = sequence;
    handleLLMAnnotate();
  }, [autoAnnotateOnMount, sequence, llmLoading, handleLLMAnnotate]);

  // Apply only the cloning-feature scan (Step 2.75) — skips the full LLM +
  // pLannotate + rule pipeline so the user can annotate an uploaded plasmid
  // with restriction sites / Gateway att / PCR warnings directly.
  const handleScanCloningFeatures = useCallback(async () => {
    if (!sequence) return;
    setCloningScanLoading(true);
    setPlannotateError("");

    try {
      const res = await fetch("/api/cloning-scan", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ sequence }),
      });
      const data = await res.json();
      if (!data.ok) {
        setPlannotateError(data.error || "Cloning feature scan failed");
        return;
      }

      const cfAnns = (data.hierarchical_annotations || []).map((a: any) => ({
        name: a.name || "cloning_feature",
        start: a.start,
        end: a.end,
        direction: a.direction === -1 ? -1 : (a.direction === 0 ? 0 : 1),
        strand: a.direction === -1 ? -1 : (a.direction === 0 ? 0 : 1),
        color: a.color,
        layer: a.layer || "cloning_feature",
        feature_family: a.feature_family,
        subtype: a.subtype,
        cut_profile: a.cut_profile,
        metadata: a.metadata,
      }));

      // Replace existing cloning-feature annotations, preserve everything else.
      setAnnotations((prev) => [
        ...prev.filter((a: any) => a.layer !== "cloning_feature"),
        ...cfAnns,
      ]);
      setCloningFeatures(data.cloning_features || null);
      setShowCloningFeatures(true);
    } catch (e: any) {
      setPlannotateError(e.message || "Network error");
    } finally {
      setCloningScanLoading(false);
    }
  }, [sequence]);

  // Analyze plasmid: POST current annotations and sequence to analysis endpoint
  const handleAnalyzeIntent = useCallback(async () => {
    if (!analyzeIntentEndpoint || !annotations.length) return;
    setAnalyzeIntentLoading(true);
    setAnalyzeIntentError("");
    setIntentAnalysis(null);
    setModuleGraphData(null);

    // Use effectivePlasmidName if available, otherwise fall back to title
    const analysisTitle = effectivePlasmidName || title;

    try {
      const res = await fetch(analyzeIntentEndpoint, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          annotations: annotations.map((a) => ({
            name: a.name,
            start: a.start,
            end: a.end,
            role: (a as any).type || (a as any).role || null,
          })),
          sequence, // Include sequence for CDS analysis
          title: analysisTitle,
          sequence_length: sequence.length,
          circular,
        }),
      });
      const data = await res.json();
      if (!data.ok) {
        setAnalyzeIntentError(data.error || "Analysis failed");
        return;
      }

      // Store module graph for download
      if (data.module_graph) {
        setModuleGraphData(data.module_graph);
      }

      // If callback is provided, send to chat; otherwise show locally
      if (onAnalysisComplete) {
        onAnalysisComplete(data.analysis, data.module_graph);
      } else {
        setIntentAnalysis(data.analysis);
      }
    } catch (e: any) {
      setAnalyzeIntentError(e.message || "Network error");
    } finally {
      setAnalyzeIntentLoading(false);
    }
  }, [analyzeIntentEndpoint, annotations, title, effectivePlasmidName, sequence, circular, onAnalysisComplete]);

  // Download module graph as JSON
  const handleDownloadModuleGraph = useCallback(() => {
    if (!moduleGraphData) return;
    const filename = effectivePlasmidName
      ? `${effectivePlasmidName}_module_graph.json`
      : `${title.replace(/[^a-zA-Z0-9_-]/g, "_")}_module_graph.json`;
    const blob = new Blob([JSON.stringify(moduleGraphData, null, 2)], { type: "application/json" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    a.remove();
    URL.revokeObjectURL(url);
  }, [moduleGraphData, effectivePlasmidName, title]);

  // CSV-escape a single cell: wrap in double quotes if it contains a comma,
  // quote, or newline; double up internal quotes per RFC 4180.
  const csvCell = useCallback((v: any): string => {
    if (v === null || v === undefined) return "";
    const s = String(v);
    if (/[",\n\r]/.test(s)) return `"${s.replace(/"/g, '""')}"`;
    return s;
  }, []);

  const handleDownloadFeaturesCsv = useCallback(() => {
    const lines: string[] = [];

    lines.push("# FEATURES");
    lines.push([
      "name", "start_1based", "end_1based", "length", "direction",
      "type", "layer", "color", "source", "description",
      "kb_source", "kb_feature_class", "kb_subclass", "kb_polymerase",
      "kb_host_scope", "kb_descriptions", "sequence",
    ].join(","));
    for (const a of annotations) {
      const start1 = a.start + 1;
      const end1 = a.end;
      const length = a.end > a.start ? a.end - a.start : (sequence.length - a.start) + a.end;
      const dir = a.direction === -1 ? "reverse" : a.direction === 0 ? "none" : "forward";
      const featSeq = (() => {
        if (!sequence) return "";
        if (a.end >= a.start) return sequence.slice(a.start, a.end);
        return sequence.slice(a.start) + sequence.slice(0, a.end);
      })();
      const kb = a.kb_data || ({} as any);
      lines.push([
        csvCell(a.name),
        csvCell(start1),
        csvCell(end1),
        csvCell(length),
        csvCell(dir),
        csvCell(a.type || ""),
        csvCell(a.layer || ""),
        csvCell(a.color || ""),
        csvCell(a.source || ""),
        csvCell(a.description || ""),
        csvCell(kb.source_type || ""),
        csvCell(kb.feature_class || ""),
        csvCell(kb.subclass || ""),
        csvCell(kb.polymerase_class || ""),
        csvCell((kb.host_scope || []).join(";")),
        csvCell((kb.descriptions || []).join(" | ")),
        csvCell(featSeq),
      ].join(","));
    }

    const moduleAnns = annotations.filter((a) => a.layer === "module");
    if (moduleAnns.length > 0) {
      lines.push("");
      lines.push("# MODULES");
      lines.push([
        "name", "start_1based", "end_1based", "module_type", "module_family",
        "payload_id", "promoter_id", "polya_id", "marker_id", "host_scope",
        "canonical_ids",
      ].join(","));
      for (const a of moduleAnns) {
        const m = a.metadata || ({} as any);
        lines.push([
          csvCell(a.name),
          csvCell(a.start + 1),
          csvCell(a.end),
          csvCell(a.module_type || ""),
          csvCell(a.module_family || ""),
          csvCell(a.payload_id || ""),
          csvCell(m.promoter_id || ""),
          csvCell(m.polya_id || ""),
          csvCell(m.marker_id || ""),
          csvCell((m.host_scope || []).join(";")),
          csvCell((m.canonical_ids || []).join(";")),
        ].join(","));
      }
    }

    if (interactions.length > 0) {
      lines.push("");
      lines.push("# INTERACTIONS");
      lines.push([
        "interaction_id", "interaction_type", "sbo_term", "rule_id",
        "confidence", "participant_names", "participant_ranges_1based", "notes",
      ].join(","));
      for (const ix of interactions) {
        const names = (ix.participants || []).map((p) => p.name).join(";");
        const ranges = (ix.participants || [])
          .map((p) => (p.start != null && p.end != null ? `${p.start + 1}..${p.end}` : ""))
          .join(";");
        lines.push([
          csvCell(ix.interaction_id),
          csvCell(ix.interaction_type),
          csvCell(ix.sbo_term || ""),
          csvCell(ix.rule_id || ""),
          csvCell(ix.confidence != null ? ix.confidence : ""),
          csvCell(names),
          csvCell(ranges),
          csvCell(ix.notes || ""),
        ].join(","));
      }
    }

    const stem = effectivePlasmidName
      ? effectivePlasmidName
      : title.replace(/[^a-zA-Z0-9_-]/g, "_");
    const blob = new Blob([lines.join("\n") + "\n"], { type: "text/csv;charset=utf-8" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = `${stem}_features.csv`;
    document.body.appendChild(a);
    a.click();
    a.remove();
    URL.revokeObjectURL(url);
  }, [annotations, interactions, sequence, effectivePlasmidName, title, csvCell]);

  // Download a CSV with one row per designed primer. The required columns
  // (number, name, sequence 5'→3', amplicon name, amplicon sequence) come
  // first; remaining characteristics (Tm, dimer/hairpin/end Tm, mispriming,
  // Sanger score) trail to the right.
  const triggerCsvDownload = useCallback((rows: string[][], filename: string) => {
    const lines = rows.map((r) => r.map(csvCell).join(","));
    const blob = new Blob([lines.join("\n") + "\n"], { type: "text/csv;charset=utf-8" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    a.remove();
    URL.revokeObjectURL(url);
  }, [csvCell]);

  const handleDownloadPrimerCsv = useCallback(() => {
    if (!primerResult) return;
    const stem = (primerResult.amplicon_name || "primers").replace(/[^a-zA-Z0-9_-]/g, "_");
    const sangerByDir: Record<string, any> = {};
    for (const s of (primerResult.sanger_scores || [])) {
      sangerByDir[s.direction] = s;
    }
    const header = [
      "primer_number", "name", "sequence_5to3", "amplicon_name", "amplicon_sequence",
      "direction", "length_nt", "tm_annealing_c", "adapter_5prime", "annealing_portion",
      "product_size_bp", "amplicon_start_1based", "amplicon_end_1based",
      "pair_penalty", "application", "selection_method",
      "hairpin_tm_c", "homodimer_tm_c", "end_stability_tm_c",
      "n_mispriming_sites", "sanger_score", "sanger_rating",
    ];
    const rows: string[][] = [header];
    const primers = [
      {
        role: "forward",
        name: primerResult.fwd_name,
        seq: primerResult.left_primer,
        adapter: primerResult.left_adapter || "",
        anneal: primerResult.left_annealing,
        tm: primerResult.left_tm,
        scores: primerResult.left_scores || {},
        misprime: primerResult.left_mispriming_sites || [],
      },
      {
        role: "reverse",
        name: primerResult.rev_name,
        seq: primerResult.right_primer,
        adapter: primerResult.right_adapter || "",
        anneal: primerResult.right_annealing,
        tm: primerResult.right_tm,
        scores: primerResult.right_scores || {},
        misprime: primerResult.right_mispriming_sites || [],
      },
    ];
    primers.forEach((p, i) => {
      const sanger = sangerByDir[p.role];
      rows.push([
        String(i + 1),
        p.name || "",
        p.seq || "",
        primerResult.amplicon_name || "",
        primerResult.amplicon || "",
        p.role,
        p.seq ? String(p.seq.length) : "",
        typeof p.tm === "number" ? p.tm.toFixed(2) : "",
        p.adapter,
        p.anneal || "",
        primerResult.product_size != null ? String(primerResult.product_size) : "",
        primerResult.amplicon_start_1based != null ? String(primerResult.amplicon_start_1based) : "",
        primerResult.amplicon_end_1based != null ? String(primerResult.amplicon_end_1based) : "",
        typeof primerResult.pair_penalty === "number" ? primerResult.pair_penalty.toFixed(3) : "",
        primerResult.application || "",
        primerResult.selection_method || "",
        p.scores && typeof p.scores.hairpin_tm === "number" ? p.scores.hairpin_tm.toFixed(2) : "",
        p.scores && typeof p.scores.homodimer_tm === "number" ? p.scores.homodimer_tm.toFixed(2) : "",
        p.scores && typeof p.scores.end_stability_tm === "number" ? p.scores.end_stability_tm.toFixed(2) : "",
        Array.isArray(p.misprime) ? String(p.misprime.length) : "",
        sanger && typeof sanger.overall_score === "number" ? sanger.overall_score.toFixed(1) : "",
        sanger ? (sanger.rating || "") : "",
      ]);
    });
    triggerCsvDownload(rows, `${stem}_primers.csv`);
  }, [primerResult, triggerCsvDownload]);

  const handleDownloadGuideCsv = useCallback(() => {
    if (!guideResult || !guideResult.guides) return;
    const stem = `sgRNAs_${(guideResult.design_region_1based || "design").replace(/[^a-zA-Z0-9_-]/g, "_")}`;
    const header = [
      "guide_number", "name", "spacer_5to3", "context_name", "context_30mer",
      "pam", "direction", "start_1based", "end_1based",
      "score", "score_method", "gc_fraction", "max_homopolymer", "n_offtargets",
    ];
    const rows: string[][] = [header];
    (guideResult.guides as any[]).forEach((g, i) => {
      const dirStr = g.direction === -1 ? "reverse" : "forward";
      rows.push([
        String(i + 1),
        g.name || "",
        g.spacer || "",
        `${g.name || "sgRNA_" + (i + 1)}_context`,
        g.context_30mer || "",
        g.pam || "",
        dirStr,
        g.start != null ? String(g.start + 1) : "",
        g.end != null ? String(g.end) : "",
        typeof g.score === "number" ? g.score.toFixed(2) : (g.score ?? ""),
        g.score_method || guideResult.summary?.score_method || "",
        typeof g.gc_fraction === "number" ? g.gc_fraction.toFixed(3) : "",
        g.max_homopolymer != null ? String(g.max_homopolymer) : "",
        g.n_offtargets != null ? String(g.n_offtargets) : "",
      ]);
    });
    triggerCsvDownload(rows, `${stem}.csv`);
  }, [guideResult, triggerCsvDownload]);

  // Determine if we have a range selection
  const hasRangeSelection = selection && selection.start !== selection.end;
  const selectionStart = hasRangeSelection ? Math.min(selection!.start, selection!.end) : null;
  const selectionEnd = hasRangeSelection ? Math.max(selection!.start, selection!.end) : null;

  return (
    <div
      ref={containerRef}
      className="mt-6 p-4"
      style={{
        background: "transparent",
        marginLeft: "calc(-50vw + 50%)",
        marginRight: "calc(-50vw + 50%)",
        width: "100vw",
        maxWidth: "100vw",
        paddingLeft: "2rem",
        paddingRight: "2rem",
        fontFamily: "var(--font-body)",
        color: "var(--mint-200)",
        borderTop: "1px solid rgba(219, 239, 231, 0.10)",
      }}
    >
      {/* Inject custom styles */}
      <style>{viewerStyles}</style>

      {/* Hidden file input for GenBank upload */}
      <input
        ref={fileInputRef}
        type="file"
        accept=".gb,.gbk,.genbank"
        onChange={handleFileUpload}
        style={{ display: "none" }}
      />

      {/* Toolbar */}
      <div className="mb-3 flex flex-col gap-2">
        <div
          className="splicify-display"
          style={{
            color: "var(--mint-200)",
            fontSize: 26,
            fontWeight: 600,
            letterSpacing: "-0.015em",
            textAlign: "center",
            margin: "4px 0 2px",
            display: effectivePlasmidName ? "block" : "none",
          }}
        >
          {effectivePlasmidName || ""}
        </div>
        <div className="splicify-toolbar">
          {/* File */}
          <div className="splicify-tb-group">
            <button onClick={handleUploadClick} className="splicify-tb-btn is-light" title="Upload a GenBank file">
              Upload .gb
            </button>
          </div>
          <div className="splicify-tb-sep" />

          {/* Annotate / Scan */}
          <div className="splicify-tb-group">
            {plannotateEndpoint && (
              <button
                onClick={handleLLMAnnotate}
                disabled={llmLoading || !sequence}
                className="splicify-tb-btn is-purple"
              >
                {llmLoading ? "Annotating…" : "Annotate"}
              </button>
            )}
            <div className="splicify-tb-sep" />
            <button
              onClick={handleScanCloningFeatures}
              disabled={cloningScanLoading || !sequence}
              className="splicify-tb-btn is-magenta"
              title="Apply cloning-feature annotations (restriction sites, Gateway att, PCR warnings) without running the full Annotate pipeline"
            >
              {cloningScanLoading ? "Scanning…" : "Scan Cloning Features"}
            </button>
            <div className="splicify-tb-sep" />
          <div style={{ position: "relative" }}>
            <button
              onClick={() => setShowAnnotationsMenu((v) => !v)}
              className={`splicify-tb-btn${(showModules && showInteractionChords) ? " is-active" : ""}`}
              title="Toggle module annotations and interaction chords"
            >
              Annotations ▼
            </button>
            {showAnnotationsMenu && (
              <div style={{
                position: "absolute",
                right: 0,
                top: "100%",
                zIndex: 50,
                background: "#1a2f25",
                border: "1px solid #2d4a3e",
                borderRadius: 6,
                minWidth: 220,
                marginTop: 4,
                padding: "8px 0",
                color: "#dbefe7",
                fontSize: "0.85rem",
              }}>
                <label style={{ display: "flex", alignItems: "center", gap: 8, padding: "4px 14px", cursor: "pointer" }}>
                  <input
                    type="checkbox"
                    checked={showModules}
                    onChange={(e) => setShowModules(e.target.checked)}
                  />
                  <span style={{ fontWeight: 600 }}>Show modules</span>
                </label>
                <label style={{ display: "flex", alignItems: "center", gap: 8, padding: "4px 14px", cursor: "pointer", opacity: interactions.length > 0 ? 1 : 0.5 }}>
                  <input
                    type="checkbox"
                    checked={showInteractionChords}
                    onChange={(e) => setShowInteractionChords(e.target.checked)}
                    disabled={interactions.length === 0}
                  />
                  <span style={{ fontWeight: 600 }}>Show interaction chords</span>
                </label>
              </div>
            )}
          </div>
          {cloningFeatures && (
            <div style={{ position: "relative" }}>
              <button
                onClick={() => setShowCloningMenu((v) => !v)}
                className={`splicify-tb-btn${showCloningFeatures ? " is-active" : ""}`}
                title="Toggle cloning-feature annotations (hidden by default)"
              >
                Cloning Features {showCloningFeatures ? "●" : "▼"}
              </button>
              {showCloningMenu && (
                <div style={{
                  position: "absolute",
                  right: 0,
                  top: "100%",
                  zIndex: 50,
                  background: "#1a2f25",
                  border: "1px solid #2d4a3e",
                  borderRadius: 6,
                  minWidth: 260,
                  marginTop: 4,
                  padding: "8px 0",
                  color: "#dbefe7",
                  fontSize: "0.85rem",
                }}>
                  <label style={{ display: "flex", alignItems: "center", gap: 8, padding: "4px 14px", cursor: "pointer" }}>
                    <input
                      type="checkbox"
                      checked={showCloningFeatures}
                      onChange={(e) => setShowCloningFeatures(e.target.checked)}
                    />
                    <span style={{ fontWeight: 600 }}>Show cloning features</span>
                  </label>
                  <div style={{ borderTop: "1px solid #2d4a3e", margin: "6px 0" }} />
                  <div style={{ padding: "2px 14px", opacity: 0.75, fontSize: "0.72rem", textTransform: "uppercase" }}>Categories</div>
                  {([
                    ["Type II enzymes", showCloningReII, setShowCloningReII],
                    ["Type IIs enzymes", showCloningReIIs, setShowCloningReIIs],
                    ["Gateway att sites", showCloningGateway, setShowCloningGateway],
                    ["PCR design warnings", showCloningPcr, setShowCloningPcr],
                  ] as [string, boolean, (v: boolean) => void][]).map(([label, val, setter]) => (
                    <label key={label} style={{ display: "flex", alignItems: "center", gap: 8, padding: "3px 14px", cursor: "pointer", opacity: showCloningFeatures ? 1 : 0.5 }}>
                      <input type="checkbox" checked={val} onChange={(e) => setter(e.target.checked)} disabled={!showCloningFeatures} />
                      <span>{label}</span>
                    </label>
                  ))}
                  <div style={{ borderTop: "1px solid #2d4a3e", margin: "6px 0" }} />
                  <div style={{ padding: "2px 14px", opacity: 0.75, fontSize: "0.72rem", textTransform: "uppercase" }}>Restriction cutter filter</div>
                  {([
                    ["unique", "Unique cutters"],
                    ["2", "2-cutters"],
                    ["3", "3-cutters"],
                    ["all", "All cutters"],
                    ["none", "None (show non-cutters list)"],
                  ] as [CutterCountFilter, string][]).map(([val, label]) => (
                    <label key={val} style={{ display: "flex", alignItems: "center", gap: 8, padding: "3px 14px", cursor: "pointer", opacity: showCloningFeatures ? 1 : 0.5 }}>
                      <input
                        type="radio"
                        name="cutter-filter"
                        checked={cutterFilter === val}
                        onChange={() => { setCutterFilter(val); if (val === "none") setShowNonCutters(true); }}
                        disabled={!showCloningFeatures}
                      />
                      <span>{label}</span>
                    </label>
                  ))}
                  {cutterFilter === "none" && (
                    <button
                      onClick={() => setShowNonCutters((v) => !v)}
                      style={{
                        margin: "6px 14px",
                        padding: "4px 8px",
                        border: "1px solid #2d4a3e",
                        borderRadius: 4,
                        background: "transparent",
                        color: "#dbefe7",
                        fontSize: "0.75rem",
                        cursor: "pointer",
                      }}
                    >
                      {showNonCutters ? "Hide" : "Show"} non-cutters ({nonCuttersList.length})
                    </button>
                  )}
                </div>
              )}
            </div>
          )}
          {restrictionSites && restrictionSites.length > 0 && (
            <div style={{ position: "relative" }}>
              <button
                onClick={() => setShowReSiteMenu((v) => !v)}
                className={`splicify-tb-btn${reSiteFilter !== "none" ? " is-active" : ""}`}
              >
                RE Sites {reSiteFilter !== "none" ? "●" : "▼"}
              </button>
              {showReSiteMenu && (
                <div style={{
                  position: "absolute",
                  right: 0,
                  top: "100%",
                  zIndex: 50,
                  background: "#1a2f25",
                  border: "1px solid #2d4a3e",
                  borderRadius: 6,
                  minWidth: 180,
                  marginTop: 4,
                }}>
                  {([
                    ["none",      "Hide All"],
                    ["type2s_re", "Type IIs (Assembly)"],
                    ["type2_re",  "Type II (Common)"],
                    ["gateway",   "Gateway Sites"],
                    ["cre_lox",   "Cre-Lox / FRT"],
                    ["all",       "Show All"],
                  ] as [ReSiteFilter, string][]).map(([value, label]) => (
                    <button
                      key={value}
                      onClick={() => { setReSiteFilter(value); setShowReSiteMenu(false); }}
                      style={{
                        background: reSiteFilter === value ? "#2d4a3e" : "transparent",
                        color: "#dbefe7",
                        display: "block",
                        width: "100%",
                        textAlign: "left",
                        padding: "6px 14px",
                        border: "none",
                        cursor: "pointer",
                        fontSize: "0.875rem",
                      }}
                    >
                      {label}
                    </button>
                  ))}
                </div>
              )}
            </div>
          )}
          </div>
          {analyzeIntentEndpoint && (
            <>
              <div className="splicify-tb-sep" />
              <div className="splicify-tb-group">
                <button
                  onClick={handleAnalyzeIntent}
                  disabled={analyzeIntentLoading || !annotations.length}
                  className="splicify-tb-btn is-outline"
                >
                  {analyzeIntentLoading ? "Analyzing…" : "Analyze Plasmid"}
                </button>
              </div>
            </>
          )}
          <div className="splicify-tb-sep" />

          {/* Edit */}
          <div className="splicify-tb-group">
            <button onClick={handleAddAnnotation} className="splicify-tb-btn">
              Add Annotation{hasRangeSelection ? " (Selected)" : ""}
            </button>
            <div className="splicify-tb-sep" />
            <button
              onClick={handleOpenTranslationModal}
              className="splicify-tb-btn"
              disabled={!hasRangeSelection}
              title={hasRangeSelection
                ? "Translate the selected region in a chosen reading frame and add it as a clickable per-AA annotation"
                : "Select a region first"}
            >
              Add Translation{hasRangeSelection ? " (Selected)" : ""}
            </button>
            <div className="splicify-tb-sep" />
            <button onClick={handleOpenAddChangeModal} className="splicify-tb-btn">
              Add / Change Sequence
            </button>
          </div>
          <div className="splicify-tb-sep" />

          {/* Import + Design */}
          <div className="splicify-tb-group">
            {(() => {
              const hasCsvImports = annotations.some((a) => a.source === "csv_import");
              return (
                <button
                  onClick={handleOpenImportModal}
                  disabled={!sequence}
                  className={`splicify-tb-btn${hasCsvImports ? " is-used" : ""}`}
                  title={
                    hasCsvImports
                      ? "CSV imports already on this plasmid — re-importing the same file will be deduped, but importing a different CSV is fine"
                      : "Upload a CSV of (name, sequence) rows to map onto this plasmid by sequence identity"
                  }
                >
                  Import Annotation(s){hasCsvImports ? " ✓" : ""}
                </button>
              );
            })()}
            <div className="splicify-tb-sep" />
            <button
              onClick={handleOpenPrimerModal}
              disabled={!sequence}
              className="splicify-tb-btn"
              title="Design Primer3 primers around the selected (or specified) region with optional exclusion"
            >
              Design Primers{hasRangeSelection ? " (Selected)" : ""}
            </button>
            <div className="splicify-tb-sep" />
            <button
              onClick={handleOpenGuideModal}
              disabled={!sequence}
              className="splicify-tb-btn"
              title="Scan a region for CRISPR sgRNA candidates ranked by efficiency score"
            >
              Design Guides{hasRangeSelection ? " (Selected)" : ""}
            </button>
            <div className="splicify-tb-sep" />
            <button
              onClick={handleOpenPegrnaModal}
              disabled={!sequence}
              className="splicify-tb-btn"
              title="Design prime-editing pegRNAs (easy_prime PE3 XGBoost model) — pick a selection then specify the edit"
            >
              Design pegRNA{hasRangeSelection ? " (Selected)" : ""}
            </button>
            <div className="splicify-tb-sep" />
            <button
              onClick={handleRevComp}
              disabled={!sequence}
              className="splicify-tb-btn"
              title="Reverse-complement the plasmid in place: flips sequence and every annotation's strand + coordinates"
            >
              Rev Comp
            </button>
          </div>
          <div className="splicify-tb-sep" />

          {/* Export */}
          <div className="splicify-tb-group">
            <GenBankDownloadButton
              sequence={sequence}
              features={annotations.map((a) => ({
                start: a.start,
                end: a.end,
                name: a.name,
                strand: a.direction === -1 ? -1 : a.direction === 1 ? 1 : undefined,
              }))}
              filename={effectivePlasmidName ? `${effectivePlasmidName}_aip.gb` : `${title.replace(/[^a-zA-Z0-9_-]/g, "_")}_aip.gb`}
              locusName={(effectivePlasmidName || title).slice(0, 16).replace(/[^a-zA-Z0-9_]/g, "_") || "PLASMID"}
              definition={effectivePlasmidName || title}
              circular={circular}
              className="splicify-tb-btn is-light"
            />
            <div className="splicify-tb-sep" />
            <button
              onClick={handleDownloadFeaturesCsv}
              disabled={!annotations.length && !interactions.length}
              className="splicify-tb-btn is-light"
              title="Download all features, modules, and interactions as a single CSV"
            >
              Download Feature List
            </button>
            {moduleGraphData && (
              <>
                <div className="splicify-tb-sep" />
                <button
                  onClick={handleDownloadModuleGraph}
                  className="splicify-tb-btn"
                >
                  Download Module Graph
                </button>
              </>
            )}
          </div>
        </div>
      </div>

      {/* Annotation error */}
      {plannotateError && (
        <div className="mb-2 text-sm px-3 py-1 rounded" style={{ backgroundColor: "rgba(239,68,68,0.2)", color: "#fca5a5" }}>
          Annotation error: {plannotateError}
        </div>
      )}

      {/* Sequence info and selection status */}
      <div className="mb-2 text-sm flex items-center gap-4 flex-wrap" style={{ color: "#dbefe7" }}>
        <span>Length: {sequence.length} bp</span>
        <span>Annotations: {annotations.length}</span>
        {hasRangeSelection && (
          <span className="px-2 py-1 rounded" style={{ backgroundColor: "#2d4a3e" }}>
            Selected: {selectionStart! + 1}..{selectionEnd}
            {" "}({selectionEnd! - selectionStart!} bp)
          </span>
        )}
      </div>

      {/* Help text */}
      <div className="mb-2 text-xs" style={{ color: "#dbefe7", opacity: 0.7 }}>
        Tip: Shift+click to extend a selection. ⌘C copy, ⌘V paste, ⌘Z undo, ⌘⇧Z redo, ⌘F find. Delete removes selected bp.
      </div>

      {/* Find panel (⌘F / Ctrl+F) — searches forward + reverse-complement on the
          current sequence. Setting findIdx scrolls the viewers to the match. */}
      {findOpen && (
        <div
          className="mb-2 rounded-lg flex items-center gap-2 p-2"
          style={{ backgroundColor: "#1a2f25", color: "#dbefe7", border: "1px solid #2d4a3e" }}
        >
          <span style={{ fontSize: "0.8rem", fontWeight: 600 }}>Find</span>
          <input
            ref={findInputRef}
            type="text"
            value={findQuery}
            onChange={(e) => {
              setFindQuery(e.target.value);
              setFindIdx(0);
            }}
            onKeyDown={(e) => {
              if (e.key === "Enter") {
                e.preventDefault();
                if (findMatches.length > 0) {
                  setFindIdx((i) => (i + (e.shiftKey ? -1 : 1) + findMatches.length) % findMatches.length);
                }
              } else if (e.key === "Escape") {
                e.preventDefault();
                setFindOpen(false);
              }
            }}
            placeholder="Type a DNA sequence (≥2 bp, ambiguity codes ok)"
            spellCheck={false}
            autoCorrect="off"
            autoCapitalize="off"
            style={{
              flex: 1,
              minWidth: 0,
              padding: "4px 8px",
              borderRadius: 4,
              border: "1px solid #46896c",
              background: "#0f1f17",
              color: "#dbefe7",
              fontFamily: "monospace",
              fontSize: "0.85rem",
              letterSpacing: "0.05em",
            }}
          />
          <span style={{ fontSize: "0.75rem", opacity: 0.85, minWidth: 70, textAlign: "center" }}>
            {findMatches.length === 0 && findQuery
              ? "No matches"
              : findMatches.length === 0
                ? ""
                : `${(findIdx % findMatches.length + findMatches.length) % findMatches.length + 1} / ${findMatches.length}`}
          </span>
          <button
            type="button"
            disabled={findMatches.length === 0}
            onClick={() =>
              setFindIdx((i) => (i - 1 + Math.max(1, findMatches.length)) % Math.max(1, findMatches.length))
            }
            className="px-2 py-1 rounded text-xs"
            style={{ background: "#2d4a3e", color: "#dbefe7", opacity: findMatches.length === 0 ? 0.4 : 1 }}
          >
            Prev
          </button>
          <button
            type="button"
            disabled={findMatches.length === 0}
            onClick={() => setFindIdx((i) => (i + 1) % Math.max(1, findMatches.length))}
            className="px-2 py-1 rounded text-xs"
            style={{ background: "#2d4a3e", color: "#dbefe7", opacity: findMatches.length === 0 ? 0.4 : 1 }}
          >
            Next
          </button>
          <button
            type="button"
            onClick={() => setFindOpen(false)}
            className="px-2 py-1 rounded text-xs"
            style={{ background: "transparent", color: "#dbefe7", border: "1px solid #46896c" }}
          >
            Close
          </button>
        </div>
      )}

      {/* Non-cutters list (shown when cutter filter = "none") */}
      {showCloningFeatures && cutterFilter === "none" && showNonCutters && cloningFeatures && (
        <div
          className="mb-2 p-3 rounded"
          style={{ backgroundColor: "#1a2f25", color: "#dbefe7", border: "1px solid #2d4a3e" }}
        >
          <div style={{ fontWeight: 600, marginBottom: 6, fontSize: "0.85rem" }}>
            Non-cutting enzymes ({nonCuttersList.length})
          </div>
          <div style={{ display: "flex", flexWrap: "wrap", gap: 6 }}>
            {nonCuttersList.map((name) => (
              <span
                key={name}
                style={{
                  padding: "3px 8px",
                  borderRadius: 4,
                  background: "#2d4a3e",
                  fontSize: "0.75rem",
                  fontFamily: "monospace",
                }}
              >
                {name}
              </span>
            ))}
            {nonCuttersList.length === 0 && (
              <span style={{ opacity: 0.7, fontSize: "0.8rem" }}>Every enabled enzyme cuts this template at least once.</span>
            )}
          </div>
        </div>
      )}

      {/* Main content: top = viewers (circular + linear), bottom = panels */}
      <div className="flex flex-col gap-3">
      <div
        ref={viewerOuterRef}
        className="rounded-2xl overflow-hidden flex relative"
        style={{
          backgroundColor: "#ffffff",
          ...(viewerLayout === "split"
            ? { height }
            : { width: "100%", aspectRatio: "1 / 1" }),
        }}
      >
        {/* View switcher: left group toggles Circular vs Flat for the left
            pane; right group toggles layout (full circular, split, full linear).
            Hoisted to the outer container so it remains visible regardless of
            which pane is rendered. */}
        <div
          style={{
            position: "absolute", top: 8, left: 8, zIndex: 10,
            display: "flex", gap: 4, alignItems: "center",
            background: "#46896c", borderRadius: 6, padding: 2,
            opacity: 0.9,
          }}
        >
          <button
            type="button"
            onClick={() => setLeftPaneView("circular")}
            aria-pressed={leftPaneView === "circular"}
            aria-label="Circular view"
            title="Circular view"
            style={{
              width: 26, height: 22,
              display: "inline-flex", alignItems: "center", justifyContent: "center",
              background: leftPaneView === "circular" ? "#dbefe7" : "transparent",
              color: leftPaneView === "circular" ? "#105b39" : "#dbefe7",
              border: "none", borderRadius: 4, cursor: "pointer", padding: 0,
            }}
          >
            <svg width="14" height="14" viewBox="0 0 16 16" fill="none" aria-hidden="true">
              <circle cx="8" cy="8" r="5.5" stroke="currentColor" strokeWidth="1.75" />
            </svg>
          </button>
          <button
            type="button"
            onClick={() => setLeftPaneView("flat")}
            aria-pressed={leftPaneView === "flat"}
            aria-label="Flat view"
            title="Flat view"
            style={{
              width: 26, height: 22,
              display: "inline-flex", alignItems: "center", justifyContent: "center",
              background: leftPaneView === "flat" ? "#dbefe7" : "transparent",
              color: leftPaneView === "flat" ? "#105b39" : "#dbefe7",
              border: "none", borderRadius: 4, cursor: "pointer", padding: 0,
            }}
          >
            <svg width="14" height="14" viewBox="0 0 16 16" fill="none" aria-hidden="true">
              <line x1="2" y1="8" x2="14" y2="8" stroke="currentColor" strokeWidth="1.75" strokeLinecap="round" />
            </svg>
          </button>
          <div style={{ width: 1, height: 16, background: "#9ca3af", margin: "0 2px" }} />
          <button
            type="button"
            onClick={() => setViewerLayout("circular")}
            aria-pressed={viewerLayout === "circular"}
            aria-label="Full circular/flat view"
            title={leftPaneView === "flat" ? "Full flat view" : "Full circular view"}
            style={{
              width: 26, height: 22,
              display: "inline-flex", alignItems: "center", justifyContent: "center",
              background: viewerLayout === "circular" ? "#dbefe7" : "transparent",
              color: viewerLayout === "circular" ? "#105b39" : "#dbefe7",
              border: "none", borderRadius: 4, cursor: "pointer", padding: 0,
            }}
          >
            {leftPaneView === "flat" ? (
              <svg width="18" height="14" viewBox="0 0 20 16" fill="none" aria-hidden="true">
                <line x1="2" y1="8" x2="18" y2="8" stroke="currentColor" strokeWidth="2.25" strokeLinecap="round" />
              </svg>
            ) : (
              <svg width="18" height="18" viewBox="0 0 20 20" fill="none" aria-hidden="true">
                <circle cx="10" cy="10" r="7.5" stroke="currentColor" strokeWidth="2" />
              </svg>
            )}
          </button>
          <button
            type="button"
            onClick={() => setViewerLayout("split")}
            aria-pressed={viewerLayout === "split"}
            aria-label="Split view"
            title="Split view (circle and line)"
            style={{
              width: 30, height: 22,
              display: "inline-flex", alignItems: "center", justifyContent: "center",
              background: viewerLayout === "split" ? "#dbefe7" : "transparent",
              color: viewerLayout === "split" ? "#105b39" : "#dbefe7",
              border: "none", borderRadius: 4, cursor: "pointer", padding: 0,
            }}
          >
            <svg width="22" height="14" viewBox="0 0 24 16" fill="none" aria-hidden="true">
              <circle cx="6" cy="8" r="4.5" stroke="currentColor" strokeWidth="1.75" />
              <line x1="13" y1="8" x2="22" y2="8" stroke="currentColor" strokeWidth="1.75" strokeLinecap="round" />
            </svg>
          </button>
          <button
            type="button"
            onClick={() => setViewerLayout("linear")}
            aria-pressed={viewerLayout === "linear"}
            aria-label="Full linear view"
            title="Full linear view"
            style={{
              width: 30, height: 22,
              display: "inline-flex", alignItems: "center", justifyContent: "center",
              background: viewerLayout === "linear" ? "#dbefe7" : "transparent",
              color: viewerLayout === "linear" ? "#105b39" : "#dbefe7",
              border: "none", borderRadius: 4, cursor: "pointer", padding: 0,
            }}
          >
            <svg width="24" height="14" viewBox="0 0 26 16" fill="none" aria-hidden="true">
              <line x1="2" y1="8" x2="24" y2="8" stroke="currentColor" strokeWidth="2.25" strokeLinecap="round" />
            </svg>
          </button>
        </div>

        {/* Left pane: Circular or Flat (hidden when layout=linear) */}
        {viewerLayout !== "linear" && (
        <div
          ref={viewerRef}
          className="p-3 overflow-hidden seqviz-viewer-container flex-1 relative"
          style={{
            backgroundColor: useCustomViewer ? "transparent" : "#ffffff",
            borderRight: viewerLayout === "split" ? "1px solid #46896c" : "none",
          }}
        >
          <div className="h-full w-full">
            {leftPaneView === "flat" ? (
              <FlatPlasmidViewer
                sequence={sequence}
                annotations={allAnnotations as SeqVizAnnotation[]}
                height={viewerOuterHeight - 24}
                selection={selection ? { start: selection.start, end: selection.end } : null}
                onSelectionChange={(sel) => handleSelectionNoRecenter(sel ? { start: sel.start, end: sel.end, clockwise: true } : null)}
                onAnnotationSelect={(sel) => handleSelection({ start: sel.start, end: sel.end, clockwise: true })}
                onAnnotationClick={(ann, idx) => openEditAnnotationModal(idx)}
              />
            ) : useCustomViewer ? (
              <CircularPlasmidViewer
                sequence={sequence}
                annotations={allAnnotations as SeqVizAnnotation[]}
                title={title}
                height={viewerOuterHeight - 24}
                onAnnotationClick={(ann, idx) => openEditAnnotationModal(idx)}
                onSelectionChange={(sel) => handleSelection(sel as Selection | null)}
                selection={selection ? { start: selection.start, end: selection.end } : null}
                interactions={interactions}
                showInteractions={showInteractionChords}
                highlightedInteractionId={highlightedInteractionId}
              />
            ) : (
              <SeqVizAny
                name={title}
                seq={sequence}
                sequence={sequence}
                annotations={allAnnotations}
                viewer={viewerMode}
                onSelection={handleSelection}
                style={{ height: "100%", width: "100%" }}
              />
            )}
          </div>
        </div>
        )}
        {/* Linear pane (hidden when layout=circular) */}
        {viewerLayout !== "circular" && (
        <div className="overflow-hidden flex-1" style={{ backgroundColor: "#ffffff" }}>
          <LinearSequenceViewer
            sequence={sequence}
            annotations={allAnnotations as SeqVizAnnotation[]}
            height={viewerOuterHeight}
            centerOnPosition={centerOnPosition}
            topOnPosition={topOnPosition}
            selection={selection ? { start: selection.start, end: selection.end } : null}
            onSelectionChange={(sel) => handleSelectionNoRecenter(sel ? { start: sel.start, end: sel.end, clockwise: true } : null)}
            onAnnotationClick={(ann, idx) => openEditAnnotationModal(idx)}
          />
        </div>
        )}
      </div>

      {/* Bottom row: Annotations + Interactions panels */}
      <div className="flex gap-4">
        {/* Interactions Panel — functional relationships inferred by rule-based
            detector: promoter→CDS→polyA, RBS→CDS, Pol3→guide-RNA, operator
            repression, insulator boundaries, recombination-flanked cassettes */}
        {interactions.length > 0 && (
          <div
            className="rounded-2xl p-3 overflow-hidden flex flex-col"
            style={{ backgroundColor: "#1e3a34", flex: 1, minWidth: 0 }}
          >
            <div className="flex items-center justify-between mb-2">
              <div className="font-medium" style={{ color: "#dbefe7" }}>
                Interactions ({interactions.length})
              </div>
              <label className="text-xs flex items-center gap-1" style={{ color: "#dbefe7" }}>
                <input
                  type="checkbox"
                  checked={showInteractionChords}
                  onChange={(e) => setShowInteractionChords(e.target.checked)}
                />
                Show chords
              </label>
            </div>

            {interactionDescLoading && (
              <div className="text-xs mb-2" style={{ color: "#dbefe7", opacity: 0.7 }}>
                Generating description...
              </div>
            )}

            {interactionDescription && (
              <div
                className="text-xs mb-3 p-2 rounded whitespace-pre-wrap"
                style={{ backgroundColor: "rgba(219, 239, 231, 0.08)", color: "#dbefe7", maxHeight: "220px", overflowY: "auto" }}
              >
                {interactionDescription}
              </div>
            )}

            <div className="flex-1 overflow-y-auto space-y-1">
              {interactions.map((ix) => {
                const INTERACTION_COLORS: Record<string, string> = {
                  "http://identifiers.org/biomodels.sbo/SBO:0000589": "#2563eb",
                  "http://identifiers.org/biomodels.sbo/SBO:0000183": "#0ea5e9",
                  "http://identifiers.org/biomodels.sbo/SBO:0000184": "#16a34a",
                  "http://identifiers.org/biomodels.sbo/SBO:0000170": "#10b981",
                  "http://identifiers.org/biomodels.sbo/SBO:0000169": "#dc2626",
                  "http://identifiers.org/biomodels.sbo/SBO:0000182": "#f59e0b",
                  "http://identifiers.org/biomodels.sbo/SBO:0000178": "#f43f5e",
                };
                const color = INTERACTION_COLORS[ix.sbo_term || ""] || "#7c3aed";
                const isHL = highlightedInteractionId === ix.interaction_id;
                const names = (ix.participants || [])
                  .map((p) => p.name)
                  .filter(Boolean)
                  .join(" → ");
                return (
                  <div
                    key={ix.interaction_id}
                    onMouseEnter={() => setHighlightedInteractionId(ix.interaction_id)}
                    onMouseLeave={() => setHighlightedInteractionId(null)}
                    className="p-2 rounded cursor-pointer text-xs"
                    style={{
                      backgroundColor: isHL ? "rgba(219, 239, 231, 0.18)" : "rgba(219, 239, 231, 0.06)",
                      color: "#dbefe7",
                      borderLeft: `3px solid ${color}`,
                    }}
                  >
                    <div className="font-medium">{ix.interaction_type}{ix.rule_id ? ` (${ix.rule_id})` : ""}</div>
                    <div style={{ opacity: 0.85 }}>{names}</div>
                    {ix.notes && <div className="mt-1" style={{ opacity: 0.7, fontSize: "0.7rem" }}>{ix.notes}</div>}
                  </div>
                );
              })}
            </div>
          </div>
        )}

        {/* Annotations Side Panel */}
        <div
          className="rounded-2xl p-3 overflow-hidden flex flex-col"
          style={{ backgroundColor: "#2d4a3e", flex: 1, minWidth: 0 }}
        >
          <div className="font-medium mb-2" style={{ color: "#dbefe7" }}>
            Annotations ({annotations.length})
          </div>

          {/* Disclaimer */}
          <div className="text-xs mb-3 p-2 rounded" style={{ backgroundColor: "rgba(219, 239, 231, 0.1)", color: "#dbefe7", opacity: 0.8 }}>
            Note: Double-clicking split annotations (spanning lines or origin) may select the wrong annotation. Use this panel for reliable editing.
          </div>

          {/* Annotations list */}
          <div className="flex-1 overflow-y-auto space-y-1">
            {annotations.map((ann, idx) => {
              const isOriginCrossing = ann.start > ann.end;
              const coords = isOriginCrossing
                ? `${ann.start + 1}..${sequence.length}, 1..${ann.end}`
                : `${ann.start + 1}..${ann.end}`;
              const direction = ann.direction === -1 ? "←" : ann.direction === 1 ? "→" : "○";

              // Determine layer badge and color
              const layer = (ann as any).layer || "feature";
              const moduleType = (ann as any).module_type || "";
              const isModule = layer === "module";
              const isMotif = layer === "motif" || ann.type === "motif";
              const isGap = layer === "gap";
              
              // Color coding by layer type
              const getBgColor = () => {
                if (isModule) return "#4a7c59";  // Darker green for modules
                if (isMotif) return "#6b5b95";   // Purple for motifs
                if (isGap) return "#5a6c7d";     // Gray-blue for gaps
                return "#46896c";                // Default green for features
              };
              
              // Layer badge
              const getLayerBadge = () => {
                if (isModule) return "M";  // Module
                if (isMotif) return "T";   // moTif
                if (isGap) return "G";     // Gap
                return "";                 // Feature (no badge)
              };

              return (
                <button
                  key={idx}
                  onClick={() => openEditAnnotationModal(idx)}
                  className="w-full text-left rounded-lg p-2 transition hover:opacity-90"
                  style={{ backgroundColor: getBgColor() }}
                >
                  <div className="flex items-center gap-1.5">
                    {getLayerBadge() && (
                      <span 
                        className="text-xs px-1.5 py-0.5 rounded font-medium"
                        style={{ 
                          backgroundColor: "rgba(0,0,0,0.2)", 
                          color: "#dbefe7",
                          fontSize: "0.65rem"
                        }}
                      >
                        {getLayerBadge()}
                      </span>
                    )}
                    <div className="text-sm font-medium truncate flex-1" style={{ color: "#dbefe7" }}>
                      {direction} {ann.name}
                    </div>
                  </div>
                  <div className="text-xs flex items-center gap-2" style={{ color: "#dbefe7", opacity: 0.7 }}>
                    <span>{coords}</span>
                    {moduleType && (
                      <span className="italic">{moduleType.replace(/_/g, " ")}</span>
                    )}
                  </div>
                </button>
              );
            })}
            {annotations.length === 0 && (
              <div className="text-sm text-center py-4" style={{ color: "#dbefe7", opacity: 0.6 }}>
                No annotations yet
              </div>
            )}
          </div>
        </div>
      </div>
      </div>

      {/* Analyze Intent error */}
      {analyzeIntentError && (
        <div className="mt-2 text-sm px-3 py-1 rounded" style={{ backgroundColor: "rgba(239,68,68,0.2)", color: "#fca5a5" }}>
          Plasmid analysis error: {analyzeIntentError}
        </div>
      )}

      {/* Plasmid Analysis Panel */}
      {intentAnalysis && (
        <div className="mt-3 rounded-2xl p-4" style={{ backgroundColor: "#2d4a3e" }}>
          <div className="flex items-center justify-between mb-3">
            <span className="text-base font-semibold" style={{ color: "#dbefe7" }}>Plasmid Analysis</span>
            <button
              onClick={() => setIntentAnalysis(null)}
              className="text-xs px-2 py-0.5 rounded hover:opacity-80"
              style={{ backgroundColor: "#46896c", color: "#dbefe7" }}
            >
              Dismiss
            </button>
          </div>
          <div className="text-sm" style={{ color: "#dbefe7", lineHeight: 1.7 }}>
            {intentAnalysis.split('\n').map((line, i) => {
              // Handle headers
              if (line.startsWith('# ')) {
                return <h2 key={i} className="text-lg font-bold mt-4 mb-2" style={{ color: "#b8e0d2" }}>{line.slice(2)}</h2>;
              }
              if (line.startsWith('## ')) {
                return <h3 key={i} className="text-base font-semibold mt-3 mb-1" style={{ color: "#9dd5c2" }}>{line.slice(3)}</h3>;
              }
              if (line.startsWith('### ')) {
                return <h4 key={i} className="text-sm font-medium mt-2 mb-1" style={{ color: "#8ccbb5" }}>{line.slice(4)}</h4>;
              }
              // Handle bullet points
              if (line.startsWith('- ')) {
                const content = line.slice(2)
                  .replace(/\*\*([^*]+)\*\*/g, '<strong>$1</strong>')
                  .replace(/\*([^*]+)\*/g, '<em>$1</em>');
                return <div key={i} className="ml-4 my-0.5" dangerouslySetInnerHTML={{ __html: `• ${content}` }} />;
              }
              // Handle bold and italic in regular lines
              if (line.trim()) {
                const content = line
                  .replace(/\*\*([^*]+)\*\*/g, '<strong>$1</strong>')
                  .replace(/\*([^*]+)\*/g, '<em>$1</em>')
                  .replace(/→/g, '<span style="color:#9dd5c2">→</span>');
                return <p key={i} className="my-1" dangerouslySetInnerHTML={{ __html: content }} />;
              }
              return <div key={i} className="h-2" />;
            })}
          </div>
        </div>
      )}

      {/* Translation direction-prompt modal */}
      {translationModalOpen && (
        <div
          className="fixed inset-0 z-50 flex items-center justify-center"
          style={{ backgroundColor: "rgba(0, 0, 0, 0.5)" }}
          onClick={() => setTranslationModalOpen(false)}
        >
          <div
            className="rounded-2xl p-6 w-full max-w-md"
            style={{ backgroundColor: "#2d4a3e" }}
            onClick={(e) => e.stopPropagation()}
          >
            <h2 className="text-lg font-medium mb-3" style={{ color: "#dbefe7" }}>
              Add Translation Annotation
            </h2>
            <p className="text-sm mb-4" style={{ color: "#9dc5b7" }}>
              Choose the reading-frame direction for the selected region
              {hasRangeSelection
                ? ` (${(selectionStart ?? 0) + 1}–${selectionEnd ?? 0}, ${(selectionEnd ?? 0) - (selectionStart ?? 0)} bp)`
                : ""}.
            </p>
            <div className="space-y-2 mb-4">
              <label className="flex items-center gap-2 text-sm" style={{ color: "#dbefe7" }}>
                <input
                  type="radio"
                  name="translation-direction"
                  value="forward"
                  checked={translationDirection === "forward"}
                  onChange={() => setTranslationDirection("forward")}
                />
                Forward (5′ → 3′ on the top strand)
              </label>
              <label className="flex items-center gap-2 text-sm" style={{ color: "#dbefe7" }}>
                <input
                  type="radio"
                  name="translation-direction"
                  value="reverse"
                  checked={translationDirection === "reverse"}
                  onChange={() => setTranslationDirection("reverse")}
                />
                Reverse (3′ → 5′ — translate the reverse complement)
              </label>
            </div>
            <div className="flex justify-end gap-2">
              <button
                className="splicify-tb-btn is-outline"
                onClick={() => setTranslationModalOpen(false)}
              >
                Cancel
              </button>
              <button
                className="splicify-tb-btn"
                onClick={handleCreateTranslation}
              >
                Create
              </button>
            </div>
          </div>
        </div>
      )}

      {/* Add/Edit Annotation Modal */}
      {modalOpen && (
        <div
          className="fixed inset-0 z-50 flex items-center justify-center"
          style={{ backgroundColor: "rgba(0, 0, 0, 0.5)" }}
          onClick={handleCloseModal}
        >
          <div
            className="rounded-2xl p-6 w-full max-w-lg"
            style={{ backgroundColor: "#2d4a3e" }}
            onClick={(e) => e.stopPropagation()}
          >
            {/* Header */}
            <h2 className="text-lg font-medium mb-4" style={{ color: "#dbefe7" }}>
              {modalMode === "add" ? "Add Annotation" : formData.name}
            </h2>

            {/* Tab Bar (edit mode only) */}
            {modalMode === "edit" && (
              <div className="flex mb-4 border-b" style={{ borderColor: "#46896c" }}>
                <button
                  onClick={() => setActiveTab("info")}
                  className="px-4 py-2 text-sm font-medium transition-colors"
                  style={{
                    color: activeTab === "info" ? "#dbefe7" : "#9dc5b7",
                    borderBottom: activeTab === "info" ? "2px solid #dbefe7" : "2px solid transparent",
                    marginBottom: "-1px"
                  }}
                >
                  Info
                </button>
                <button
                  onClick={() => setActiveTab("edit")}
                  className="px-4 py-2 text-sm font-medium transition-colors"
                  style={{
                    color: activeTab === "edit" ? "#dbefe7" : "#9dc5b7",
                    borderBottom: activeTab === "edit" ? "2px solid #dbefe7" : "2px solid transparent",
                    marginBottom: "-1px"
                  }}
                >
                  Edit
                </button>
              </div>
            )}

            {/* Tab Content */}
            {(modalMode === "add" || activeTab === "edit") ? (
              /* Edit Form */
              <>
                <div className="space-y-4">
                  <div>
                    <label className="block text-sm mb-1" style={{ color: "#dbefe7" }}>Name *</label>
                    <input
                      type="text"
                      value={formData.name}
                      onChange={(e) => setFormData({ ...formData, name: e.target.value })}
                      className="w-full rounded-lg px-3 py-2 text-sm focus:outline-none focus:ring-2"
                      style={{
                        backgroundColor: "#46896c",
                        color: "#dbefe7",
                        border: formErrors.name ? "2px solid #ff6b6b" : "none",
                      }}
                      placeholder="e.g., GFP, Promoter, Overlap_1"
                    />
                    {formErrors.name && <div className="text-xs mt-1" style={{ color: "#ff6b6b" }}>{formErrors.name}</div>}
                  </div>

                  <div>
                    <label className="block text-sm mb-1" style={{ color: "#dbefe7" }}>Start</label>
                    <input
                      type="number"
                      value={formData.start}
                      onChange={(e) => setFormData({ ...formData, start: e.target.value })}
                      className="w-full rounded-lg px-3 py-2 text-sm focus:outline-none focus:ring-2"
                      style={{
                        backgroundColor: "#46896c",
                        color: "#dbefe7",
                        border: formErrors.start ? "2px solid #ff6b6b" : "none",
                      }}
                      placeholder="1"
                      min="1"
                      max={sequence.length}
                    />
                    {formErrors.start && <div className="text-xs mt-1" style={{ color: "#ff6b6b" }}>{formErrors.start}</div>}
                  </div>

                  <div>
                    <label className="block text-sm mb-1" style={{ color: "#dbefe7" }}>End</label>
                    <input
                      type="number"
                      value={formData.end}
                      onChange={(e) => setFormData({ ...formData, end: e.target.value })}
                      className="w-full rounded-lg px-3 py-2 text-sm focus:outline-none focus:ring-2"
                      style={{
                        backgroundColor: "#46896c",
                        color: "#dbefe7",
                        border: formErrors.end ? "2px solid #ff6b6b" : "none",
                      }}
                      placeholder={String(sequence.length)}
                      min="1"
                      max={sequence.length}
                    />
                    {formErrors.end && <div className="text-xs mt-1" style={{ color: "#ff6b6b" }}>{formErrors.end}</div>}
                    {parseInt(formData.start, 10) > parseInt(formData.end, 10) && formData.start && formData.end && (
                      <div className="text-xs mt-1" style={{ color: "#dbefe7", opacity: 0.7 }}>
                        Origin-crossing annotation
                      </div>
                    )}
                  </div>

                  <div>
                    <label className="block text-sm mb-1" style={{ color: "#dbefe7" }}>Direction</label>
                    <select
                      value={formData.direction}
                      onChange={(e) => setFormData({ ...formData, direction: e.target.value as "forward" | "reverse" | "none" })}
                      className="w-full rounded-lg px-3 py-2 text-sm focus:outline-none focus:ring-2"
                      style={{ backgroundColor: "#46896c", color: "#dbefe7" }}
                    >
                      <option value="forward">Forward (+)</option>
                      <option value="reverse">Reverse (-)</option>
                      <option value="none">None (no arrow)</option>
                    </select>
                  </div>
                </div>

                <div className="mt-6 flex items-center justify-between">
                  <div>
                    {modalMode === "edit" && editingIndex !== null && (
                      <button
                        onClick={() => handleDeleteAnnotation(editingIndex)}
                        className="rounded-lg px-4 py-2 text-sm font-medium hover:opacity-90"
                        style={{ backgroundColor: "#ff6b6b", color: "#fff" }}
                      >
                        Delete
                      </button>
                    )}
                  </div>
                  <div className="flex items-center gap-2">
                    <button
                      onClick={handleCloseModal}
                      className="rounded-lg px-4 py-2 text-sm font-medium hover:opacity-90"
                      style={{ backgroundColor: "#46896c", color: "#dbefe7" }}
                    >
                      Cancel
                    </button>
                    <button
                      onClick={handleSaveAnnotation}
                      className="rounded-lg px-4 py-2 text-sm font-medium hover:opacity-90"
                      style={{ backgroundColor: "#dbefe7", color: "#105b39" }}
                    >
                      Save
                    </button>
                  </div>
                </div>
              </>
            ) : (
              /* Info Tab - Gene Card Content */
              (() => {
                const annotation = editingIndex !== null ? annotations[editingIndex] : null;
                if (!annotation) return null;

                const strandLabel = annotation.direction === 1 ? "+" : annotation.direction === -1 ? "-" : "none";
                const kb = annotation.kb_data;

                const InfoRow = ({ label, value, link }: { label: string; value: string; link?: string }) => (
                  <div className="flex items-start mb-1">
                    <span className="text-xs font-medium" style={{ color: "#dbefe7", opacity: 0.7, minWidth: "80px" }}>
                      {label}:
                    </span>
                    {link ? (
                      <a href={link} target="_blank" rel="noopener noreferrer"
                         className="ml-2 text-sm underline hover:opacity-80" style={{ color: "#9dd5c2" }}>
                        {value}
                      </a>
                    ) : (
                      <span className="ml-2 text-sm" style={{ color: "#dbefe7" }}>{value}</span>
                    )}
                  </div>
                );

                return (
                  <div className="space-y-4">
                    {/* Feature Info Section */}
                    <div className="p-3 rounded-lg" style={{ backgroundColor: "#46896c" }}>
                      <h3 className="text-sm font-semibold mb-2" style={{ color: "#dbefe7" }}>
                        Feature
                      </h3>
                      <InfoRow label="Name" value={annotation.name} />
                      <InfoRow label="Type" value={annotation.type || "misc_feature"} />
                      <InfoRow label="Location" value={`${annotation.start + 1}..${annotation.end} (${strandLabel})`} />
                      <InfoRow label="Length" value={`${annotation.end - annotation.start} bp`} />
                      {annotation.description && <InfoRow label="Description" value={annotation.description} />}
                      {/* sgRNA-flanking primer-design shortcut */}
                      {(() => {
                        const t = (annotation.type || "").toLowerCase();
                        const isSgrna = t.includes("sgrna") || t.includes("guide") || t === "grna";
                        if (!isSgrna) return null;
                        return (
                          <div className="mt-3 pt-3" style={{ borderTop: "1px solid rgba(219, 239, 231, 0.2)" }}>
                            <button
                              onClick={() => handleOpenPrimerModalForSgrna(annotation)}
                              className="rounded-lg px-3 py-2 text-xs font-medium hover:opacity-90"
                              style={{ backgroundColor: "#dbefe7", color: "#105b39" }}
                            >
                              Design Primers Around sgRNA (±250 bp, exclude ±75 bp)
                            </button>
                          </div>
                        );
                      })()}
                    </div>

                    {/* pegRNA design section (if this is a pegRNA from /design_pegrnas) */}
                    {annotation.metadata?.pegrna_design && (
                      <div className="p-3 rounded-lg" style={{ backgroundColor: "#46896c" }}>
                        <h3 className="text-sm font-semibold mb-2" style={{ color: "#dbefe7" }}>
                          Prime Editing pegRNA{" "}
                          <span style={{ opacity: 0.7, fontWeight: 400 }}>
                            (#{annotation.metadata.rank} · eff {annotation.metadata.predicted_efficiency})
                          </span>
                        </h3>
                        <InfoRow label="Edit" value={`${annotation.metadata.edit_type} ${annotation.metadata.edit_ref || ""}->${annotation.metadata.edit_alt || ""} @ ${annotation.metadata.edit_start_1based}..${annotation.metadata.edit_end_1based}`} />
                        <InfoRow label="dPAM" value={annotation.metadata.is_dpam ? "yes" : "no"} />
                        <InfoRow label="PE3b" value={annotation.metadata.is_pe3b ? "yes" : "no"} />
                        <InfoRow label="Cas9 score" value={String(annotation.metadata.cas9_score)} />
                        <div className="flex items-start gap-2 py-1">
                          <span className="text-xs font-medium" style={{ color: "#dbefe7", opacity: 0.7, minWidth: 80 }}>Spacer:</span>
                          <code className="text-xs flex-1" style={{ color: "#dbefe7", fontFamily: '"Hack", monospace', wordBreak: "break-all" }}>{annotation.metadata.spacer}</code>
                          <button onClick={() => handleCopyToClipboard(annotation.metadata!.spacer)} className="text-xs px-2 py-0.5 rounded hover:opacity-90" style={{ backgroundColor: "#dbefe7", color: "#105b39" }}>Copy</button>
                        </div>
                        <div className="flex items-start gap-2 py-1">
                          <span className="text-xs font-medium" style={{ color: "#dbefe7", opacity: 0.7, minWidth: 80 }}>RTT:</span>
                          <code className="text-xs flex-1" style={{ color: "#dbefe7", fontFamily: '"Hack", monospace', wordBreak: "break-all" }}>{annotation.metadata.rtt}</code>
                          <span className="text-xs" style={{ color: "#dbefe7", opacity: 0.7 }}>{annotation.metadata.rtt_length}nt</span>
                        </div>
                        <div className="flex items-start gap-2 py-1">
                          <span className="text-xs font-medium" style={{ color: "#dbefe7", opacity: 0.7, minWidth: 80 }}>PBS:</span>
                          <code className="text-xs flex-1" style={{ color: "#dbefe7", fontFamily: '"Hack", monospace', wordBreak: "break-all" }}>{annotation.metadata.pbs}</code>
                          <span className="text-xs" style={{ color: "#dbefe7", opacity: 0.7 }}>{annotation.metadata.pbs_length}nt</span>
                        </div>
                        <div className="flex items-start gap-2 py-1">
                          <span className="text-xs font-medium" style={{ color: "#dbefe7", opacity: 0.7, minWidth: 80 }}>Full pegRNA:</span>
                          <code className="text-xs flex-1" style={{ color: "#dbefe7", fontFamily: '"Hack", monospace', wordBreak: "break-all" }}>{annotation.metadata.full_pegrna}</code>
                          <button onClick={() => handleCopyToClipboard(annotation.metadata!.full_pegrna)} className="text-xs px-2 py-0.5 rounded hover:opacity-90" style={{ backgroundColor: "#dbefe7", color: "#105b39" }}>Copy</button>
                        </div>
                        <div className="text-xs mt-1" style={{ color: "#dbefe7", opacity: 0.6 }}>5{'2192'}3 order: spacer + scaffold ({annotation.metadata.scaffold ? annotation.metadata.scaffold.length : 0} nt) + RTT + PBS · total {annotation.metadata.full_pegrna_length} nt.</div>
                        {annotation.metadata.ngrna && (
                          <div className="mt-3 pt-3" style={{ borderTop: "1px solid rgba(219, 239, 231, 0.2)" }}>
                            <h4 className="text-xs font-semibold mb-1" style={{ color: "#dbefe7" }}>
                              ngRNA{annotation.metadata.is_pe3b ? " (PE3b)" : " (PE3)"}
                            </h4>
                            <div className="flex items-start gap-2 py-1">
                              <span className="text-xs font-medium" style={{ color: "#dbefe7", opacity: 0.7, minWidth: 80 }}>Spacer:</span>
                              <code className="text-xs flex-1" style={{ color: "#dbefe7", fontFamily: '"Hack", monospace', wordBreak: "break-all" }}>{annotation.metadata.ngrna.spacer}</code>
                              <button onClick={() => handleCopyToClipboard(annotation.metadata!.ngrna.spacer)} className="text-xs px-2 py-0.5 rounded hover:opacity-90" style={{ backgroundColor: "#dbefe7", color: "#105b39" }}>Copy</button>
                            </div>
                            <InfoRow label="nick-to-peg" value={`${annotation.metadata.ngrna.nick_to_pegRNA} bp`} />
                            <InfoRow label="ngRNA Cas9" value={String(annotation.metadata.ngrna.cas9_score)} />
                          </div>
                        )}
                      </div>
                    )}

                    {/* Guide design section (if this is an sgRNA from /design_guides) */}
                    {annotation.metadata?.guide_design && (
                      <div className="p-3 rounded-lg" style={{ backgroundColor: "#46896c" }}>
                        <h3 className="text-sm font-semibold mb-2" style={{ color: "#dbefe7" }}>
                          CRISPR Guide
                        </h3>
                        {annotation.metadata.spacer && (
                          <div className="flex items-start gap-2 py-1">
                            <span className="text-xs font-medium" style={{ color: "#dbefe7", opacity: 0.7, minWidth: 80 }}>
                              Spacer:
                            </span>
                            <code className="text-xs flex-1" style={{ color: "#dbefe7", fontFamily: '"Hack", monospace', wordBreak: "break-all" }}>
                              {annotation.metadata.spacer}
                            </code>
                            <button
                              onClick={() => handleCopyToClipboard(annotation.metadata!.spacer)}
                              className="text-xs px-2 py-0.5 rounded hover:opacity-90"
                              style={{ backgroundColor: "#dbefe7", color: "#105b39" }}
                            >
                              Copy
                            </button>
                          </div>
                        )}
                        {annotation.metadata.pam && (
                          <InfoRow label="PAM" value={`${annotation.metadata.pam} (${annotation.metadata.pam_setting || "?"} · ${annotation.metadata.pam_position || "?"})`} />
                        )}
                        {typeof annotation.metadata.score === "number" && (
                          <InfoRow
                            label="Score"
                            value={`${annotation.metadata.score.toFixed(1)} / 100${
                              annotation.metadata.score_method
                                ? ` · ${annotation.metadata.score_method}`
                                : ""
                            }`}
                          />
                        )}
                        {annotation.metadata.context_30mer && (
                          <div className="flex items-start gap-2 py-1">
                            <span className="text-xs font-medium" style={{ color: "#dbefe7", opacity: 0.7, minWidth: 80 }}>
                              30-mer:
                            </span>
                            <code className="text-xs flex-1" style={{ color: "#dbefe7", fontFamily: '"Hack", monospace', wordBreak: "break-all" }}>
                              {annotation.metadata.context_30mer}
                            </code>
                          </div>
                        )}
                        {typeof annotation.metadata.gc_fraction === "number" && (
                          <InfoRow label="GC content" value={`${(annotation.metadata.gc_fraction * 100).toFixed(0)} %`} />
                        )}
                        {typeof annotation.metadata.max_homopolymer === "number" && (
                          <InfoRow label="Max homopolymer" value={`${annotation.metadata.max_homopolymer} nt`} />
                        )}
                        {typeof annotation.metadata.n_offtargets === "number" && (
                          <InfoRow
                            label="Off-targets"
                            value={
                              annotation.metadata.n_offtargets === 0
                                ? "0 (unique on this plasmid)"
                                : `${annotation.metadata.n_offtargets} additional hit(s) on this plasmid`
                            }
                          />
                        )}
                        {annotation.metadata.score_components && (
                          <div className="mt-2 pt-2" style={{ borderTop: "1px solid rgba(219, 239, 231, 0.2)" }}>
                            <div className="text-xs font-medium mb-1" style={{ color: "#dbefe7", opacity: 0.7 }}>
                              Score breakdown
                            </div>
                            <div className="text-xs grid grid-cols-2 gap-x-3 gap-y-0.5" style={{ color: "#dbefe7" }}>
                              {Object.entries(annotation.metadata.score_components).map(([k, v]) => (
                                <div key={k} className="flex justify-between">
                                  <span style={{ opacity: 0.75 }}>{k.replace(/_/g, " ")}</span>
                                  <span>{(v as number).toFixed(1)}</span>
                                </div>
                              ))}
                            </div>
                          </div>
                        )}
                        {annotation.metadata.design_region_1based && (
                          <div className="text-xs mt-2" style={{ color: "#dbefe7", opacity: 0.6 }}>
                            Designed against region {annotation.metadata.design_region_1based}
                          </div>
                        )}
                      </div>
                    )}

                    {/* Primer3 Design Section (if this is a primer designed by /design-primers) */}
                    {annotation.metadata?.primer_design && (
                      <div className="p-3 rounded-lg" style={{ backgroundColor: "#46896c" }}>
                        <h3 className="text-sm font-semibold mb-2" style={{ color: "#dbefe7" }}>
                          Primer3 Design
                        </h3>
                        <InfoRow label="Role" value={annotation.metadata.primer_role || ""} />
                        {annotation.metadata.primer_sequence && (
                          <div className="flex items-start gap-2 py-1">
                            <span className="text-xs font-medium" style={{ color: "#dbefe7", opacity: 0.7, minWidth: 80 }}>
                              Sequence:
                            </span>
                            <code className="text-xs flex-1" style={{ color: "#dbefe7", fontFamily: '"Hack", monospace', wordBreak: "break-all" }}>
                              {annotation.metadata.primer_sequence}
                            </code>
                            <button
                              onClick={() => handleCopyToClipboard(annotation.metadata!.primer_sequence)}
                              className="text-xs px-2 py-0.5 rounded hover:opacity-90"
                              style={{ backgroundColor: "#dbefe7", color: "#105b39" }}
                            >
                              Copy
                            </button>
                          </div>
                        )}
                        {typeof annotation.metadata.primer_tm === "number" && (
                          <InfoRow label="Tm" value={`${annotation.metadata.primer_tm.toFixed(1)} °C`} />
                        )}
                        {typeof annotation.metadata.product_size === "number" && (
                          <InfoRow label="Product size" value={`${annotation.metadata.product_size} bp`} />
                        )}
                        {typeof annotation.metadata.pair_penalty === "number" && (
                          <InfoRow label="Pair penalty" value={annotation.metadata.pair_penalty.toFixed(3)} />
                        )}
                        {annotation.metadata.amplicon_start_1based && annotation.metadata.amplicon_end_1based && (
                          <InfoRow
                            label="Amplicon"
                            value={`${annotation.metadata.amplicon_start_1based}..${annotation.metadata.amplicon_end_1based} (${annotation.metadata.amplicon_length} bp)`}
                          />
                        )}
                        {annotation.metadata.thermo_scores && (
                          <div className="text-xs mt-1" style={{ color: "#dbefe7", opacity: 0.85 }}>
                            Thermo: hairpin Tm{" "}
                            {annotation.metadata.thermo_scores.hairpin_th != null
                              ? `${annotation.metadata.thermo_scores.hairpin_th.toFixed(1)}°C`
                              : "—"}
                            {", self-dimer Tm "}
                            {annotation.metadata.thermo_scores.any_th != null
                              ? `${annotation.metadata.thermo_scores.any_th.toFixed(1)}°C`
                              : "—"}
                          </div>
                        )}
                        {Array.isArray(annotation.metadata.mispriming_sites) &&
                          annotation.metadata.mispriming_sites.length > 1 && (
                            <div className="text-xs mt-1" style={{ color: "#fbbf24" }}>
                              Mispriming: {annotation.metadata.mispriming_sites.length} exact hits on the template
                            </div>
                          )}
                        {annotation.metadata.amplicon && (
                          <div className="mt-2 pt-2" style={{ borderTop: "1px solid rgba(219, 239, 231, 0.2)" }}>
                            <div className="flex items-center justify-between mb-1">
                              <span className="text-xs font-medium" style={{ color: "#dbefe7", opacity: 0.7 }}>
                                Amplicon sequence:
                              </span>
                              <button
                                onClick={() => handleCopyToClipboard(annotation.metadata!.amplicon)}
                                className="text-xs px-2 py-0.5 rounded hover:opacity-90"
                                style={{ backgroundColor: "#dbefe7", color: "#105b39" }}
                              >
                                Copy amplicon
                              </button>
                            </div>
                            <textarea
                              readOnly
                              value={annotation.metadata.amplicon}
                              className="w-full text-xs rounded p-2"
                              style={{
                                backgroundColor: "#2d4a3e",
                                color: "#dbefe7",
                                fontFamily: '"Hack", monospace',
                                height: 80,
                                resize: "vertical",
                              }}
                            />
                          </div>
                        )}
                      </div>
                    )}

                    {/* Module Heuristics Section (if module with metadata) */}
                    {annotation.layer === "module" && annotation.metadata && (
                      <div className="p-3 rounded-lg" style={{ backgroundColor: "#46896c" }}>
                        <h3 className="text-sm font-semibold mb-2" style={{ color: "#dbefe7" }}>
                          Module Heuristics
                        </h3>
                        
                        {/* Module Classification */}
                        {annotation.module_type && <InfoRow label="Module Type" value={annotation.module_type.replace(/_/g, " ")} />}
                        {annotation.module_family && <InfoRow label="Module Family" value={annotation.module_family.replace(/_/g, " ")} />}
                        {annotation.payload_id && <InfoRow label="Payload ID" value={annotation.payload_id} />}
                        
                        {/* Component IDs - Show what features were used to call this module */}
                        {annotation.metadata.promoter_id && <InfoRow label="Promoter" value={annotation.metadata.promoter_id} />}
                        {annotation.metadata.polya_id && <InfoRow label="PolyA Signal" value={annotation.metadata.polya_id} />}
                        {annotation.metadata.marker_id && <InfoRow label="Marker" value={annotation.metadata.marker_id} />}
                        {annotation.metadata.payload_family && <InfoRow label="Payload Family" value={annotation.metadata.payload_family} />}
                        
                        {/* Pol3-specific metadata */}
                        {annotation.metadata.promoter_type && <InfoRow label="Promoter Type" value={annotation.metadata.promoter_type} />}
                        {annotation.metadata.terminator_type && <InfoRow label="Terminator Type" value={annotation.metadata.terminator_type.replace(/_/g, " ")} />}
                        {annotation.metadata.terminator_detected !== undefined && (
                          <InfoRow label="Terminator Detected" value={annotation.metadata.terminator_detected ? "Yes" : "No"} />
                        )}
                        {annotation.metadata.scaffold_present !== undefined && (
                          <InfoRow label="Scaffold Present" value={annotation.metadata.scaffold_present ? "Yes" : "No"} />
                        )}
                        
                        {/* Host and scope information */}
                        {annotation.metadata.host_scope && annotation.metadata.host_scope.length > 0 && (
                          <InfoRow label="Host Scope" value={annotation.metadata.host_scope.join(", ")} />
                        )}
                        
                        {/* Canonical IDs - Shows all feature IDs in this module */}
                        {annotation.metadata.canonical_ids && annotation.metadata.canonical_ids.length > 0 && (
                          <div className="mt-2 pt-2" style={{ borderTop: "1px solid rgba(219, 239, 231, 0.2)" }}>
                            <span className="text-xs font-medium" style={{ color: "#dbefe7", opacity: 0.7 }}>
                              Component Features:
                            </span>
                            <div className="mt-1 text-xs" style={{ color: "#dbefe7", opacity: 0.9 }}>
                              {annotation.metadata.canonical_ids.slice(0, 10).join(", ")}
                              {annotation.metadata.canonical_ids.length > 10 && ` (+${annotation.metadata.canonical_ids.length - 10} more)`}
                            </div>
                          </div>
                        )}
                        
                        {/* Flags */}
                        {annotation.metadata.has_viral_elements && (
                          <div className="mt-2 text-xs" style={{ color: "#fbbf24" }}>
                            Contains viral elements
                          </div>
                        )}
                        {annotation.metadata.has_guide_components && (
                          <div className="mt-1 text-xs" style={{ color: "#9dd5c2" }}>
                            Contains guide RNA components
                          </div>
                        )}
                        
                        {/* Detection notes */}
                        <div className="mt-2 pt-2 text-xs" style={{ borderTop: "1px solid rgba(219, 239, 231, 0.2)", color: "#dbefe7", opacity: 0.6 }}>
                          {annotation.source === "hierarchical_annotator" && "Detected by heuristic pipeline"}
                          {annotation.source === "grammar_annotator" && "Detected by grammar-based pipeline"}
                        </div>
                      </div>
                    )}


                    {/* KB Metadata Section (if available) */}
                    {kb && (
                      <div className="p-3 rounded-lg" style={{ backgroundColor: "#46896c" }}>
                        <h3 className="text-sm font-semibold mb-2" style={{ color: "#dbefe7" }}>
                          Knowledge Base {kb.source_type === "swissprot" ? "(SwissProt)" : kb.source_type === "feature_kb" ? "(Feature DB)" : ""}
                        </h3>
                        
                        {/* SwissProt fields */}
                        {kb.protein_name && <InfoRow label="Protein" value={kb.protein_name} />}
                        {kb.gene_name && <InfoRow label="Gene" value={kb.gene_name} />}
                        {kb.organism && <InfoRow label="Organism" value={kb.organism} />}
                        {kb.taxonomy_id && (
                          <InfoRow
                            label="Taxonomy"
                            value={kb.taxonomy_id}
                            link={`https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=${kb.taxonomy_id}`}
                          />
                        )}
                        {kb.protein_existence && <InfoRow label="Evidence" value={kb.protein_existence} />}
                        {kb.entry_name && (
                          <InfoRow
                            label="UniProt"
                            value={kb.entry_name}
                            link={`https://www.uniprot.org/uniprotkb/${kb.entry_name}`}
                          />
                        )}
                        
                        {/* Feature KB fields */}
                        {kb.feature_class && <InfoRow label="Class" value={kb.feature_class.replace(/_/g, " ")} />}
                        {kb.subclass && <InfoRow label="Subclass" value={kb.subclass.replace(/_/g, " ")} />}
                        {kb.polymerase_class && <InfoRow label="Polymerase" value={kb.polymerase_class.replace(/_/g, " ").toUpperCase()} />}
                        {kb.host_scope && kb.host_scope.length > 0 && kb.host_scope[0] !== "unknown" && (
                          <InfoRow label="Host Scope" value={kb.host_scope.join(", ")} />
                        )}
                        {kb.delivery_scope && kb.delivery_scope.length > 0 && (
                          <InfoRow label="Delivery" value={kb.delivery_scope.join(", ")} />
                        )}
                        {kb.annotation_source && <InfoRow label="Source" value={kb.annotation_source} />}
                        {kb.orientation_requirements && <InfoRow label="Orientation" value={kb.orientation_requirements} />}
                        
                        {/* Descriptions array */}
                        {kb.descriptions && kb.descriptions.length > 0 && (
                          <div className="mt-2 pt-2" style={{ borderTop: "1px solid rgba(219, 239, 231, 0.2)" }}>
                            <span className="text-xs font-medium" style={{ color: "#dbefe7", opacity: 0.7 }}>
                              Description:
                            </span>
                            <div className="mt-1 text-sm" style={{ color: "#dbefe7" }}>
                              {kb.descriptions.map((desc, i) => (
                                <p key={i} className="mb-1">{desc}</p>
                              ))}
                            </div>
                          </div>
                        )}
                      </div>
                    )}

                    {/* No KB data message */}
                    {!kb && !(annotation.layer === "module" && annotation.metadata) && (
                      <div className="p-3 rounded-lg text-sm" style={{ backgroundColor: "#46896c", color: "#dbefe7", opacity: 0.7 }}>
                        No additional knowledge base information available for this feature.
                      </div>
                    )}

                    {/* Close button */}
                    <div className="mt-4 flex justify-end">
                      <button
                        onClick={handleCloseModal}
                        className="rounded-lg px-4 py-2 text-sm font-medium hover:opacity-90"
                        style={{ backgroundColor: "#dbefe7", color: "#105b39" }}
                      >
                        Close
                      </button>
                    </div>
                  </div>
                );
              })()
            )}
          </div>
        </div>
      )}

      {/* Sequence Editor Modal */}
      {sequenceEditorOpen && (
        <div
          className="fixed inset-0 z-50 flex items-center justify-center"
          style={{ backgroundColor: "rgba(0, 0, 0, 0.5)" }}
          onClick={handleCloseSequenceEditor}
        >
          <div
            className="rounded-2xl p-6 w-full max-w-2xl"
            style={{ backgroundColor: "#2d4a3e" }}
            onClick={(e) => e.stopPropagation()}
          >
            <h2 className="text-lg font-medium mb-4" style={{ color: "#dbefe7" }}>Edit Sequence</h2>

            <div className="mb-2 text-sm" style={{ color: "#dbefe7", opacity: 0.8 }}>
              Current length: {sequence.length} bp | New length: {editedSequence.replace(/\s/g, "").length} bp
            </div>

            {editedSequence.replace(/\s/g, "").length !== sequence.length && (
              <div className="mb-2 text-sm" style={{ color: "#ffc107" }}>
                Warning: Changing sequence length may affect annotation coordinates.
              </div>
            )}

            <textarea
              value={editedSequence}
              onChange={(e) => {
                setEditedSequence(e.target.value);
                setSequenceError(validateDnaSequence(e.target.value.replace(/\s/g, "").toUpperCase()));
              }}
              className="w-full rounded-lg px-3 py-2 text-sm focus:outline-none focus:ring-2 resize-none"
              style={{
                backgroundColor: "#46896c",
                color: "#72777e",
                height: "300px",
                border: sequenceError ? "2px solid #ff6b6b" : "none",
                fontFamily: '"Hack", monospace',
              }}
              placeholder="Enter DNA sequence (A, T, G, C, N only)"
              spellCheck={false}
            />

            {sequenceError && <div className="text-sm mt-2" style={{ color: "#ff6b6b" }}>{sequenceError}</div>}

            <div className="mt-4 flex items-center justify-end gap-2">
              <button
                onClick={handleCloseSequenceEditor}
                className="rounded-lg px-4 py-2 text-sm font-medium hover:opacity-90"
                style={{ backgroundColor: "#46896c", color: "#dbefe7" }}
              >
                Cancel
              </button>
              <button
                onClick={handleSaveSequence}
                disabled={!!sequenceError}
                className="rounded-lg px-4 py-2 text-sm font-medium hover:opacity-90 disabled:opacity-50"
                style={{ backgroundColor: "#dbefe7", color: "#105b39" }}
              >
                Save
              </button>
            </div>
          </div>
        </div>
      )}

      {/* Add / Change Sequence Modal */}
      {addChangeModalOpen && (
        <div
          className="fixed inset-0 z-50 flex items-center justify-center"
          style={{ backgroundColor: "rgba(0, 0, 0, 0.5)" }}
          onClick={handleCloseAddChangeModal}
        >
          <div
            className="rounded-2xl p-6 w-full max-w-lg"
            style={{ backgroundColor: "#2d4a3e" }}
            onClick={(e) => e.stopPropagation()}
          >
            <h2 className="text-lg font-medium mb-4" style={{ color: "#dbefe7" }}>
              {hasRangeSelection ? "Modify Selected Sequence" : "Insert Sequence"}
            </h2>

            {hasRangeSelection ? (
              <div className="mb-3 text-sm" style={{ color: "#dbefe7", opacity: 0.8 }}>
                Selected region: {selectionStart! + 1}..{selectionEnd} ({selectionEnd! - selectionStart!} bp)
              </div>
            ) : (
              <div className="mb-3 text-sm" style={{ color: "#dbefe7", opacity: 0.8 }}>
                {selection ? `Insert at position ${selection.start + 1}` : `Insert at end of sequence (position ${sequence.length + 1})`}
              </div>
            )}

            <div className="mb-4">
              <label className="block text-sm mb-1" style={{ color: "#dbefe7" }}>
                {hasRangeSelection ? "Replace with sequence:" : "Sequence to insert:"}
              </label>
              <textarea
                value={newSequenceText}
                onChange={(e) => {
                  setNewSequenceText(e.target.value.toUpperCase());
                  const cleaned = e.target.value.replace(/\s/g, "").toUpperCase();
                  setAddChangeError(cleaned.length > 0 ? validateDnaSequence(cleaned) : "");
                }}
                className="w-full rounded-lg px-3 py-2 text-sm focus:outline-none focus:ring-2 resize-none"
                style={{
                  backgroundColor: "#46896c",
                  color: "#72777e",
                  height: "120px",
                  border: addChangeError ? "2px solid #ff6b6b" : "none",
                  fontFamily: '"Hack", monospace',
                }}
                placeholder="Enter DNA sequence (A, T, G, C, N only)"
                spellCheck={false}
              />
              {addChangeError && <div className="text-xs mt-1" style={{ color: "#ff6b6b" }}>{addChangeError}</div>}
              <div className="text-xs mt-1" style={{ color: "#dbefe7", opacity: 0.7 }}>
                Length: {newSequenceText.replace(/\s/g, "").length} bp
                {hasRangeSelection && ` (original: ${selectionEnd! - selectionStart!} bp)`}
              </div>
            </div>

            <div className="flex items-center justify-between">
              <div>
                {hasRangeSelection && (
                  <button
                    onClick={handleDeleteSelected}
                    className="rounded-lg px-4 py-2 text-sm font-medium hover:opacity-90"
                    style={{ backgroundColor: "#ff6b6b", color: "#fff" }}
                  >
                    Delete Selected
                  </button>
                )}
              </div>
              <div className="flex items-center gap-2">
                <button
                  onClick={handleCloseAddChangeModal}
                  className="rounded-lg px-4 py-2 text-sm font-medium hover:opacity-90"
                  style={{ backgroundColor: "#46896c", color: "#dbefe7" }}
                >
                  Cancel
                </button>
                <button
                  onClick={handleApplyAddChange}
                  disabled={!!addChangeError}
                  className="rounded-lg px-4 py-2 text-sm font-medium hover:opacity-90 disabled:opacity-50"
                  style={{ backgroundColor: "#dbefe7", color: "#105b39" }}
                >
                  {hasRangeSelection ? "Replace" : "Insert"}
                </button>
              </div>
            </div>
          </div>
        </div>
      )}

      {/* Import Annotations Modal */}
      {importModalOpen && (
        <div
          className="fixed inset-0 z-50 flex items-center justify-center"
          style={{ backgroundColor: "rgba(0, 0, 0, 0.5)" }}
          onClick={handleCloseImportModal}
        >
          <div
            className="rounded-2xl p-6 w-full max-w-xl"
            style={{ backgroundColor: "#2d4a3e", maxHeight: "85vh", overflowY: "auto" }}
            onClick={(e) => e.stopPropagation()}
          >
            <h2 className="text-lg font-medium mb-2" style={{ color: "#dbefe7" }}>
              Import Annotations from CSV
            </h2>
            <p className="text-sm mb-4" style={{ color: "#dbefe7", opacity: 0.85 }}>
              Upload a CSV of features and we'll map each <code>sequence</code> onto the
              loaded plasmid by exact identity (forward + reverse strand, with circular
              wrap-around). Required headers: <strong>name</strong>, <strong>sequence</strong>.
              Optional headers: <em>type, location, length, description</em> — these are
              displayed in the same gene-card popup as knowledge-base features.
            </p>

            <div className="mb-4">
              <button
                onClick={handleDownloadImportTemplate}
                className="rounded-lg px-3 py-2 text-xs font-medium hover:opacity-90"
                style={{ backgroundColor: "#46896c", color: "#dbefe7" }}
              >
                Download Template CSV
              </button>
            </div>

            <div className="mb-4">
              <input
                ref={importFileInputRef}
                type="file"
                accept=".csv,text/csv"
                onChange={handleImportFileSelected}
                style={{ display: "none" }}
              />
              <button
                onClick={handleClickImportFile}
                className="rounded-lg px-3 py-2 text-xs font-medium hover:opacity-90"
                style={{ backgroundColor: "#46896c", color: "#dbefe7" }}
              >
                Upload Annotations CSV
              </button>
              <span className="ml-3 text-xs" style={{ color: "#dbefe7", opacity: 0.85 }}>
                {importFileName || "No file selected."}
              </span>
            </div>

            <div className="mb-4">
              <label className="block text-sm mb-1" style={{ color: "#dbefe7" }}>
                Max mismatches per query (0 = exact only)
              </label>
              <input
                type="number"
                min={0}
                max={5}
                value={importMaxMismatches}
                onChange={(e) => setImportMaxMismatches(Math.max(0, Math.min(5, parseInt(e.target.value || "0", 10))))}
                className="rounded-lg px-3 py-1 text-sm"
                style={{ backgroundColor: "#46896c", color: "#dbefe7", width: "80px", border: "none" }}
              />
              <span className="ml-2 text-xs" style={{ color: "#dbefe7", opacity: 0.7 }}>
                Slower at &gt;0; useful for noisy guide / primer panels.
              </span>
            </div>

            {importError && (
              <div className="mb-3 text-sm rounded-lg px-3 py-2"
                   style={{ backgroundColor: "#ff6b6b", color: "#fff" }}>
                {importError}
              </div>
            )}

            {importResult && (
              <div className="mb-3 rounded-lg px-3 py-2 text-sm"
                   style={{ backgroundColor: "#46896c", color: "#dbefe7" }}>
                <div>
                  Matched <strong>{importResult.summary.n_matched}</strong> of{" "}
                  <strong>{importResult.summary.n_input}</strong> entries →{" "}
                  <strong>{importResult.summary.n_annotations}</strong> annotation
                  {importResult.summary.n_annotations === 1 ? "" : "s"} placed.
                </div>
                {importResult.unmatched.length > 0 && (
                  <details className="mt-2">
                    <summary className="cursor-pointer">
                      {importResult.unmatched.length} unmatched
                    </summary>
                    <ul className="mt-1 ml-4 list-disc text-xs" style={{ opacity: 0.85 }}>
                      {importResult.unmatched.map((u, i) => (
                        <li key={i}>
                          <code>{u.name || "(no name)"}</code> — {u.reason}
                        </li>
                      ))}
                    </ul>
                  </details>
                )}
              </div>
            )}

            <div className="flex items-center justify-end gap-2">
              <button
                onClick={handleCloseImportModal}
                className="rounded-lg px-4 py-2 text-sm font-medium hover:opacity-90"
                style={{ backgroundColor: "#46896c", color: "#dbefe7" }}
              >
                {importResult ? "Close" : "Cancel"}
              </button>
              <button
                onClick={handleSubmitImport}
                disabled={importLoading || !sequence}
                className="rounded-lg px-4 py-2 text-sm font-medium hover:opacity-90 disabled:opacity-50"
                style={{ backgroundColor: "#dbefe7", color: "#105b39" }}
              >
                {importLoading ? "Importing…" : "Import"}
              </button>
            </div>
          </div>
        </div>
      )}

      {/* Design Primers Modal */}
      {primerModalOpen && (
        <div
          className="fixed inset-0 z-50 flex items-center justify-center"
          style={{ backgroundColor: "rgba(0, 0, 0, 0.5)" }}
          onClick={handleClosePrimerModal}
        >
          <div
            className="rounded-2xl p-6 w-full max-w-2xl"
            style={{ backgroundColor: "#2d4a3e", maxHeight: "90vh", overflowY: "auto" }}
            onClick={(e) => e.stopPropagation()}
          >
            <h2 className="text-lg font-medium mb-2" style={{ color: "#dbefe7" }}>
              Design Primers (Primer3)
            </h2>
            <p className="text-sm mb-4" style={{ color: "#dbefe7", opacity: 0.85 }}>
              Designs a forward + reverse primer pair flanking the chosen region. Coordinates are
              1-indexed and inclusive. Pre-filled from the current selection if any.
            </p>

            {/* Application selector — drives adapter injection, default name
                suffix, and amplicon-size guards. */}
            <div className="mb-4">
              <label className="block text-xs mb-1" style={{ color: "#dbefe7", opacity: 0.85 }}>
                Application
              </label>
              <select
                value={primerForm.application}
                onChange={handlePrimerFormChange("application") as any}
                className="rounded-lg px-3 py-1 text-sm w-full"
                style={{ backgroundColor: "#46896c", color: "#dbefe7", border: "none" }}
              >
                <option value="fragment">Fragment (cloning / general PCR)</option>
                <option value="sanger">Sanger sequencing</option>
                <option value="illumina">Illumina (Nextera adapters, &le; 600 bp)</option>
              </select>
              <span className="text-xs" style={{ color: "#dbefe7", opacity: 0.7 }}>
                {primerForm.application === "illumina"
                  ? "Nextera adapters are prepended to each primer; Tm is calculated for the annealing portion only. Amplicons over 600 bp are rejected."
                  : primerForm.application === "sanger"
                  ? "Designed primers are scored for Sanger trace quality (Tm, dimers, read window, template structure, mispriming)."
                  : "Standard PCR primers, no adapter overhang."}
              </span>
            </div>

            <div className="grid grid-cols-2 gap-3 mb-4">
              <div>
                <label className="block text-xs mb-1" style={{ color: "#dbefe7", opacity: 0.85 }}>
                  Region start
                </label>
                <input
                  type="number" min={1}
                  value={primerForm.region_start}
                  onChange={handlePrimerFormChange("region_start")}
                  className="rounded-lg px-3 py-1 text-sm w-full"
                  style={{ backgroundColor: "#46896c", color: "#dbefe7", border: "none" }}
                />
              </div>
              <div>
                <label className="block text-xs mb-1" style={{ color: "#dbefe7", opacity: 0.85 }}>
                  Region end
                </label>
                <input
                  type="number" min={1}
                  value={primerForm.region_end}
                  onChange={handlePrimerFormChange("region_end")}
                  className="rounded-lg px-3 py-1 text-sm w-full"
                  style={{ backgroundColor: "#46896c", color: "#dbefe7", border: "none" }}
                />
              </div>
              <div>
                <label className="block text-xs mb-1" style={{ color: "#dbefe7", opacity: 0.85 }}>
                  Excluded start (optional)
                </label>
                <input
                  type="number" min={1}
                  value={primerForm.excluded_start}
                  onChange={handlePrimerFormChange("excluded_start")}
                  placeholder="e.g. 12100"
                  className="rounded-lg px-3 py-1 text-sm w-full"
                  style={{ backgroundColor: "#46896c", color: "#dbefe7", border: "none" }}
                />
              </div>
              <div>
                <label className="block text-xs mb-1" style={{ color: "#dbefe7", opacity: 0.85 }}>
                  Excluded end (optional)
                </label>
                <input
                  type="number" min={1}
                  value={primerForm.excluded_end}
                  onChange={handlePrimerFormChange("excluded_end")}
                  placeholder="e.g. 12180"
                  className="rounded-lg px-3 py-1 text-sm w-full"
                  style={{ backgroundColor: "#46896c", color: "#dbefe7", border: "none" }}
                />
              </div>
              <div>
                <label className="block text-xs mb-1" style={{ color: "#dbefe7", opacity: 0.85 }}>
                  Product size min (bp)
                </label>
                <input
                  type="number" min={50}
                  value={primerForm.product_size_min}
                  onChange={handlePrimerFormChange("product_size_min")}
                  className="rounded-lg px-3 py-1 text-sm w-full"
                  style={{ backgroundColor: "#46896c", color: "#dbefe7", border: "none" }}
                />
              </div>
              <div>
                <label className="block text-xs mb-1" style={{ color: "#dbefe7", opacity: 0.85 }}>
                  Product size max (bp)
                </label>
                <input
                  type="number" min={50}
                  value={primerForm.product_size_max}
                  onChange={handlePrimerFormChange("product_size_max")}
                  className="rounded-lg px-3 py-1 text-sm w-full"
                  style={{ backgroundColor: "#46896c", color: "#dbefe7", border: "none" }}
                />
              </div>
              <div>
                <label className="block text-xs mb-1" style={{ color: "#dbefe7", opacity: 0.85 }}>
                  Tm min / opt / max (°C, optional)
                </label>
                <div className="flex gap-1">
                  <input type="number" step="0.1" placeholder="min"
                    value={primerForm.primer_min_tm}
                    onChange={handlePrimerFormChange("primer_min_tm")}
                    className="rounded-lg px-2 py-1 text-sm w-1/3"
                    style={{ backgroundColor: "#46896c", color: "#dbefe7", border: "none" }}/>
                  <input type="number" step="0.1" placeholder="opt"
                    value={primerForm.primer_opt_tm}
                    onChange={handlePrimerFormChange("primer_opt_tm")}
                    className="rounded-lg px-2 py-1 text-sm w-1/3"
                    style={{ backgroundColor: "#46896c", color: "#dbefe7", border: "none" }}/>
                  <input type="number" step="0.1" placeholder="max"
                    value={primerForm.primer_max_tm}
                    onChange={handlePrimerFormChange("primer_max_tm")}
                    className="rounded-lg px-2 py-1 text-sm w-1/3"
                    style={{ backgroundColor: "#46896c", color: "#dbefe7", border: "none" }}/>
                </div>
              </div>
              <div>
                <label className="block text-xs mb-1" style={{ color: "#dbefe7", opacity: 0.85 }}>
                  Length min / opt / max (optional)
                </label>
                <div className="flex gap-1">
                  <input type="number" placeholder="min"
                    value={primerForm.primer_min_size}
                    onChange={handlePrimerFormChange("primer_min_size")}
                    className="rounded-lg px-2 py-1 text-sm w-1/3"
                    style={{ backgroundColor: "#46896c", color: "#dbefe7", border: "none" }}/>
                  <input type="number" placeholder="opt"
                    value={primerForm.primer_opt_size}
                    onChange={handlePrimerFormChange("primer_opt_size")}
                    className="rounded-lg px-2 py-1 text-sm w-1/3"
                    style={{ backgroundColor: "#46896c", color: "#dbefe7", border: "none" }}/>
                  <input type="number" placeholder="max"
                    value={primerForm.primer_max_size}
                    onChange={handlePrimerFormChange("primer_max_size")}
                    className="rounded-lg px-2 py-1 text-sm w-1/3"
                    style={{ backgroundColor: "#46896c", color: "#dbefe7", border: "none" }}/>
                </div>
              </div>
              <div>
                <label className="block text-xs mb-1" style={{ color: "#dbefe7", opacity: 0.85 }}>
                  Pair label (optional, used as prefix)
                </label>
                <input
                  type="text"
                  value={primerForm.pair_label}
                  onChange={handlePrimerFormChange("pair_label")}
                  placeholder="e.g. cirbp_exon4"
                  className="rounded-lg px-3 py-1 text-sm w-full"
                  style={{ backgroundColor: "#46896c", color: "#dbefe7", border: "none" }}
                />
              </div>
              <div>
                <label className="block text-xs mb-1" style={{ color: "#dbefe7", opacity: 0.85 }}>
                  Candidate pairs to evaluate
                </label>
                <input
                  type="number" min={1} max={20}
                  value={primerForm.num_return}
                  onChange={handlePrimerFormChange("num_return")}
                  className="rounded-lg px-3 py-1 text-sm w-full"
                  style={{ backgroundColor: "#46896c", color: "#dbefe7", border: "none" }}
                />
              </div>
              <div>
                <label className="block text-xs mb-1" style={{ color: "#dbefe7", opacity: 0.85 }}>
                  Forward primer name (overrides label_F)
                </label>
                <input
                  type="text"
                  value={primerForm.primer_fwd_name}
                  onChange={handlePrimerFormChange("primer_fwd_name")}
                  placeholder="auto"
                  className="rounded-lg px-3 py-1 text-sm w-full"
                  style={{ backgroundColor: "#46896c", color: "#dbefe7", border: "none" }}
                />
              </div>
              <div>
                <label className="block text-xs mb-1" style={{ color: "#dbefe7", opacity: 0.85 }}>
                  Reverse primer name (overrides label_R)
                </label>
                <input
                  type="text"
                  value={primerForm.primer_rev_name}
                  onChange={handlePrimerFormChange("primer_rev_name")}
                  placeholder="auto"
                  className="rounded-lg px-3 py-1 text-sm w-full"
                  style={{ backgroundColor: "#46896c", color: "#dbefe7", border: "none" }}
                />
              </div>
            </div>

            {primerError && (
              <div className="mb-3 text-sm rounded-lg px-3 py-2"
                   style={{ backgroundColor: "#ff6b6b", color: "#fff" }}>
                {primerError}
              </div>
            )}

            {primerResult && (
              <div className="mb-3 rounded-lg px-3 py-2 text-sm"
                   style={{ backgroundColor: "#46896c", color: "#dbefe7" }}>
                <div className="font-medium mb-1">
                  Designed pair · {primerResult.application || "fragment"} ·
                  product {primerResult.product_size} bp
                  {typeof primerResult.pair_penalty === "number" && (
                    <span style={{ opacity: 0.75 }}>
                      {" "}· penalty {primerResult.pair_penalty.toFixed(2)}
                    </span>
                  )}
                </div>
                {primerResult.selection_rationale && (
                  <div className="text-xs mb-1" style={{ opacity: 0.85 }}>
                    {primerResult.selection_method === "sanger_aware" ? "\u2728 " : ""}
                    {primerResult.selection_rationale}
                  </div>
                )}
                {/* Each primer: render adapter (if any) in muted text and
                    the annealing portion in bold so the user can see what
                    actually binds the template. */}
                {[
                  { role: "Forward", adapter: primerResult.left_adapter, anneal: primerResult.left_annealing, full: primerResult.left_primer, tm: primerResult.left_tm },
                  { role: "Reverse", adapter: primerResult.right_adapter, anneal: primerResult.right_annealing, full: primerResult.right_primer, tm: primerResult.right_tm },
                ].map((p) => (
                  <div key={p.role} className="text-xs" style={{ wordBreak: "break-all" }}>
                    {p.role}:{" "}
                    <code>
                      {p.adapter ? (
                        <span style={{ opacity: 0.7 }}>{p.adapter}</span>
                      ) : null}
                      <strong>{p.anneal}</strong>
                    </code>
                    {typeof p.tm === "number" && (
                      <span> · Tm {p.tm.toFixed(1)}°C (annealing)</span>
                    )}
                  </div>
                ))}
                <div className="text-xs mt-1" style={{ opacity: 0.85 }}>
                  Amplicon: {primerResult.amplicon_start_1based}..{primerResult.amplicon_end_1based}
                  {" "}({primerResult.amplicon.length} bp). Click either primer on the plasmid to
                  copy the amplicon sequence.
                </div>

                {/* Sanger-quality breakdown when scoring is available. */}
                {Array.isArray(primerResult.sanger_scores) && primerResult.sanger_scores.length > 0 && (
                  <div className="text-xs mt-2 pt-2" style={{ borderTop: "1px solid rgba(255,255,255,0.18)" }}>
                    <div style={{ fontWeight: 600, marginBottom: 2 }}>Sanger sequencing quality</div>
                    {primerResult.sanger_scores.map((s: any) => (
                      <div key={s.primer_index} style={{ marginBottom: 2 }}>
                        {s.direction === "forward" ? "Forward" : "Reverse"}: {Math.round(s.overall_score)}/100 · <em>{s.rating}</em>
                        {Array.isArray(s.warnings) && s.warnings.length > 0 && (
                          <ul style={{ marginTop: 2, marginLeft: 12, listStyle: "disc", opacity: 0.9 }}>
                            {s.warnings.slice(0, 4).map((w: string, i: number) => (
                              <li key={i}>{w}</li>
                            ))}
                          </ul>
                        )}
                      </div>
                    ))}
                  </div>
                )}
              </div>
            )}

            <div className="flex items-center justify-end gap-2">
              {primerResult && (
                <button
                  onClick={handleDownloadPrimerCsv}
                  className="rounded-lg px-4 py-2 text-sm font-medium hover:opacity-90"
                  style={{ backgroundColor: "#dbefe7", color: "#105b39" }}
                >
                  Download CSV
                </button>
              )}
              <button
                onClick={handleClosePrimerModal}
                className="rounded-lg px-4 py-2 text-sm font-medium hover:opacity-90"
                style={{ backgroundColor: "#46896c", color: "#dbefe7" }}
              >
                {primerResult ? "Close" : "Cancel"}
              </button>
              <button
                onClick={handleSubmitPrimerDesign}
                disabled={primerLoading || !sequence}
                className="rounded-lg px-4 py-2 text-sm font-medium hover:opacity-90 disabled:opacity-50"
                style={{ backgroundColor: "#dbefe7", color: "#105b39" }}
              >
                {primerLoading ? "Designing…" : "Design"}
              </button>
            </div>
          </div>
        </div>
      )}

      {/* Design Guides Modal */}
      {guideModalOpen && (
        <div
          className="fixed inset-0 z-50 flex items-center justify-center"
          style={{ backgroundColor: "rgba(0, 0, 0, 0.5)" }}
          onClick={handleCloseGuideModal}
        >
          <div
            className="rounded-2xl p-6 w-full max-w-2xl"
            style={{ backgroundColor: "#2d4a3e", maxHeight: "90vh", overflowY: "auto" }}
            onClick={(e) => e.stopPropagation()}
          >
            <h2 className="text-lg font-medium mb-2" style={{ color: "#dbefe7" }}>
              Design CRISPR Guides
            </h2>
            <p className="text-sm mb-4" style={{ color: "#dbefe7", opacity: 0.85 }}>
              Scans the chosen region on both strands for protospacers adjacent to the
              specified PAM, scores each candidate (0–100), and adds the top results as
              sgRNA annotations coloured by score (red → green). The default scorer is a
              pure-Python re-implementation of <strong>Doench 2014 Rule&nbsp;Set&nbsp;1</strong>
              (Nature Biotechnol. 32:1262), calibrated for SpCas9 + NGG. Off-target counts
              consider only the loaded plasmid, not a genome index.
            </p>

            <div className="grid grid-cols-2 gap-3 mb-4">
              <div>
                <label className="block text-xs mb-1" style={{ color: "#dbefe7", opacity: 0.85 }}>
                  Region start
                </label>
                <input type="number" min={1}
                  value={guideForm.region_start}
                  onChange={handleGuideFormChange("region_start")}
                  className="rounded-lg px-3 py-1 text-sm w-full"
                  style={{ backgroundColor: "#46896c", color: "#dbefe7", border: "none" }} />
              </div>
              <div>
                <label className="block text-xs mb-1" style={{ color: "#dbefe7", opacity: 0.85 }}>
                  Region end
                </label>
                <input type="number" min={1}
                  value={guideForm.region_end}
                  onChange={handleGuideFormChange("region_end")}
                  className="rounded-lg px-3 py-1 text-sm w-full"
                  style={{ backgroundColor: "#46896c", color: "#dbefe7", border: "none" }} />
              </div>
              <div>
                <label className="block text-xs mb-1" style={{ color: "#dbefe7", opacity: 0.85 }}>
                  PAM (IUPAC)
                </label>
                <input type="text"
                  value={guideForm.pam}
                  onChange={handleGuideFormChange("pam")}
                  placeholder="NGG"
                  className="rounded-lg px-3 py-1 text-sm w-full"
                  style={{ backgroundColor: "#46896c", color: "#dbefe7", border: "none" }} />
                <span className="text-xs" style={{ color: "#dbefe7", opacity: 0.6 }}>
                  Cas9 = NGG, Cas12a = TTTV, SaCas9 = NNGRRT
                </span>
              </div>
              <div>
                <label className="block text-xs mb-1" style={{ color: "#dbefe7", opacity: 0.85 }}>
                  PAM position
                </label>
                <select
                  value={guideForm.pam_position}
                  onChange={handleGuideFormChange("pam_position") as any}
                  className="rounded-lg px-3 py-1 text-sm w-full"
                  style={{ backgroundColor: "#46896c", color: "#dbefe7", border: "none" }}
                >
                  <option value="3prime">3' of protospacer (Cas9)</option>
                  <option value="5prime">5' of protospacer (Cas12a)</option>
                </select>
              </div>
              <div>
                <label className="block text-xs mb-1" style={{ color: "#dbefe7", opacity: 0.85 }}>
                  Guide length (nt)
                </label>
                <input type="number" min={16} max={25}
                  value={guideForm.guide_length}
                  onChange={handleGuideFormChange("guide_length")}
                  className="rounded-lg px-3 py-1 text-sm w-full"
                  style={{ backgroundColor: "#46896c", color: "#dbefe7", border: "none" }} />
              </div>
              <div>
                <label className="block text-xs mb-1" style={{ color: "#dbefe7", opacity: 0.85 }}>
                  Max guides to return
                </label>
                <input type="number" min={1} max={200}
                  value={guideForm.max_guides}
                  onChange={handleGuideFormChange("max_guides")}
                  className="rounded-lg px-3 py-1 text-sm w-full"
                  style={{ backgroundColor: "#46896c", color: "#dbefe7", border: "none" }} />
              </div>
              <div>
                <label className="block text-xs mb-1" style={{ color: "#dbefe7", opacity: 0.85 }}>
                  Minimum score (0–100, optional)
                </label>
                <input type="number" min={0} max={100}
                  value={guideForm.min_score}
                  onChange={handleGuideFormChange("min_score")}
                  className="rounded-lg px-3 py-1 text-sm w-full"
                  style={{ backgroundColor: "#46896c", color: "#dbefe7", border: "none" }} />
              </div>
              <div className="col-span-2">
                <label className="block text-xs mb-1" style={{ color: "#dbefe7", opacity: 0.85 }}>
                  Score method
                </label>
                <select
                  value={guideForm.score_method}
                  onChange={handleGuideFormChange("score_method") as any}
                  className="rounded-lg px-3 py-1 text-sm w-full"
                  style={{ backgroundColor: "#46896c", color: "#dbefe7", border: "none" }}
                >
                  <option value="doench2014">Doench 2014 Rule Set 1 (literature, SpCas9 + NGG)</option>
                  <option value="heuristic">Heuristic (GC + Pol III + homopolymer + seed + off-target)</option>
                </select>
                <span className="text-xs block mt-1" style={{ color: "#dbefe7", opacity: 0.6 }}>
                  Doench requires PAM=NGG, length=20, PAM 3'; otherwise the heuristic is used as a
                  fallback automatically. Base-editor scoring (BE-HIVE) requires a separate model bundle —
                  not yet integrated.
                </span>
              </div>
            </div>

            {guideError && (
              <div className="mb-3 text-sm rounded-lg px-3 py-2"
                   style={{ backgroundColor: "#ff6b6b", color: "#fff" }}>
                {guideError}
              </div>
            )}
            {guideResult && (
              <div className="mb-3 rounded-lg px-3 py-2 text-sm"
                   style={{ backgroundColor: "#46896c", color: "#dbefe7" }}>
                <div>
                  Scanned <strong>{guideResult.summary.n_candidates}</strong> protospacers ·{" "}
                  added <strong>{guideResult.n_added}</strong> sgRNA annotation
                  {guideResult.n_added === 1 ? "" : "s"} · scored with{" "}
                  <strong>{guideResult.summary.score_method}</strong>
                  {guideResult.summary.score_method !== guideResult.summary.score_method_requested && (
                    <span style={{ opacity: 0.7 }}>
                      {" "}(requested {guideResult.summary.score_method_requested}, fell back —
                      Doench requires PAM=NGG and length=20)
                    </span>
                  )}.
                </div>
                <div className="text-xs mt-1" style={{ opacity: 0.85 }}>
                  Click any guide on the plasmid to view its full score breakdown and copy the spacer.
                </div>
              </div>
            )}

            <div className="flex items-center justify-end gap-2">
              {guideResult && guideResult.guides && guideResult.guides.length > 0 && (
                <button
                  onClick={handleDownloadGuideCsv}
                  className="rounded-lg px-4 py-2 text-sm font-medium hover:opacity-90"
                  style={{ backgroundColor: "#dbefe7", color: "#105b39" }}
                >
                  Download CSV
                </button>
              )}
              <button
                onClick={handleCloseGuideModal}
                className="rounded-lg px-4 py-2 text-sm font-medium hover:opacity-90"
                style={{ backgroundColor: "#46896c", color: "#dbefe7" }}
              >
                {guideResult ? "Close" : "Cancel"}
              </button>
              <button
                onClick={handleSubmitGuideDesign}
                disabled={guideLoading || !sequence}
                className="rounded-lg px-4 py-2 text-sm font-medium hover:opacity-90 disabled:opacity-50"
                style={{ backgroundColor: "#dbefe7", color: "#105b39" }}
              >
                {guideLoading ? "Designing…" : "Design"}
              </button>
            </div>
          </div>
        </div>
      )}

      {/* Design pegRNA Modal (prime editing — easy_prime PE3 XGBoost port) */}
      {pegrnaModalOpen && (
        <div
          className="fixed inset-0 z-50 flex items-center justify-center"
          style={{ backgroundColor: "rgba(0, 0, 0, 0.5)" }}
          onClick={handleClosePegrnaModal}
        >
          <div
            className="rounded-2xl p-6 w-full max-w-2xl"
            style={{ backgroundColor: "#2d4a3e", maxHeight: "90vh", overflowY: "auto" }}
            onClick={(e) => e.stopPropagation()}
          >
            <h2 className="text-lg font-medium mb-2" style={{ color: "#dbefe7" }}>
              Design pegRNAs (Prime Editing)
            </h2>
            <p className="text-sm mb-4" style={{ color: "#dbefe7", opacity: 0.85 }}>
              Designs a prime-editing pegRNA for the specified edit using a pure-Python
              port of <strong>easy_prime</strong> (Li et al. 2021, Nat. Commun. 12:5121).
              Sweeps PBS length 10-15 nt and RTT length 10-20 nt over every nearby NGG
              PAM, then scores each (sgRNA, PBS, RTT, ngRNA) candidate with the bundled
              <strong> PE3 XGBoost</strong> regressor (23 features incl. 10 RNAplfold
              base-pair-probabilities for scaffold/RTT folding). Returns the top 3
              re-ranked by easy_prime's dPAM / PE3b precedence. The spacer is annotated
              on the plasmid; click the annotation for the full pegRNA component
              breakdown.
            </p>

            <div className="grid grid-cols-2 gap-3 mb-4">
              <div>
                <label className="block text-xs mb-1" style={{ color: "#dbefe7", opacity: 0.85 }}>Edit start (1-indexed)</label>
                <input type="number" min={1}
                  value={pegrnaForm.edit_start}
                  onChange={handlePegrnaFormChange("edit_start")}
                  className="rounded-lg px-3 py-1 text-sm w-full"
                  style={{ backgroundColor: "#46896c", color: "#dbefe7", border: "none" }} />
              </div>
              <div>
                <label className="block text-xs mb-1" style={{ color: "#dbefe7", opacity: 0.85 }}>Edit end (1-indexed)</label>
                <input type="number" min={1}
                  value={pegrnaForm.edit_end}
                  onChange={handlePegrnaFormChange("edit_end")}
                  className="rounded-lg px-3 py-1 text-sm w-full"
                  style={{ backgroundColor: "#46896c", color: "#dbefe7", border: "none" }} />
              </div>
              <div>
                <label className="block text-xs mb-1" style={{ color: "#dbefe7", opacity: 0.85 }}>Edit type</label>
                <select
                  value={pegrnaForm.edit_type}
                  onChange={handlePegrnaFormChange("edit_type") as any}
                  className="rounded-lg px-3 py-1 text-sm w-full"
                  style={{ backgroundColor: "#46896c", color: "#dbefe7", border: "none" }}>
                  <option value="substitution">Substitution (len(alt) == range)</option>
                  <option value="insertion">Insertion (alt inserted after range)</option>
                  <option value="deletion">Deletion (range removed; alt ignored)</option>
                </select>
              </div>
              <div>
                <label className="block text-xs mb-1" style={{ color: "#dbefe7", opacity: 0.85 }}>Top-N pegRNAs</label>
                <input type="number" min={1} max={20}
                  value={pegrnaForm.n_results}
                  onChange={handlePegrnaFormChange("n_results")}
                  className="rounded-lg px-3 py-1 text-sm w-full"
                  style={{ backgroundColor: "#46896c", color: "#dbefe7", border: "none" }} />
              </div>
              <div className="col-span-2">
                <label className="block text-xs mb-1" style={{ color: "#dbefe7", opacity: 0.85 }}>Alt (replacement / insertion) sequence</label>
                <input type="text"
                  value={pegrnaForm.alt}
                  onChange={handlePegrnaFormChange("alt")}
                  placeholder="e.g. A (sub G->A); or AGCT for an insertion; ignored for deletion"
                  className="rounded-lg px-3 py-1 text-sm w-full font-mono"
                  style={{ backgroundColor: "#46896c", color: "#dbefe7", border: "none" }} />
              </div>
              <div className="col-span-2">
                <label className="flex items-center gap-2 text-xs" style={{ color: "#dbefe7", opacity: 0.85 }}>
                  <input type="checkbox"
                    checked={pegrnaForm.use_pe3}
                    onChange={handlePegrnaFormChange("use_pe3")} />
                  Require PE3 ngRNA (recommended — also enables PE3b auto-detection on pure substitutions)
                </label>
              </div>
            </div>

            {pegrnaError && (
              <div className="mb-3 text-sm rounded-lg px-3 py-2"
                   style={{ backgroundColor: "#ff6b6b", color: "#fff" }}>
                {pegrnaError}
              </div>
            )}
            {pegrnaResult && (
              <div className="mb-3 rounded-lg px-3 py-2 text-sm"
                   style={{ backgroundColor: "#46896c", color: "#dbefe7" }}>
                <div>
                  Scanned <strong>{pegrnaResult.summary.n_sgrnas_scanned}</strong> sgRNAs ·
                  <strong> {pegrnaResult.summary.n_valid_pegRNAs}</strong> valid pegRNA spacers ·
                  scored <strong>{pegrnaResult.summary.n_candidates}</strong> candidates ·
                  added <strong>{pegrnaResult.n_added}</strong> pegRNA annotation
                  {pegrnaResult.n_added === 1 ? "" : "s"} to the plasmid.
                </div>
                <div className="text-xs mt-1" style={{ opacity: 0.85 }}>
                  Edit: {pegrnaResult.summary.edit_type} {pegrnaResult.summary.edit_ref || "(empty)"} {'2192'} {pegrnaResult.summary.edit_alt || "(empty)"}. Click any pegRNA on the plasmid for the full design (spacer + scaffold + RTT + PBS + ngRNA).
                </div>
              </div>
            )}

            <div className="flex items-center justify-end gap-2">
              <button
                onClick={handleClosePegrnaModal}
                className="rounded-lg px-4 py-2 text-sm font-medium hover:opacity-90"
                style={{ backgroundColor: "#46896c", color: "#dbefe7" }}
              >
                {pegrnaResult ? "Close" : "Cancel"}
              </button>
              <button
                onClick={handleSubmitPegrnaDesign}
                disabled={pegrnaLoading || !sequence}
                className="rounded-lg px-4 py-2 text-sm font-medium hover:opacity-90 disabled:opacity-50"
                style={{ backgroundColor: "#dbefe7", color: "#105b39" }}
              >
                {pegrnaLoading ? "Designing…" : "Design pegRNAs"}
              </button>
            </div>
          </div>
        </div>
      )}
    </div>
  );
}
