"use client";

import React, { useCallback, useEffect, useMemo, useRef, useState } from "react";
import dynamic from "next/dynamic";
import InteractiveSequenceViewer, { SeqVizAnnotation as InteractiveAnnotation, RestrictionSiteAnnotation } from "./InteractiveSequenceViewer";

// Loading animation component - Sheep running around plasmid (runs inside ring, jumps onto features)
function SplicifyLoader() {
  const sheepRef = useRef<HTMLImageElement>(null);
  const feat0Ref = useRef<SVGPathElement>(null);
  const feat1Ref = useRef<SVGPathElement>(null);
  const feat2Ref = useRef<SVGPathElement>(null);

  useEffect(() => {
    const TAU = Math.PI * 2;

    // Plasmid sizing (smaller map, sheep runs inside)
    const R_MAP = 56;
    const R_RUN = R_MAP - 12;
    const R_FEATURE = R_MAP;
    const RUN_BOB = 1.3;

    // Motion tuning
    const runSpeed = 1.35;
    const jumpDuration = 0.44;
    const jumpHeight = 22;
    const triggerWindow = 0.12;
    const sheepScale = 1.0;

    // Feature placement (3 features)
    const featureAngles = [
      (25 * Math.PI) / 180,
      (155 * Math.PI) / 180,
      (275 * Math.PI) / 180,
    ];

    const featureSpanDeg = 22;
    const featEls = [feat0Ref.current, feat1Ref.current, feat2Ref.current];

    function arcPath(radius: number, startAngle: number, endAngle: number) {
      const sx = radius * Math.cos(startAngle);
      const sy = radius * Math.sin(startAngle);
      const ex = radius * Math.cos(endAngle);
      const ey = radius * Math.sin(endAngle);
      const delta = (endAngle - startAngle) % TAU;
      const largeArc = delta > Math.PI ? 1 : 0;
      const sweep = 1;
      return `M ${sx.toFixed(3)} ${sy.toFixed(3)} A ${radius} ${radius} 0 ${largeArc} ${sweep} ${ex.toFixed(3)} ${ey.toFixed(3)}`;
    }

    featureAngles.forEach((a, i) => {
      if (!featEls[i]) return;
      const start = a - (featureSpanDeg / 2) * (Math.PI / 180);
      const end = a + (featureSpanDeg / 2) * (Math.PI / 180);
      featEls[i]!.setAttribute("d", arcPath(R_FEATURE, start, end));
    });

    let lastT = performance.now();
    let angle = 0;
    let jumping = false;
    let jumpT = 0;
    let lastTriggeredFeature = -1;
    let animationId: number;

    function angDist(a: number, b: number) {
      let d = Math.abs(a - b) % TAU;
      return d > Math.PI ? TAU - d : d;
    }

    function whichFeatureTriggered(a: number) {
      for (let i = 0; i < featureAngles.length; i++) {
        const d = angDist(a, featureAngles[i]);
        if (d < triggerWindow) return i;
      }
      return -1;
    }

    function tick(now: number) {
      const dt = Math.min(0.033, (now - lastT) / 1000);
      lastT = now;

      // Move clockwise
      angle = (angle + runSpeed * dt) % TAU;

      // Trigger jump only when near a feature
      if (!jumping) {
        const hitIdx = whichFeatureTriggered(angle);
        if (hitIdx !== -1 && hitIdx !== lastTriggeredFeature) {
          jumping = true;
          jumpT = 0;
          lastTriggeredFeature = hitIdx;
        }
      }

      // Base running radius is inside the plasmid ring
      let radial = R_RUN + Math.sin(now * 0.02) * RUN_BOB;

      // Jump curve: run inside -> jump outward onto feature -> return inside
      if (jumping) {
        jumpT += dt;
        const p = Math.min(1, jumpT / jumpDuration);
        const j = Math.sin(Math.PI * p);
        radial = R_RUN + (jumpHeight * j);

        if (p >= 1) {
          jumping = false;
          jumpT = 0;
        }
      }

      // Polar -> cartesian
      const x = radial * Math.cos(angle);
      const y = radial * Math.sin(angle);

      // Orient sheep along tangent
      const rot = (angle + Math.PI / 2) * (180 / Math.PI);

      // Small run waggle (reduced during jump)
      const waggle = Math.sin(now * 0.03) * (jumping ? 2 : 5);

      if (sheepRef.current) {
        sheepRef.current.style.transform =
          `translate(calc(-50% + ${x}px), calc(-50% + ${y}px)) ` +
          `rotate(${rot + waggle}deg) ` +
          `scale(${sheepScale})`;
      }

      animationId = requestAnimationFrame(tick);
    }

    animationId = requestAnimationFrame(tick);
    return () => cancelAnimationFrame(animationId);
  }, []);

  return (
    <div
      style={{
        width: "240px",
        height: "240px",
        position: "relative",
        margin: "0 auto",
      }}
    >
      <svg viewBox="-120 -120 240 240" aria-label="Loading" style={{ width: "100%", height: "100%", display: "block" }}>
        {/* Backbone (smaller ring) */}
        <circle r="56" fill="none" stroke="#dbefe7" strokeWidth="9" strokeLinecap="round" opacity="0.95" />
        {/* Features (3 arcs) */}
        <g stroke="#dbefe7" fill="none" strokeLinecap="round" opacity="1">
          <path ref={feat0Ref} d="" strokeWidth="13" />
          <path ref={feat1Ref} d="" strokeWidth="13" />
          <path ref={feat2Ref} d="" strokeWidth="13" />
        </g>
      </svg>
      <img
        ref={sheepRef}
        src="/splicify-logo.png"
        alt="Loading"
        style={{
          position: "absolute",
          left: "50%",
          top: "50%",
          width: "92px",
          height: "92px",
          transform: "translate(-50%, -50%)",
          transformOrigin: "50% 60%",
          userSelect: "none",
          pointerEvents: "none",
        }}
      />
    </div>
  );
}

const SeqVizAny = dynamic(async () => {
  const mod: any = await import("seqviz");
  return (mod.default ?? mod.SeqViz ?? mod) as React.ComponentType<any>;
}, { ssr: false });

type SeqVizAnnotation = {
  name: string;
  start: number;
  end: number;
  direction?: 1 | -1;
  color?: string;
  role?: string;
  origin?: string;
  source?: string;
};

type VizMeta = {
  template_index?: number;
  template_name?: string;
  amplicon_name?: string;
};

type VizPayload =
  | {
      type: "pcr" | "gibson";
      sequence: string;
      annotations: SeqVizAnnotation[];
      meta?: VizMeta;
    }
  | {
      type: "design";
      sequence: string;
      annotations: SeqVizAnnotation[];
      topology: string;
      title: string;
      total_length: number;
      restriction_sites?: RestrictionSiteAnnotation[];
    }
  | {
      type: "annotation";
      sequence: string;
      annotations: SeqVizAnnotation[];
      circular?: boolean;
      title?: string;
    }
  | null;

type ApiFile = {
  fileName?: string;
  mimeType?: string;
  dataBase64?: string;
};

type ApiResponseShape = {
  ok?: boolean;
  reply?: string;
  sessionId?: string;
  viz?: VizPayload;
  viz_list?: VizPayload[];
  files?: ApiFile[];
  response?: {
    ok?: boolean;
    reply?: string;
    sessionId?: string;
    viz?: VizPayload;
    viz_list?: VizPayload[];
    files?: ApiFile[];
  };
};

type Msg = {
  role: "user" | "assistant";
  content: string;
  files?: ApiFile[];
};

function uuidv4() {
  if (typeof crypto !== "undefined" && "randomUUID" in crypto) return crypto.randomUUID();
  return "xxxxxxxxyxxx4xxxyxxxxyxxxxxxxxxxx".replace(/[xy]/g, (c) => {
    const r = (Math.random() * 16) | 0;
    const v = c === "x" ? r : (r & 0x3) | 0x8;
    return v.toString(16);
  });
}

function unwrapN8nPayload(x: any) {
  if (Array.isArray(x)) return x[0] ?? null;
  return x ?? null;
}

function computeAmpliconLength(viz: VizPayload): number | null {
  if (!viz || viz.type !== "pcr") return null;
  const anns = viz.annotations ?? [];
  const left = anns.find((a) => (a.name ?? "").toLowerCase().includes("left"));
  const right = anns.find((a) => (a.name ?? "").toLowerCase().includes("right"));
  if (!left || !right) return null;

  const start = Math.min(left.start, left.end, right.start, right.end);
  const end = Math.max(left.start, left.end, right.start, right.end);
  const len = Math.max(0, end - start);
  return len > 0 ? len : null;
}

function normalizeAnnotations(annotations: SeqVizAnnotation[]) {
  return (annotations ?? []).map((a: any) => {
    const start = Number(a.start ?? 0);
    const end = Number(a.end ?? start);
    const length = Math.max(0, end - start);
    const strand = a.direction === -1 ? -1 : 1;

    return {
      name: a.name ?? "Annotation",
      start,
      end,
      length,
      strand,
      direction: strand,
      ...(a.color ? { color: a.color } : {}),
    };
  });
}

/**
 * Merge annotations that cross the origin in circular sequences.
 * When an annotation crosses the origin, it may be split into two parts:
 * one ending near the sequence end and one starting near the beginning.
 * This function merges them back into a single annotation with start > end.
 */
function mergeOriginCrossingAnnotations(
  annotations: ReturnType<typeof normalizeAnnotations>,
  seqLength: number
): ReturnType<typeof normalizeAnnotations> {
  if (!annotations || annotations.length === 0 || seqLength <= 0) {
    return annotations;
  }

  // Group annotations by name
  const byName: Record<string, typeof annotations> = {};
  for (const ann of annotations) {
    const key = ann.name;
    if (!byName[key]) byName[key] = [];
    byName[key].push(ann);
  }

  const result: typeof annotations = [];

  for (const name of Object.keys(byName)) {
    const group = byName[name];

    if (group.length === 2) {
      // Check if these two annotations might be parts of one crossing the origin
      const [a, b] = group;

      // One should end near the sequence end, other should start near the beginning
      const aEndsNearEnd = a.end >= seqLength * 0.8;
      const aStartsNearStart = a.start <= seqLength * 0.2;
      const bEndsNearEnd = b.end >= seqLength * 0.8;
      const bStartsNearStart = b.start <= seqLength * 0.2;

      // Check if they have the same strand/direction
      const sameStrand = a.strand === b.strand;

      if (sameStrand) {
        if (aEndsNearEnd && bStartsNearStart && b.start < a.start) {
          // a is the "end" part, b is the "start" part - merge them
          // For origin-crossing: start is where the annotation begins (in a), end is where it ends (in b)
          result.push({
            name: a.name,
            start: a.start,
            end: b.end,
            length: (seqLength - a.start) + b.end,
            strand: a.strand,
            direction: a.direction,
          });
          continue;
        } else if (bEndsNearEnd && aStartsNearStart && a.start < b.start) {
          // b is the "end" part, a is the "start" part - merge them
          result.push({
            name: b.name,
            start: b.start,
            end: a.end,
            length: (seqLength - b.start) + a.end,
            strand: b.strand,
            direction: b.direction,
          });
          continue;
        }
      }
    }

    // No merge - add all annotations from this group
    result.push(...group);
  }

  return result;
}

function pcrTabLabel(v: VizPayload, idx: number) {
  const meta = v && typeof v === "object" ? (v as any).meta : null;
  return meta?.amplicon_name || meta?.template_name || `PCR_${idx + 1}`;
}

function base64ToBytes(b64: string): Uint8Array {
  const cleaned = String(b64).replace(/^data:.*;base64,/, "");
  const binStr = atob(cleaned);
  const bytes = new Uint8Array(binStr.length);
  for (let i = 0; i < binStr.length; i++) bytes[i] = binStr.charCodeAt(i);
  return bytes;
}

function fileToObjectUrl(file: ApiFile): string | null {
  const b64 = file?.dataBase64;
  if (!b64 || typeof b64 !== "string" || b64.length < 8) return null;

  try {
    const bytes = base64ToBytes(b64);
    const ab = new ArrayBuffer(bytes.byteLength);
    new Uint8Array(ab).set(bytes);

    const blob = new Blob([ab], {
      type: file.mimeType || "application/octet-stream",
    });

    return URL.createObjectURL(blob);
  } catch {
    return null;
  }
}

function formatBytes(n: number) {
  if (!Number.isFinite(n) || n <= 0) return "0 B";
  const units = ["B", "KB", "MB", "GB"];
  let i = 0;
  let x = n;
  while (x >= 1024 && i < units.length - 1) {
    x /= 1024;
    i++;
  }
  return `${x.toFixed(i === 0 ? 0 : 1)} ${units[i]}`;
}

function isGenbankFile(f: File) {
  const name = (f?.name || "").toLowerCase();
  return name.endsWith(".gb") || name.endsWith(".gbk") || name.endsWith(".genbank");
}

// Helper to detect and format DNA sequences in message content
function formatMessageContent(content: string): React.ReactNode {
  // Match DNA sequences: 20+ consecutive A, T, G, C, N (case insensitive)
  const dnaPattern = /([ATGCNatgcn]{20,})/g;

  const parts = content.split(dnaPattern);
  if (parts.length === 1) {
    // No DNA sequences found, return as-is
    return content;
  }

  return parts.map((part, idx) => {
    if (dnaPattern.test(part)) {
      // Reset lastIndex since we're reusing the regex
      dnaPattern.lastIndex = 0;
      return (
        <span key={idx} className="dna-sequence-inline">
          {part}
        </span>
      );
    }
    return part;
  });
}

export default function Chat() {
  const [input, setInput] = useState("");
  const [messages, setMessages] = useState<Msg[]>([]);
  const [sessionId, setSessionId] = useState<string>("");
  const [viz, setViz] = useState<VizPayload>(null);
  const [vizList, setVizList] = useState<VizPayload[]>([]);
  const [activePcrIdx, setActivePcrIdx] = useState<number>(0);
  const [designVizList, setDesignVizList] = useState<VizPayload[]>([]);
  const [activeDesignIdx, setActiveDesignIdx] = useState<number>(0);
  const [loading, setLoading] = useState(false);

  const [uploadedGb, setUploadedGb] = useState<File | null>(null);
  const [uploadedInventoryGbs, setUploadedInventoryGbs] = useState<File[]>([]);
  const [dragActive, setDragActive] = useState(false);
  const [uploadError, setUploadError] = useState<string>("");
  const [inventoryUploadError, setInventoryUploadError] = useState<string>("");

  // Track the current plasmid name (from uploaded file, without extension)
  const [currentPlasmidName, setCurrentPlasmidName] = useState<string | null>(null);

  // Workflow options
  const [includeAiExplanation, setIncludeAiExplanation] = useState(true);
  const describePlasmidIntent = false;

  const [fileUrlMap, setFileUrlMap] = useState<Record<string, string>>({});

  const chatEndRef = useRef<HTMLDivElement>(null);
  const textareaRef = useRef<HTMLTextAreaElement>(null);

  const examples = useMemo(
    () => [
      "Design a Gibson assembly to build my target plasmid from inventory",
      "Create PCR primers to amplify and modify my templates",
      "Clone my gene into a vector using Gateway recombination",
      "Generate sgRNA oligos for Golden Gate cloning into lentiCRISPR",
      "Design site-directed mutagenesis primers to modify my plasmid",
      "Plan a restriction digest cloning strategy for my insert",
    ],
    []
  );

  // Example prompts for quick-action buttons
  const examplePrompts = useMemo(
    () => ({
      plasmidAnnotation: "Annotate this plasmid",
      gibsonFragments: `Design gibson assembly primers for assembling a plasmid with these fragments: Frag1: GCCTCCTGCTGGTCCCAAGTTGTGAAATCTTTATCGTGTTTGGTCAGTTCCAGGCGATGTTCAACGAAATGAAATTCTGGCAGCTTGACAGGTTTTTTGGCTTTATACGTCGTTTTGAAATCAACACGATGATGACCGCCTCCTTTCA, Frag2: CTTGCCTTCACCTTTACCCTCGATCGTAAAATCATGGCCATTGACAGTGCCCTCAAGGTGCAACTTGGTCTTCATTACCTGCTTAATCACTGACATAGATCCTTTCTCCTCTTTAGATCTTTTGAATTCACTAGTATTATACCTAGGACTGAGCTAGCTGTCAAGCGCAACGCAATTAATGTAAGTTAGCTCACTCATTAGGCACCGACGTCAGGTGGCACTTTT, Frag3: TCCACTATTCGAGGCCGTTCGTTAATACTTGTTGCGTTCCTAGCCGCTATATTTGTCTCTTTGCCGACTAATGTGGACAAGCACACCATAGCCATTTATCGGAGCGCCTCGGAATACGGTATGAGCAGGCGCCTCGTGAGACCATTGCGA, frag4: GCCACCATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAAAGCGGCCGCACTCCTCAGG`,
      restrictionWorkflow: "Run the test restriction enzyme workflow for cloning GFP into pUC19.",
      sdmWorkflow: "delete the His-tag from the test plasmid.",
      sgRNAGoldenGate: "Design oligos to clone gRNA GAGTCCGAGCAGAAGAAGAA (EMX1) into lentiCRISPR v2 using Golden Gate assembly.",
      goldenGateAssembly: "Design Golden Gate primers to assemble EF1a promoter + eGFP + bGH polyA terminator.",
      gatewayCloning: "Design primers for insertion of GFP into pDONR221 using Gateway BP recombination.",
    }),
    []
  );

  const [exampleIdx, setExampleIdx] = useState(0);
  const [prevIdx, setPrevIdx] = useState<number | null>(null);
  const [phase, setPhase] = useState<"idle" | "animating">("idle");

  // Auto-scroll to bottom when messages change
  useEffect(() => {
    chatEndRef.current?.scrollIntoView({ behavior: "smooth" });
  }, [messages]);

  // Auto-resize textarea up to 6 lines
  const adjustTextareaHeight = useCallback(() => {
    const textarea = textareaRef.current;
    if (!textarea) return;

    // Reset height to auto to get the actual scrollHeight
    textarea.style.height = "auto";

    // Calculate line height (approximately 1.5rem = 24px for text-base)
    const lineHeight = 24;
    const maxLines = 6;
    const maxHeight = lineHeight * maxLines;
    const minHeight = lineHeight; // Single line minimum

    // Set height based on content, capped at max
    const newHeight = Math.min(Math.max(textarea.scrollHeight, minHeight), maxHeight);
    textarea.style.height = `${newHeight}px`;

    // Enable/disable scrolling based on content
    textarea.style.overflowY = textarea.scrollHeight > maxHeight ? "auto" : "hidden";
  }, []);

  // Adjust height when input changes
  useEffect(() => {
    adjustTextareaHeight();
  }, [input, adjustTextareaHeight]);

  useEffect(() => {
    const interval = window.setInterval(() => {
      setPrevIdx(exampleIdx);
      setPhase("animating");

      window.setTimeout(() => {
        setExampleIdx((i) => (i + 1) % examples.length);
        setPhase("idle");
        setPrevIdx(null);
      }, 420);
    }, 5000);

    return () => window.clearInterval(interval);
  }, [exampleIdx, examples.length]);

  useEffect(() => {
    const key = "splicify_session_id";
    const existing = window.localStorage.getItem(key);
    if (existing) {
      setSessionId(existing);
    } else {
      const fresh = uuidv4();
      window.localStorage.setItem(key, fresh);
      setSessionId(fresh);
    }
  }, []);

  useEffect(() => {
    return () => {
      Object.values(fileUrlMap).forEach((u) => {
        try {
          URL.revokeObjectURL(u);
        } catch {}
      });
    };
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  const canSend = useMemo(
    () => !loading && (input.trim().length > 0 || !!uploadedGb),
    [loading, input, uploadedGb]
  );

  const activePcrViz: VizPayload = useMemo(() => {
    const pcrs = (vizList || []).filter((v) => v && v.type === "pcr") as VizPayload[];
    if (pcrs.length > 0) {
      const idx = Math.max(0, Math.min(activePcrIdx, pcrs.length - 1));
      return pcrs[idx] ?? null;
    }
    return viz?.type === "pcr" ? viz : null;
  }, [viz, vizList, activePcrIdx]);

  const pcrAmpliconLen = computeAmpliconLength(activePcrViz);
  const pcrSeq = activePcrViz?.type === "pcr" ? activePcrViz.sequence : "";
  const pcrAnn = activePcrViz?.type === "pcr" ? normalizeAnnotations(activePcrViz.annotations ?? []) : [];

  const gibsonSeq = viz?.type === "gibson" ? viz.sequence : "";
  const gibsonAnnNormalized = viz?.type === "gibson" ? normalizeAnnotations(viz.annotations ?? []) : [];
  const gibsonAnnMerged = viz?.type === "gibson"
    ? mergeOriginCrossingAnnotations(gibsonAnnNormalized, gibsonSeq.length)
    : [];

  // Remove direction from overlap annotations so they display without arrows
  const gibsonAnn = gibsonAnnMerged.map((a) => {
    const name = (a.name ?? "").toLowerCase();
    if (name.includes("overlap") || name.includes("homology")) {
      return { ...a, direction: 0, strand: 0 };
    }
    return a;
  });

  const onPickFile = useCallback((file: File | null) => {
    setUploadError("");
    if (!file) return;

    if (!isGenbankFile(file)) {
      setUploadError("Please upload a GenBank file ending in .gb, .gbk, or .genbank.");
      return;
    }

    const maxBytes = 10 * 1024 * 1024;
    if (file.size > maxBytes) {
      setUploadError(`That file is ${formatBytes(file.size)}. Please upload a GenBank file under 10 MB.`);
      return;
    }

    setUploadedGb(file);

    // Extract plasmid name from filename (without extension)
    const fileName = file.name;
    const nameWithoutExt = fileName.replace(/\.(gb|gbk|genbank)$/i, "");
    setCurrentPlasmidName(nameWithoutExt);
  }, []);

  // Inventory file handler - accepts multiple files
  const onPickInventoryFiles = useCallback((files: FileList | File[] | null) => {
    setInventoryUploadError("");
    if (!files) return;

    const arr: File[] = Array.isArray(files) ? files : Array.from(files);
    if (!arr.length) return;

    const maxBytes = 10 * 1024 * 1024;

    setUploadedInventoryGbs((prev) => {
      let next = [...prev];

      for (const file of arr) {
        if (!file) continue;

        if (!isGenbankFile(file)) {
          setInventoryUploadError("Inventory files must be GenBank (.gb, .gbk, .genbank).");
          continue;
        }

        if (file.size > maxBytes) {
          setInventoryUploadError(`One inventory file is ${formatBytes(file.size)}. Max 10 MB each.`);
          continue;
        }

        if (next.length >= 10) {
          setInventoryUploadError("Maximum of 10 inventory plasmids per session.");
          break;
        }

        // Avoid duplicates by name+size
        const sig = `${file.name}::${file.size}`;
        const exists = next.some((f) => `${f.name}::${f.size}` === sig);
        if (exists) continue;

        next.push(file);
      }

      return next;
    });
  }, []);

  const removeInventoryFile = useCallback((idx: number) => {
    setUploadedInventoryGbs((prev) => prev.filter((_, i) => i !== idx));
  }, []);

  const clearInventory = useCallback(() => {
    setUploadedInventoryGbs([]);
    setInventoryUploadError("");
  }, []);

  const onDrop = useCallback(
    (e: React.DragEvent) => {
      e.preventDefault();
      e.stopPropagation();
      setDragActive(false);

      const f = e.dataTransfer?.files?.[0] ?? null;
      if (f) onPickFile(f);
    },
    [onPickFile]
  );

  const onDragOver = useCallback((e: React.DragEvent) => {
    e.preventDefault();
    e.stopPropagation();
    setDragActive(true);
  }, []);

  const onDragLeave = useCallback((e: React.DragEvent) => {
    e.preventDefault();
    e.stopPropagation();
    setDragActive(false);
  }, []);

  // Agent V2 toggle: when ON, send() posts to /agent_v2/chat (Claude
  // Sonnet 4.6 tool-use loop). Default is OFF — the classic /api/chat
  // intent dispatcher handles the request. V2 is opt-in until the
  // production cutover decision is made (see agent_v2_summary.md).
  const [useAgentV2, setUseAgentV2] = useState<boolean>(false);

  async function send() {
    const text = input.trim();
    if (loading) return;
    if (!text && !uploadedGb) return;

    setInput("");
    setLoading(true);

    // Build user message content
    const uploadedFiles: string[] = [];
    if (uploadedGb) uploadedFiles.push(`target: ${uploadedGb.name}`);
    if (uploadedInventoryGbs.length > 0) {
      uploadedFiles.push(`inventory: ${uploadedInventoryGbs.map(f => f.name).join(", ")}`);
    }
    const filesSummary = uploadedFiles.length > 0 ? ` [${uploadedFiles.join("; ")}]` : "";
    const messageContent = text ? `${text}${filesSummary}` : (filesSummary ? filesSummary.trim() : "");

    setMessages((m) => [
      ...m,
      {
        role: "user",
        content: messageContent || "(uploaded files)",
      },
    ]);

    // Agent V2 mode (opt-in — when useAgentV2 is ON): short-circuit to
    // /agent_v2/chat (Claude Sonnet 4.6 tool-use loop). Always sends
    // multipart so any uploaded .gb files come along; the Vercel proxy at
    // app/api/agent_v2/chat/route.ts forwards to the VPS backend.
    if (useAgentV2) {
      try {
        const apiBase = process.env.NEXT_PUBLIC_API_BASE_URL || "";
        const form = new FormData();
        form.append("message", text);
        form.append("session_id", sessionId);
        if (uploadedGb) form.append("file", uploadedGb, uploadedGb.name);
        for (const inv of uploadedInventoryGbs) {
          form.append("inventory_files", inv, inv.name);
        }
        const agentRes = await fetch(apiBase + "/api/agent_v2/chat", {
          method: "POST",
          body: form,
        });
        const ct = agentRes.headers.get("content-type") || "";
        const raw = await agentRes.text();
        let data: any = null;
        if (ct.includes("application/json")) {
          try { data = JSON.parse(raw); } catch { /* fall through */ }
        }
        if (!data) {
          const preview = raw.slice(0, 200).replace(/<[^>]+>/g, "").trim();
          throw new Error(`Non-JSON response (HTTP ${agentRes.status}): ${preview || "(empty)"}`);
        }
        const reply = (data && data.reply) || "(no reply)";
        const traceNote = data && typeof data.n_tool_calls === "number"
          ? `\n\n_(${data.n_tool_calls} tool call${data.n_tool_calls === 1 ? "" : "s"})_`
          : "";
        // Surface files + viz exactly like the regular /api/chat path so
        // assembled.gb downloads and CircularPlasmidViewer render.
        const incomingFiles = Array.isArray((data as any)?.files) ? (data as any).files : [];
        if (incomingFiles.length > 0) {
          setFileUrlMap((prev) => {
            const next = { ...prev };
            for (const f of incomingFiles) {
              const fn = f?.fileName || "";
              if (!fn) continue;
              if (next[fn]) {
                try { URL.revokeObjectURL(next[fn]); } catch {}
                delete next[fn];
              }
              const url = fileToObjectUrl(f);
              if (url) next[fn] = url;
            }
            return next;
          });
        }
        setMessages((m) => [
          ...m,
          { role: "assistant", content: reply + traceNote,
            files: incomingFiles.length ? incomingFiles : undefined },
        ]);
        const incomingViz = (data as any)?.viz ?? null;
        if (incomingViz && (incomingViz as any).sequence) {
          setViz(incomingViz);
          setVizList([]);
          setActivePcrIdx(0);
          setDesignVizList([]);
          setActiveDesignIdx(0);
        } else {
          // V2's chat envelope returns files but no viz. If assembled.gb came
          // back, hydrate the viewer by POSTing it to /agent_v2/annotate-on-upload
          // (same endpoint Chat already uses for eager upload annotation).
          const assembled = incomingFiles.find((f: any) => f?.fileName === "assembled.gb");
          if (assembled && assembled.dataBase64) {
            try {
              const bin = atob(assembled.dataBase64);
              const bytes = new Uint8Array(bin.length);
              for (let i = 0; i < bin.length; i++) bytes[i] = bin.charCodeAt(i);
              const fileObj = new File([bytes], "assembled.gb", { type: "chemical/seq-na-genbank" });
              const annForm = new FormData();
              annForm.append("file", fileObj, "assembled.gb");
              const annRes = await fetch(apiBase + "/agent_v2/annotate-on-upload", {
                method: "POST",
                body: annForm,
              });
              if (annRes.ok) {
                const annData = await annRes.json();
                const annViz = (annData as any)?.viz ?? annData;
                if (annViz && (annViz as any).sequence) {
                  setViz(annViz);
                  setVizList([]);
                  setActivePcrIdx(0);
                  setDesignVizList([]);
                  setActiveDesignIdx(0);
                }
              }
            } catch {
              // viewer hydration is best-effort — agent reply still renders
            }
          }
        }
      } catch (err: any) {
        const msg = err && err.message ? err.message : "agent request failed";
        setMessages((m) => [...m, { role: "assistant", content: "Error: " + msg }]);
      } finally {
        setLoading(false);
      }
      return;
    }

    try {
      let res: Response;

      // Use FormData when there are files (target or inventory)
      if (uploadedGb || uploadedInventoryGbs.length > 0) {
        const form = new FormData();
        form.append("message", text);
        form.append("session_id", sessionId);
        form.append("include_ai_explanation", String(includeAiExplanation));
        form.append("describe_plasmid_intent", String(describePlasmidIntent));

        if (uploadedGb) {
          form.append("file", uploadedGb, uploadedGb.name);
        }

        // Add inventory files with the key "inventory_files" (repeated for each file)
        for (const inv of uploadedInventoryGbs) {
          form.append("inventory_files", inv, inv.name);
        }

        res = await fetch((process.env.NEXT_PUBLIC_API_BASE_URL || "") + "/api/chat", {
          method: "POST",
          body: form,
        });
      } else {
        res = await fetch((process.env.NEXT_PUBLIC_API_BASE_URL || "") + "/api/chat", {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({
            message: text,
            session_id: sessionId,
            include_ai_explanation: includeAiExplanation,
            describe_plasmid_intent: describePlasmidIntent,
          }),
        });
      }

      const dataRawAny = await res.json().catch(() => null);
      const dataRaw = unwrapN8nPayload(dataRawAny) as ApiResponseShape | null;

      if (!res.ok) {
        const errMsg =
          (dataRaw as any)?.error ||
          (dataRaw as any)?.message ||
          `Request failed with status ${res.status}`;
        setMessages((m) => [...m, { role: "assistant", content: `Error: ${errMsg}` }]);
        setViz(null);
        setVizList([]);
        setActivePcrIdx(0);
        return;
      }

      const replyText =
        (dataRaw as any)?.response?.reply ??
        (dataRaw as any)?.reply ??
        (dataRaw as any)?.output ??
        (dataRaw as any)?.message ??
        "No reply returned.";

      const incomingFiles =
        (Array.isArray((dataRaw as any)?.response?.files) ? (dataRaw as any).response.files : null) ??
        (Array.isArray((dataRaw as any)?.files) ? (dataRaw as any).files : null) ??
        [];

      if (incomingFiles.length > 0) {
        setFileUrlMap((prev) => {
          const next = { ...prev };
          for (const f of incomingFiles) {
            const fn = f?.fileName || "";
            if (!fn) continue;
            if (next[fn]) {
              try { URL.revokeObjectURL(next[fn]); } catch {}
              delete next[fn];
            }
            const url = fileToObjectUrl(f);
            if (url) next[fn] = url;
          }
          return next;
        });
      }

      setMessages((m) => [
        ...m,
        { role: "assistant", content: replyText, files: incomingFiles.length ? incomingFiles : undefined },
      ]);

      const newSessionId =
        typeof (dataRaw as any)?.response?.sessionId === "string" && (dataRaw as any).response.sessionId.length > 0
          ? (dataRaw as any).response.sessionId
          : typeof (dataRaw as any)?.sessionId === "string" && (dataRaw as any).sessionId.length > 0
            ? (dataRaw as any).sessionId
            : "";

      if (newSessionId) {
        setSessionId(newSessionId);
        window.localStorage.setItem("splicify_session_id", newSessionId);
      }

      const incomingViz = (dataRaw as any)?.response?.viz ?? (dataRaw as any)?.viz ?? null;
      const incomingVizList = (dataRaw as any)?.response?.viz_list ?? (dataRaw as any)?.viz_list ?? null;

      if (Array.isArray(incomingVizList) && incomingVizList.length > 0) {
        const designVizzes = incomingVizList.filter((v: any) => v?.type === "design");
        const pcrVizzes = incomingVizList.filter((v: any) => v?.type !== "design");
        if (designVizzes.length > 0) {
          setDesignVizList(designVizzes);
          setActiveDesignIdx(0);
          setVizList([]);
        } else {
          setDesignVizList([]);
          setVizList(pcrVizzes);
          setActivePcrIdx(0);
        }
        if (incomingViz && (incomingViz as any).sequence) setViz(incomingViz);
        else setViz(incomingVizList[0] ?? null);
      } else {
        setVizList([]);
        setActivePcrIdx(0);
        setDesignVizList([]);
        setActiveDesignIdx(0);

        if (incomingViz && (incomingViz as any).sequence) setViz(incomingViz);
        else setViz(null);
      }

      setUploadedGb(null);
      setUploadedInventoryGbs([]);
      setUploadError("");
      setInventoryUploadError("");
    } catch (e: any) {
      setMessages((m) => [...m, { role: "assistant", content: `Error: ${e?.message ?? "Unknown error"}` }]);
      setViz(null);
      setVizList([]);
      setActivePcrIdx(0);
      setDesignVizList([]);
      setActiveDesignIdx(0);
    } finally {
      setLoading(false);
    }
  }

  // Send with explicit files - bypasses state to avoid timing issues
  async function sendWithFiles(
    message: string,
    targetFile: File | null,
    inventoryFiles: File[]
  ) {
    if (loading) return;

    setInput("");
    setLoading(true);

    // Build user message content
    const uploadedFilesList: string[] = [];
    if (targetFile) uploadedFilesList.push(`target: ${targetFile.name}`);
    if (inventoryFiles.length > 0) {
      uploadedFilesList.push(`inventory: ${inventoryFiles.map(f => f.name).join(", ")}`);
    }
    const filesSummary = uploadedFilesList.length > 0 ? ` [${uploadedFilesList.join("; ")}]` : "";
    const messageContent = message ? `${message}${filesSummary}` : (filesSummary ? filesSummary.trim() : "");

    setMessages((m) => [
      ...m,
      {
        role: "user",
        content: messageContent || "(uploaded files)",
      },
    ]);

    // Update UI to show files (for visual feedback)
    if (targetFile) setUploadedGb(targetFile);
    if (inventoryFiles.length > 0) setUploadedInventoryGbs(inventoryFiles);

    try {
      let res: Response;

      if (targetFile || inventoryFiles.length > 0) {
        const form = new FormData();
        form.append("message", message);
        form.append("session_id", sessionId);
        form.append("include_ai_explanation", String(includeAiExplanation));
        form.append("describe_plasmid_intent", String(describePlasmidIntent));

        if (targetFile) {
          form.append("file", targetFile, targetFile.name);
        }

        for (const inv of inventoryFiles) {
          form.append("inventory_files", inv, inv.name);
        }

        res = await fetch((process.env.NEXT_PUBLIC_API_BASE_URL || "") + "/api/chat", {
          method: "POST",
          body: form,
        });
      } else {
        res = await fetch((process.env.NEXT_PUBLIC_API_BASE_URL || "") + "/api/chat", {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({
            message: message,
            session_id: sessionId,
            include_ai_explanation: includeAiExplanation,
            describe_plasmid_intent: describePlasmidIntent,
          }),
        });
      }

      const dataRawAny = await res.json().catch(() => null);
      const dataRaw = unwrapN8nPayload(dataRawAny) as ApiResponseShape | null;

      if (!res.ok) {
        const errMsg =
          (dataRaw as any)?.error ||
          (dataRaw as any)?.message ||
          `Request failed with status ${res.status}`;
        setMessages((m) => [...m, { role: "assistant", content: `Error: ${errMsg}` }]);
        setViz(null);
        setVizList([]);
        setActivePcrIdx(0);
        return;
      }

      const replyText =
        (dataRaw as any)?.response?.reply ??
        (dataRaw as any)?.reply ??
        (dataRaw as any)?.output ??
        (dataRaw as any)?.message ??
        "No reply returned.";

      const incomingFiles =
        (Array.isArray((dataRaw as any)?.response?.files) ? (dataRaw as any).response.files : null) ??
        (Array.isArray((dataRaw as any)?.files) ? (dataRaw as any).files : null) ??
        [];

      if (incomingFiles.length > 0) {
        setFileUrlMap((prev) => {
          const next = { ...prev };
          for (const f of incomingFiles) {
            const fn = f?.fileName || "";
            if (!fn) continue;
            if (next[fn]) {
              try { URL.revokeObjectURL(next[fn]); } catch {}
              delete next[fn];
            }
            const url = fileToObjectUrl(f);
            if (url) next[fn] = url;
          }
          return next;
        });
      }

      setMessages((m) => [
        ...m,
        { role: "assistant", content: replyText, files: incomingFiles.length ? incomingFiles : undefined },
      ]);

      const newSessionId =
        typeof (dataRaw as any)?.response?.sessionId === "string" && (dataRaw as any).response.sessionId.length > 0
          ? (dataRaw as any).response.sessionId
          : typeof (dataRaw as any)?.sessionId === "string" && (dataRaw as any).sessionId.length > 0
            ? (dataRaw as any).sessionId
            : "";

      if (newSessionId) {
        setSessionId(newSessionId);
        window.localStorage.setItem("splicify_session_id", newSessionId);
      }

      const incomingViz = (dataRaw as any)?.response?.viz ?? (dataRaw as any)?.viz ?? null;
      const incomingVizList = (dataRaw as any)?.response?.viz_list ?? (dataRaw as any)?.viz_list ?? null;

      if (Array.isArray(incomingVizList) && incomingVizList.length > 0) {
        const designVizzes = incomingVizList.filter((v: any) => v?.type === "design");
        const pcrVizzes = incomingVizList.filter((v: any) => v?.type !== "design");
        if (designVizzes.length > 0) {
          setDesignVizList(designVizzes);
          setActiveDesignIdx(0);
          setVizList([]);
        } else {
          setDesignVizList([]);
          setVizList(pcrVizzes);
          setActivePcrIdx(0);
        }
        if (incomingViz && (incomingViz as any).sequence) setViz(incomingViz);
        else setViz(incomingVizList[0] ?? null);
      } else {
        setVizList([]);
        setActivePcrIdx(0);
        setDesignVizList([]);
        setActiveDesignIdx(0);

        if (incomingViz && (incomingViz as any).sequence) setViz(incomingViz);
        else setViz(null);
      }

      setUploadedGb(null);
      setUploadedInventoryGbs([]);
      setUploadError("");
      setInventoryUploadError("");
    } catch (e: any) {
      setMessages((m) => [...m, { role: "assistant", content: `Error: ${e?.message ?? "Unknown error"}` }]);
      setViz(null);
      setVizList([]);
      setActivePcrIdx(0);
      setDesignVizList([]);
      setActiveDesignIdx(0);
    } finally {
      setLoading(false);
    }
  }

  async function handleExampleClick(exampleType: "plasmidAnnotation" | "gibsonFragments" | "restrictionWorkflow" | "sdmWorkflow" | "sgRNAGoldenGate" | "goldenGateAssembly" | "gatewayCloning") {
    if (loading) return;

    if (exampleType === "restrictionWorkflow") {
      try {
        const fileRes = await fetch("/pUC19_empty.gb");
        if (!fileRes.ok) throw new Error(`Failed to fetch file: ${fileRes.status}`);
        const blob = await fileRes.blob();
        const vectorFile = new File([blob], "pUC19_empty.gb", { type: "application/octet-stream" });
        await sendWithFiles(examplePrompts.restrictionWorkflow, vectorFile, []);
      } catch (err) {
        console.error("[restrictionWorkflow] Failed to load sample file:", err);
        setInput(examplePrompts.restrictionWorkflow);
      }
      return;
    }

    if (exampleType === "sdmWorkflow") {
      try {
        const fileRes = await fetch("/pCMV_6xHis_GFP_demo.gb");
        if (!fileRes.ok) throw new Error(`Failed to fetch file: ${fileRes.status}`);
        const blob = await fileRes.blob();
        const targetFile = new File([blob], "pCMV_6xHis_GFP_demo.gb", { type: "application/octet-stream" });
        await sendWithFiles(examplePrompts.sdmWorkflow, targetFile, []);
      } catch (err) {
        console.error("[sdmWorkflow] Failed to load sample file:", err);
        setInput(examplePrompts.sdmWorkflow);
      }
      return;
    }

    if (exampleType === "sgRNAGoldenGate") {
      try {
        const fileRes = await fetch("/lentiCRISPR_v2_unannotated.gb");
        if (!fileRes.ok) throw new Error(`Failed to fetch file: ${fileRes.status}`);
        const blob = await fileRes.blob();
        const targetFile = new File([blob], "lentiCRISPR_v2_unannotated.gb", { type: "application/octet-stream" });
        await sendWithFiles(examplePrompts.sgRNAGoldenGate, targetFile, []);
      } catch (err) {
        console.error("[sgRNAGoldenGate] Failed to load sample file:", err);
        setInput(examplePrompts.sgRNAGoldenGate);
      }
      return;
    }

    if (exampleType === "gatewayCloning") {
      try {
        const fileRes = await fetch("/pDONR221.gb");
        if (!fileRes.ok) throw new Error(`Failed to fetch file: ${fileRes.status}`);
        const blob = await fileRes.blob();
        const invFile = new File([blob], "pDONR221.gb", { type: "application/octet-stream" });
        await sendWithFiles(examplePrompts.gatewayCloning, null, [invFile]);
      } catch (err) {
        console.error("[gatewayCloning] Failed to load sample file:", err);
        setInput(examplePrompts.gatewayCloning);
      }
      return;
    }

    if (exampleType === "plasmidAnnotation") {
      // Bypass /api/chat: load the sample .gb client-side and let the
      // viewer's own /api/plannotate-llm flow run via autoAnnotateOnMount.
      // That path returns and renders interactions; the chat path drops them.
      try {
        const fileRes = await fetch("/lentiCRISPR_v2_unannotated.gb");
        if (!fileRes.ok) throw new Error(`Failed to fetch file: ${fileRes.status}`);
        const gbText = await fileRes.text();
        const locusLine = gbText.split(/\r?\n/, 1)[0] || "";
        const isCircular = /\bcircular\b/i.test(locusLine);
        const originIdx = gbText.indexOf("ORIGIN");
        const tail = originIdx >= 0 ? gbText.slice(originIdx) : "";
        const sequence = tail
          .replace(/^ORIGIN[^\n]*\n/, "")
          .replace(/\/\/[\s\S]*$/, "")
          .replace(/[^A-Za-z]/g, "")
          .toUpperCase();
        if (!sequence) throw new Error("Could not extract sequence from GenBank file");
        setCurrentPlasmidName("lentiCRISPR_v2_unannotated");
        setMessages((m) => [
          ...m,
          { role: "user", content: examplePrompts.plasmidAnnotation },
          { role: "assistant", content: "Annotating lentiCRISPR_v2_unannotated.gb…" },
        ]);
        setViz({
          type: "annotation",
          sequence,
          annotations: [],
          circular: isCircular,
          title: "lentiCRISPR_v2_unannotated.gb",
          auto_annotate: true,
        } as any);
        setVizList([]);
      } catch (err) {
        console.error("[plasmidAnnotation] Error:", err);
        setInput(examplePrompts.plasmidAnnotation);
      }
      return;
    }

    // For text-only examples, send directly
    await sendWithFiles(examplePrompts[exampleType], null, []);
  }

  // Handler for plasmid analysis completion - adds result to chat
  const handleAnalysisComplete = useCallback((analysis: string, moduleGraph: any) => {
    // Add the analysis as an assistant message in the chat
    setMessages((prev) => [
      ...prev,
      {
        role: "assistant" as const,
        content: analysis,
      },
    ]);
  }, []);

  const vizHeight = 500;
  const pcrTabs = (vizList || []).filter((v) => v && v.type === "pcr") as VizPayload[];
  const hasMessages = messages.length > 0;

  // Splicify Preview data - initial visualization
  const previewSequence = "GCCTCCTGCTGGTCCCAAGTTGTGAAATCTTTATCGTGTTTGGTCAGTTCCAGGCGATGTTCAACGAAATGAAATTCTGGCAGCTTGACAGGTTTTTTGGCTTTATACGTCGTTTTGAAATCAACACGATGATGACCGCCTCCTTTCACTTGCCTTCACCTTTACCCTCGATCGTAAAATCATGGCCATTGACAGTGCCCTCAAGGTGCAACTTGGTCTTCATTACCTGCTTAATCACTGACATAGATCCTTTCTCCTCTTTAGATCTTTTGAATTCACTAGTATTATACCTAGGACTGAGCTAGCTGTCAAGCGCAACGCAATTAATGTAAGTTAGCTCACTCATTAGGCACCGACGTCAGGTGGCACTTTTTCCACTATTCGAGGCCGTTCGTTAATACTTGTTGCGTTCCTAGCCGCTATATTTGTCTCTTTGCCGACTAATGTGGACAAGCACACCATAGCCATTTATCGGAGCGCCTCGGAATACGGTATGAGCAGGCGCCTCGTGAGACCATTGCGA";
  const previewAnnotationsRaw: SeqVizAnnotation[] = [
    { name: "Fragment_1", start: 0, end: 148, direction: 1 },
    { name: "Fragment_2", start: 148, end: 373, direction: 1 },
    { name: "Fragment_3", start: 373, end: 523, direction: 1 },
    { name: "Fragment_1→Fragment_2 overlap (30bp, designed)", start: 139, end: 169, direction: 1 },
    { name: "Fragment_2→Fragment_3 overlap (28bp, designed)", start: 365, end: 393, direction: 1 },
    { name: "Fragment_3→Fragment_1 overlap (26bp, designed)", start: 505, end: 523, direction: 1 },
    { name: "Fragment_3→Fragment_1 overlap (26bp, designed)", start: 0, end: 8, direction: 1 },
    { name: "Fragment_1 FWD primer (38bp)", start: 505, end: 523, direction: 1 },
    { name: "Fragment_1 FWD primer (38bp)", start: 0, end: 20, direction: 1 },
    { name: "Fragment_1 REV primer (44bp)", start: 125, end: 169, direction: -1 },
    { name: "Fragment_2 FWD primer (34bp)", start: 139, end: 173, direction: 1 },
    { name: "Fragment_2 REV primer (41bp)", start: 352, end: 393, direction: -1 },
    { name: "Fragment_3 FWD primer (34bp)", start: 365, end: 399, direction: 1 },
    { name: "Fragment_3 REV primer (28bp)", start: 503, end: 523, direction: -1 },
    { name: "Fragment_3 REV primer (28bp)", start: 0, end: 8, direction: -1 },
  ];
  const previewAnnNormalized = normalizeAnnotations(previewAnnotationsRaw);
  const previewAnnMerged = mergeOriginCrossingAnnotations(previewAnnNormalized, previewSequence.length);
  const previewAnnotations = previewAnnMerged.map((a) => {
    const name = (a.name ?? "").toLowerCase();
    if (name.includes("overlap") || name.includes("homology")) {
      return { ...a, direction: 0, strand: 0 };
    }
    return a;
  });

  // Show preview only when no active visualization (design, gibson, or pcr)
  const showPreview = !viz && !activePcrViz;

  return (
    <div
      className="w-full relative"
      onDrop={onDrop}
      onDragOver={onDragOver}
      onDragLeave={onDragLeave}
    >
      {/* Drag overlay */}
      {dragActive && (
        <div className="absolute inset-0 z-20 rounded-lg border-2 border-dashed flex items-center justify-center" style={{ borderColor: "#dbefe7", backgroundColor: "rgba(16, 91, 57, 0.95)" }}>
          <div className="text-center">
            <div className="text-lg font-medium" style={{ color: "#dbefe7" }}>Drop GenBank file to attach</div>
            <div className="text-sm" style={{ color: "#dbefe7" }}>.gb / .gbk / .genbank</div>
          </div>
        </div>
      )}

      {/* Centered Chat Section - vertically centered on page */}
      <div className="mx-auto flex flex-col justify-center" style={{ maxWidth: "50rem", minHeight: "35vh", fontFamily: "var(--font-body)" }}>
        {/* Hero: logo + display title + subtitle (mirrors engine-demo) */}
        <div style={{ textAlign: "center", padding: "8px 4px 24px" }}>
          <div style={{ display: "flex", alignItems: "center", justifyContent: "center", gap: 18, marginBottom: 14 }}>
            <img
              src="/splicify-logo.png"
              alt=""
              style={{ width: 64, height: 64, objectFit: "contain", flexShrink: 0 }}
            />
            <h1
              className="splicify-display"
              style={{
                fontSize: 52, fontWeight: 600, letterSpacing: "-0.03em",
                color: "var(--mint-200)", margin: 0, lineHeight: 1.05,
              }}
            >
              Dream a plasmid.
            </h1>
          </div>
          <p
            style={{
              fontSize: 16, lineHeight: 1.55,
              color: "rgba(219,239,231,0.7)",
              margin: "0 auto", maxWidth: 720,
              fontFamily: "var(--font-body)",
            }}
          >
            Annotate modules and the interactions between them. Assemble plasmids from inventory parts with auto-routed cloning workflows, or describe a plasmid and the closest match is picked for you. Verify designs against module-interaction rules, and design CRISPR/Cas9 guides and primers around any feature on the map.
          </p>
        </div>
      </div>

      {/* Chat history (only shown when there are messages) */}
      {hasMessages && (
        <div
          className="mx-auto"
          style={{
            maxWidth: "50rem",
            maxHeight: 380, overflowY: "auto",
            padding: "4px 2px",
            display: "flex", flexDirection: "column", gap: 14,
            marginBottom: 20, fontFamily: "var(--font-body)",
          }}
        >
          {messages.map((m, i) => {
            const isUser = m.role === "user";
            return (
              <div key={i} style={{ display: "flex", justifyContent: isUser ? "flex-end" : "flex-start" }}>
                <div
                  style={{
                    maxWidth: "88%",
                    background: isUser ? "var(--forest-600)" : "rgba(219,239,231,0.06)",
                    border: isUser ? "none" : "1px solid rgba(219,239,231,0.1)",
                    color: "var(--mint-200)",
                    padding: "11px 14px", borderRadius: 14,
                    fontSize: 14, lineHeight: 1.55,
                    whiteSpace: "pre-wrap", wordBreak: "break-word",
                  }}
                >
                  {!isUser && (
                    <div
                      className="splicify-mono"
                      style={{
                        fontSize: 9.5, letterSpacing: "0.18em", textTransform: "uppercase",
                        color: "var(--brass-400)", marginBottom: 6,
                      }}
                    >
                      Splicify
                    </div>
                  )}
                  {formatMessageContent(m.content)}

                  {!isUser && Array.isArray(m.files) && m.files.length > 0 && (
                    <div
                      style={{
                        marginTop: 10,
                        paddingTop: 8,
                        borderTop: "1px solid rgba(219,239,231,0.18)",
                        display: "flex", flexDirection: "column", gap: 4,
                      }}
                    >
                      <div
                        className="splicify-mono"
                        style={{ fontSize: 9.5, letterSpacing: "0.18em", textTransform: "uppercase", color: "rgba(219,239,231,0.55)" }}
                      >
                        Downloads
                      </div>
                      {m.files.map((f, idx) => {
                        const fileName = f.fileName || `file_${idx + 1}`;
                        const url = fileUrlMap[fileName] || null;
                        return (
                          <div key={`${fileName}-${idx}`} style={{ fontSize: 13 }}>
                            {url ? (
                              <a href={url} download={fileName} style={{ color: "var(--brass-400)", textDecoration: "underline" }}>
                                {fileName}
                              </a>
                            ) : (
                              <span style={{ color: "rgba(219,239,231,0.55)" }}>
                                {fileName} <span>(no data)</span>
                              </span>
                            )}
                          </div>
                        );
                      })}
                    </div>
                  )}
                </div>
              </div>
            );
          })}
          <div ref={chatEndRef} />
        </div>
      )}

      {/* Compose section: dark slim input panel + one-click demos (mirrors engine-demo) */}
      <div className="mx-auto" style={{ maxWidth: "50rem", marginTop: "1vh", fontFamily: "var(--font-body)" }}>
      <div
        style={{
          background: "rgba(0,0,0,0.22)",
          border: "1px solid rgba(219,239,231,0.14)",
          borderRadius: 14,
          padding: 6,
          width: "100%",
        }}
      >
        {/* Top: textarea + Design button */}
        <div style={{ display: "flex", gap: 0, alignItems: "stretch" }}>
          <div className="relative" style={{ flex: 1 }}>
            {input.length === 0 && (
              <div
                className="pointer-events-none"
                style={{
                  position: "absolute", left: 14, right: 14, top: 12,
                  color: "rgba(219,239,231,0.55)", fontSize: 15, lineHeight: 1.5,
                  fontFamily: "var(--font-body)",
                }}
              >
                <div className="placeholder-marquee" style={{ overflow: "visible" }}>
                  {prevIdx !== null && (
                    <div className={"placeholder-marquee-line placeholder-out-up"}>{examples[prevIdx]}</div>
                  )}
                  <div
                    className={
                      "placeholder-marquee-line " +
                      (phase === "animating" ? "placeholder-in-from-down" : "placeholder-in")
                    }
                  >
                    {examples[exampleIdx]}
                  </div>
                </div>
              </div>
            )}
            <textarea
              ref={textareaRef}
              className="chat-textarea"
              style={{
                width: "100%", resize: "none", border: "none", outline: "none",
                background: "transparent", color: "var(--mint-200)",
                padding: "12px 14px", fontSize: 15, lineHeight: 1.5,
                fontFamily: "var(--font-body)",
                maxHeight: 160, minHeight: 44,
              }}
              value={input}
              onChange={(e) => setInput(e.target.value)}
              onKeyDown={(e) => {
                if (e.key === "Enter" && !e.shiftKey) {
                  e.preventDefault();
                  send();
                }
              }}
              disabled={loading}
              rows={1}
            />
          </div>
          <button
            data-send-btn
            onClick={send}
            disabled={!canSend}
            style={{
              background: !canSend ? "rgba(216,169,58,0.35)" : "var(--brass-500)",
              color: "var(--forest-900)", border: "none",
              padding: "10px 18px", margin: 4, borderRadius: 10,
              fontFamily: "var(--font-display)", fontSize: 14, fontWeight: 600,
              cursor: !canSend ? "not-allowed" : "pointer", whiteSpace: "nowrap",
            }}
          >
            {loading ? "Designing…" : "Design it →"}
          </button>
        </div>

        {/* Inline tools row: attach + Classic backend toggle */}
        <div
          style={{
            display: "flex", flexWrap: "wrap", alignItems: "center", gap: 8,
            padding: "6px 8px 4px", borderTop: "1px solid rgba(219,239,231,0.08)",
            marginTop: 4,
          }}
        >
          <label
            style={{
              display: "inline-flex", alignItems: "center", gap: 6,
              fontSize: 12, color: "var(--mint-200)",
              cursor: loading ? "not-allowed" : "pointer", opacity: loading ? 0.6 : 1,
            }}
          >
            <span
              style={{
                border: "1px solid rgba(219,239,231,0.2)",
                padding: "5px 11px", borderRadius: 999, fontSize: 12.5,
              }}
            >
              + Inventory .gb
            </span>
            <input
              type="file"
              style={{ display: "none" }}
              accept=".gb,.gbk,.genbank"
              multiple
              onChange={(e) => {
                onPickInventoryFiles(e.target.files);
                e.currentTarget.value = "";
              }}
              disabled={loading}
            />
          </label>

          <label
            style={{
              display: "inline-flex", alignItems: "center", gap: 6,
              fontSize: 12, color: "var(--mint-200)",
              cursor: loading ? "not-allowed" : "pointer", opacity: loading ? 0.6 : 1,
            }}
          >
            <span
              style={{
                border: "1px solid rgba(219,239,231,0.2)",
                padding: "5px 11px", borderRadius: 999, fontSize: 12.5,
              }}
            >
              + Target .gb
            </span>
            <input
              type="file"
              style={{ display: "none" }}
              accept=".gb,.gbk,.genbank"
              onChange={(e) => {
                const f = e.target.files?.[0] ?? null;
                onPickFile(f);
                e.currentTarget.value = "";
              }}
              disabled={loading}
            />
          </label>

          {uploadedInventoryGbs.length > 0 && (
            <button
              onClick={clearInventory}
              disabled={loading}
              style={{
                background: "transparent",
                border: "1px solid rgba(219,239,231,0.18)",
                color: "rgba(219,239,231,0.75)",
                padding: "4px 10px", borderRadius: 999, fontSize: 11.5,
                cursor: loading ? "not-allowed" : "pointer",
              }}
            >
              Clear inventory
            </button>
          )}

          <label
            title="OFF (default): your message goes to the classic backend at /api/chat (deterministic intent dispatcher; fast, no LLM). ON: Agent V2 at /agent_v2/chat (Claude Sonnet 4.6 with AIPlasmidDesign tools + plan.md ReAct loop, four-file output set)."
            style={{
              display: "inline-flex", alignItems: "center", gap: 8,
              marginLeft: "auto", cursor: "pointer", userSelect: "none",
            }}
          >
            <span
              className="splicify-mono"
              style={{
                fontSize: 10.5, letterSpacing: "0.16em", textTransform: "uppercase",
                color: useAgentV2 ? "var(--brass-400)" : "rgba(219,239,231,0.55)",
              }}
            >
              Agent V2
            </span>
            <span
              style={{
                position: "relative", display: "inline-block",
                height: 18, width: 34, borderRadius: 999,
                background: useAgentV2 ? "rgba(216,169,58,0.5)" : "rgba(219,239,231,0.12)",
                border: "1px solid rgba(219,239,231,0.18)",
                transition: "background 160ms ease",
              }}
            >
              <input
                type="checkbox"
                style={{ position: "absolute", opacity: 0, pointerEvents: "none" }}
                checked={useAgentV2}
                onChange={(e) => setUseAgentV2(e.target.checked)}
                disabled={loading}
              />
              <span
                style={{
                  position: "absolute", top: 1, left: 1,
                  height: 14, width: 14, borderRadius: 999,
                  background: useAgentV2 ? "var(--brass-500)" : "rgba(219,239,231,0.85)",
                  transform: useAgentV2 ? "translateX(16px)" : "translateX(0)",
                  transition: "transform 160ms ease, background 160ms ease",
                }}
              />
            </span>
          </label>
        </div>

        {/* File chips */}
        {(uploadedGb || uploadedInventoryGbs.length > 0 || uploadError || inventoryUploadError) && (
          <div style={{ padding: "4px 8px 8px", display: "flex", flexDirection: "column", gap: 6 }}>
            {uploadedGb && (
              <div
                style={{
                  display: "inline-flex", alignItems: "center", gap: 6,
                  border: "1px solid rgba(219,239,231,0.2)", borderRadius: 999,
                  padding: "3px 10px", fontSize: 11.5,
                  color: "var(--mint-200)", width: "fit-content",
                }}
              >
                <span style={{ fontWeight: 600 }}>Target:</span>
                <span>{uploadedGb.name}</span>
                <span style={{ opacity: 0.6 }}>({formatBytes(uploadedGb.size)})</span>
                <button
                  onClick={() => setUploadedGb(null)}
                  title="Remove file"
                  disabled={loading}
                  style={{ background: "transparent", border: "none", color: "var(--mint-200)", cursor: "pointer", padding: 0, marginLeft: 2 }}
                >✕</button>
              </div>
            )}
            {uploadedInventoryGbs.length > 0 && (
              <div style={{ display: "flex", flexDirection: "column", gap: 4 }}>
                <div className="splicify-mono" style={{ fontSize: 9.5, letterSpacing: "0.18em", textTransform: "uppercase", color: "rgba(219,239,231,0.5)" }}>
                  Inventory ({uploadedInventoryGbs.length}/10)
                </div>
                <div style={{ display: "flex", flexWrap: "wrap", gap: 6 }}>
                  {uploadedInventoryGbs.map((f, idx) => (
                    <div
                      key={`${f.name}-${idx}`}
                      style={{
                        display: "inline-flex", alignItems: "center", gap: 6,
                        border: "1px solid rgba(219,239,231,0.2)", borderRadius: 999,
                        padding: "3px 10px", fontSize: 11.5,
                        color: "var(--mint-200)",
                      }}
                    >
                      <span>{f.name}</span>
                      <span style={{ opacity: 0.6 }}>({formatBytes(f.size)})</span>
                      <button
                        onClick={() => removeInventoryFile(idx)}
                        title="Remove file"
                        disabled={loading}
                        style={{ background: "transparent", border: "none", color: "var(--mint-200)", cursor: "pointer", padding: 0, marginLeft: 2 }}
                      >✕</button>
                    </div>
                  ))}
                </div>
              </div>
            )}
            {uploadError && <div style={{ fontSize: 11.5, color: "#fca5a5" }}>{uploadError}</div>}
            {inventoryUploadError && <div style={{ fontSize: 11.5, color: "#fca5a5" }}>{inventoryUploadError}</div>}
          </div>
        )}
      </div>{/* end input panel */}

      {/* One-click demos (pill-style, mirrors engine-demo) */}
      <div style={{ marginTop: 22 }}>
        <div
          className="splicify-mono"
          style={{
            fontSize: 10, letterSpacing: "0.18em", textTransform: "uppercase",
            color: "rgba(219,239,231,0.45)", marginBottom: 12, textAlign: "center",
          }}
        >
          One-click demos
        </div>
        <div style={{ display: "flex", flexWrap: "wrap", gap: 7, justifyContent: "center" }}>
          {([
            ["plasmidAnnotation", "Modular plasmid annotation"],
            ["gibsonFragments", "Gibson fragments"],
            ["sdmWorkflow", "Site-directed mutagenesis"],
            ["restrictionWorkflow", "Restriction cloning"],
            ["sgRNAGoldenGate", "sgRNA Golden Gate"],
            ["goldenGateAssembly", "Golden Gate assembly"],
            ["gatewayCloning", "Gateway cloning"],
          ] as const).map(([key, label]) => (
            <button
              key={key}
              onClick={() => handleExampleClick(key)}
              disabled={loading}
              style={{
                background: "transparent",
                border: "1px solid rgba(219,239,231,0.2)",
                color: "var(--mint-200)",
                padding: "7px 12px", borderRadius: 999,
                fontSize: 12.5, fontFamily: "var(--font-body)",
                cursor: loading ? "not-allowed" : "pointer",
                opacity: loading ? 0.5 : 1,
                transition: "all 160ms ease",
              }}
              onMouseEnter={(e) => { if (!loading) { e.currentTarget.style.borderColor = "var(--brass-400)"; e.currentTarget.style.color = "var(--brass-400)"; } }}
              onMouseLeave={(e) => { e.currentTarget.style.borderColor = "rgba(219,239,231,0.2)"; e.currentTarget.style.color = "var(--mint-200)"; }}
            >
              {label}
            </button>
          ))}
        </div>
      </div>
      </div>{/* End of Compose Section */}

      {/* Loading Animation - shown when sending message */}
      {loading && (
        <div className="mt-12 py-8">
          <SplicifyLoader />
          <div className="text-center mt-4" style={{ color: "#dbefe7" }}>
            Ruminating (can take 30-60 seconds)...
          </div>
        </div>
      )}

      {/* PCR Visualization with tabs - Full Width */}
      {activePcrViz?.type === "pcr" && (
        <div
          className="mt-6 rounded-2xl p-4"
          style={{
            backgroundColor: "#46896c",
            marginLeft: "calc(-50vw + 50%)",
            marginRight: "calc(-50vw + 50%)",
            width: "100vw",
            maxWidth: "100vw",
            paddingLeft: "2rem",
            paddingRight: "2rem",
          }}
        >
          <div className="mb-2 font-medium" style={{ color: "#dbefe7" }}>PCR Visualization</div>

          {pcrTabs.length > 1 && (
            <div className="mb-3 flex w-full gap-2 overflow-x-auto pb-1">
              {pcrTabs.map((v, idx) => {
                const isActive = idx === activePcrIdx;
                const label = pcrTabLabel(v, idx);

                return (
                  <button
                    key={idx}
                    onClick={() => setActivePcrIdx(idx)}
                    className="shrink-0 rounded-full px-4 py-2 text-sm font-medium transition"
                    style={{
                      backgroundColor: isActive ? "#dbefe7" : "#2d4a3e",
                      color: isActive ? "#105b39" : "#dbefe7",
                    }}
                  >
                    {label}
                  </button>
                );
              })}
            </div>
          )}

          <div className="mb-2 text-base" style={{ color: "#dbefe7" }}>
            Template length: {pcrSeq.length} bp
            {typeof pcrAmpliconLen === "number" ? ` • Amplicon: ${pcrAmpliconLen} bp` : ""}
          </div>

          <div className="rounded-2xl p-3 overflow-hidden" style={{ backgroundColor: "#ffffff", height: vizHeight }}>
            <div className="h-full w-full">
              <SeqVizAny
                name="PCR Template"
                seq={pcrSeq}
                sequence={pcrSeq}
                annotations={pcrAnn}
                viewer="linear"
                style={{ height: "100%", width: "100%" }}
              />
            </div>
          </div>
        </div>
      )}

      {/* Plasmid Design Visualization - Full Width */}
      {(designVizList.length > 0 || viz?.type === "design") && (() => {
        const list = designVizList.length > 0
          ? designVizList
          : (viz?.type === "design" ? [viz] : []);
        if (list.length === 0) return null;
        const clampedIdx = Math.max(0, Math.min(activeDesignIdx, list.length - 1));
        const activeViz = list[clampedIdx] as any;
        return (
          <div style={{ marginTop: "calc(2rem - 1vh)" }}>
            {list.length > 1 && (
              <div className="mb-3 flex gap-2 overflow-x-auto pb-1 flex-wrap">
                {list.map((v, idx) => {
                  const isActive = idx === clampedIdx;
                  const tabTitle = (v as any)?.title || `Design ${idx + 1}`;
                  return (
                    <button
                      key={idx}
                      onClick={() => setActiveDesignIdx(idx)}
                      className="shrink-0 rounded-full px-4 py-2 text-sm font-medium transition"
                      style={{
                        backgroundColor: isActive ? "#dbefe7" : "#2d4a3e",
                        color: isActive ? "#105b39" : "#dbefe7",
                      }}
                    >
                      {tabTitle}
                    </button>
                  );
                })}
              </div>
            )}
            <InteractiveSequenceViewer
              key={`design:${clampedIdx}:${(activeViz?.sequence || "").length}:${(activeViz?.sequence || "").slice(0, 40)}`}
              initialSequence={activeViz?.sequence || ""}
              initialAnnotations={normalizeAnnotations(activeViz?.annotations ?? []) as InteractiveAnnotation[]}
              title={activeViz?.title || "Proposed Plasmid"}
              plasmidName={currentPlasmidName || undefined}
              viewerMode="both"
              height={vizHeight}
              circular={activeViz?.topology !== "linear"}
              plannotateEndpoint="/api/plannotate"
              restrictionSites={activeViz?.restriction_sites ?? []}
              analyzeIntentEndpoint="/api/analyze-intent"
              onAnalysisComplete={handleAnalysisComplete}
            />
          </div>
        );
      })()}

      {/* Annotation Visualization */}
      {viz?.type === "annotation" && (
        <div style={{ marginTop: "calc(2rem - 1vh)" }}>
          <InteractiveSequenceViewer
            key={`ann:${(viz as any).title || ""}:${((viz as any).sequence || "").length}:${((viz as any).sequence || "").slice(0, 40)}`}
            initialSequence={(viz as any).sequence || ""}
            initialAnnotations={(viz as any).annotations || []}
            title={(viz as any).title || "Annotated Plasmid"}
            plasmidName={currentPlasmidName || (viz as any).title || undefined}
            viewerMode="both"
            height={vizHeight}
            circular={(viz as any).circular !== false}
            plannotateEndpoint="/api/plannotate"
            analyzeIntentEndpoint="/api/analyze-intent"
            onAnalysisComplete={handleAnalysisComplete}
            autoAnnotateOnMount={(viz as any).auto_annotate === true}
          />
        </div>
      )}

      {/* Gibson Visualization - Full Width */}
      {viz?.type === "gibson" && (
        <div style={{ marginTop: "calc(2rem - 1vh)" }}>
          <InteractiveSequenceViewer
            key={`gibson:${gibsonSeq.length}:${gibsonSeq.slice(0, 40)}`}
            initialSequence={gibsonSeq}
            initialAnnotations={gibsonAnn as InteractiveAnnotation[]}
            title="Assembly Visualization"
            plasmidName={currentPlasmidName || undefined}
            viewerMode="both"
            height={vizHeight}
            circular={true}
            plannotateEndpoint="/api/plannotate"
            analyzeIntentEndpoint="/api/analyze-intent"
            onAnalysisComplete={handleAnalysisComplete}
          />
        </div>
      )}

      {/* Interactive Splicify Preview - Full Width (shown when no active visualization and not loading) */}
      {showPreview && !loading && (
        <div style={{ marginTop: "calc(2rem - 1vh)" }}>
          <InteractiveSequenceViewer
            key={`preview:${previewSequence.length}:${previewSequence.slice(0, 40)}`}
            initialSequence={previewSequence}
            initialAnnotations={previewAnnotations as InteractiveAnnotation[]}
            title="Plasmid Visualizer"
            plasmidName={currentPlasmidName || undefined}
            viewerMode="both"
            height={vizHeight}
            plannotateEndpoint="/api/plannotate"
            analyzeIntentEndpoint="/api/analyze-intent"
            onAnalysisComplete={handleAnalysisComplete}
          />
        </div>
      )}
    </div>
  );
}
