"use client";

import { useEffect, useMemo, useRef, useState } from "react";
import dynamic from "next/dynamic";
import { writeZip, ZipEntry } from "../lib/zip";

// Heavy viewer — dynamic-import so the chat page stays light until a
// plasmid is uploaded or a .gb output file lands.
const InteractiveSequenceViewer = dynamic(
  () => import("../components/InteractiveSequenceViewer"),
  { ssr: false, loading: () => (
      <div style={{ padding: 12, fontSize: 12, opacity: 0.7 }}>Loading viewer…</div>
    ) },
);


type FileEnvelope = { fileName: string; dataBase64: string };

type Viz = {
  title?: string;
  sequence: string;
  circular?: boolean;
  annotations?: any[];
  hierarchical_annotations?: any[];
  interactions?: any[];
  cloning_features?: any;
};

type Citation = { plasmid_id?: string; via?: string; ref?: string };

type ToolCall = {
  iteration: number;
  tool: string;
  input?: Record<string, any>;
  status: "running" | "done";
  startedAt: number;
  elapsed_ms?: number;
  n_results?: number | null;
  ok?: boolean | null;
};

type Turn = {
  id: string;
  question: string;
  fileNames: string[];
  // Streaming/pending state
  status: "pending" | "done" | "error";
  shorthand?: string;
  intent?: string;
  toolCalls?: ToolCall[];
  // Final assistant state
  reply?: string;
  outputFileNames?: string[];
  outputFiles?: FileEnvelope[];
  nToolCalls?: number;
  error?: string;
  citations?: Citation[];
};


function briefArgs(input: Record<string, any> | undefined): string {
  if (!input) return "";
  const pairs = Object.entries(input)
    .filter(([_, v]) => v !== undefined && v !== null && v !== "")
    .map(([k, v]) => {
      const s = typeof v === "string" ? v : JSON.stringify(v);
      return `${k}=${s.length > 40 ? s.slice(0, 38) + "…" : s}`;
    });
  return pairs.join(" · ");
}


const PLASMID_EXTENSIONS = new Set([".gb", ".gbk", ".genbank", ".dna"]);
const ACCEPT_TYPES =
  ".gb,.gbk,.genbank,.dna,.fa,.fasta,.fna,.fastq,.fq,.ab1,.csv,.tsv,.txt,.json";

function fileExt(name: string): string {
  const ix = name.lastIndexOf(".");
  return ix < 0 ? "" : name.slice(ix).toLowerCase();
}
function isPlasmidFile(name: string): boolean {
  return PLASMID_EXTENSIONS.has(fileExt(name));
}

function base64ToBytes(b64: string): Uint8Array<ArrayBuffer> {
  const bin = atob(b64);
  const bytes = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) bytes[i] = bin.charCodeAt(i);
  return bytes;
}

function base64ToBlob(b64: string, mime = "application/octet-stream"): Blob {
  return new Blob([base64ToBytes(b64)], { type: mime });
}

function downloadAllAsZip(
  files: FileEnvelope[],
  reply: string | undefined,
  zipName: string,
) {
  const entries: ZipEntry[] = files.map((f) => ({
    name: f.fileName,
    bytes: base64ToBytes(f.dataBase64),
  }));
  if (reply && reply.trim()) {
    entries.push({ name: "reply.txt", bytes: new TextEncoder().encode(reply) });
  }
  const blob = writeZip(entries);
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = zipName;
  document.body.appendChild(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}


// Map the hierarchical annotation payload from /agent_v2/annotate-on-upload
// into the viewer's annotation shape.
function combineVizAnnotations(viz: Viz): any[] {
  // Forward the server's color when set; otherwise leave the field
  // undefined so the viewer's featureColor logic can pick a hue from
  // the design / preserved-gray / feature-alt palettes. The previous
  // hardcoded #7C3AED / #6b5b95 fallbacks made every uncolored row
  // paint as the same purple — which masked added_by_design entirely.
  // added_by_design is preserved so the design palette can light up
  // pegRNAs / ngRNAs / primers on agent-emitted .gb files.
  const base = (viz.annotations || []).map((a: any) => ({
    name: a.name || "annotation",
    start: a.start, end: a.end,
    direction: a.direction === -1 ? -1 : a.direction === 0 ? 0 : 1,
    strand: a.direction === -1 ? -1 : a.direction === 0 ? 0 : 1,
    color: a.color || undefined,
    type: a.type || "misc_feature",
    description: a.description || "",
    layer: a.layer || "feature",
    added_by_design: a.added_by_design || undefined,
    sseqid: a.sseqid, db: a.db, kb_data: a.kb_data ?? null,
  }));
  const hier = (viz.hierarchical_annotations || []).map((a: any) => ({
    name: a.name || "annotation",
    start: a.start, end: a.end,
    direction: a.direction === -1 ? -1 : a.direction === 0 ? 0 : 1,
    strand: a.direction === -1 ? -1 : a.direction === 0 ? 0 : 1,
    color: a.color || undefined,
    type: a.type || a.motif_type || "misc_feature",
    description: a.description || "",
    layer: a.layer || "feature",
    added_by_design: a.added_by_design || undefined,
    motif_type: a.motif_type, module_type: a.module_type,
    source: a.source, payload_id: a.payload_id, metadata: a.metadata,
  }));
  return [...base, ...hier];
}

// Parse the interpreter's citation refs (e.g. "Cas9 5109-9214") into a
// numeric coordinate range. Returns the union of every numeric range
// found in the citations as {start, end}, or null when none parse.
// The viewer uses this to recenter on the region of interest.
function extractFocusRegion(citations: Citation[] | undefined): { start: number; end: number } | null {
  if (!citations || citations.length === 0) return null;
  const ranges: [number, number][] = [];
  for (const c of citations) {
    if (!c?.ref) continue;
    const matches = c.ref.matchAll(/(\d+)\s*[-–]\s*(\d+)/g);
    for (const m of matches) {
      const a = parseInt(m[1], 10);
      const b = parseInt(m[2], 10);
      if (!Number.isFinite(a) || !Number.isFinite(b)) continue;
      ranges.push([Math.min(a, b), Math.max(a, b)]);
    }
  }
  if (ranges.length === 0) return null;
  return {
    start: Math.min(...ranges.map((r) => r[0])),
    end: Math.max(...ranges.map((r) => r[1])),
  };
}


async function annotatePlasmidFile(file: File): Promise<Viz | null> {
  try {
    const form = new FormData();
    form.append("file", file, file.name);
    const r = await fetch("/api/agent_v2/annotate-on-upload", { method: "POST", body: form });
    if (!r.ok) return null;
    const data = await r.json();
    if (!data?.ok || !data?.viz?.sequence) return null;
    return data.viz as Viz;
  } catch { return null; }
}

async function annotatePlasmidEnvelope(env: FileEnvelope): Promise<Viz | null> {
  try {
    const bin = atob(env.dataBase64);
    const bytes = new Uint8Array(bin.length);
    for (let i = 0; i < bin.length; i++) bytes[i] = bin.charCodeAt(i);
    const file = new File([bytes], env.fileName, { type: "chemical/seq-na-genbank" });
    return await annotatePlasmidFile(file);
  } catch { return null; }
}


export default function SplicifyPage() {
  // Active draft (what the user is typing right now)
  const [message, setMessage] = useState("");
  const [files, setFiles] = useState<File[]>([]);
  const [sending, setSending] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [sessionId, setSessionId] = useState<string>("");

  // Chat history — every turn is a user→assistant exchange.
  const [history, setHistory] = useState<Turn[]>([]);

  // Per-file viz hydrations. Inputs accumulate across the session;
  // outputs are keyed per-turn so an earlier turn's assembled.gb keeps
  // rendering even after the next turn lands.
  const [inputViz, setInputViz] = useState<Record<string, Viz>>({});
  const [outputViz, setOutputViz] = useState<Record<string, Viz>>({});
  const [annotatingNames, setAnnotatingNames] = useState<Set<string>>(new Set());

  // Blob URLs for output-file download chips.
  const urlMapRef = useRef<Record<string, string>>({});
  useEffect(() => () => {
    Object.values(urlMapRef.current).forEach((u) => URL.revokeObjectURL(u));
  }, []);
  function urlFor(file: FileEnvelope): string {
    const cached = urlMapRef.current[file.fileName];
    if (cached) return cached;
    const url = URL.createObjectURL(base64ToBlob(file.dataBase64));
    urlMapRef.current[file.fileName] = url;
    return url;
  }

  // Auto-scroll the chat history to the bottom whenever a new turn lands.
  const historyEndRef = useRef<HTMLDivElement | null>(null);
  useEffect(() => {
    historyEndRef.current?.scrollIntoView({ behavior: "smooth", block: "end" });
  }, [history]);

  // The viewer should focus on the most recent answer's coordinates so
  // the user sees the region the interpreter cited. Derive from the
  // last `done` turn that has citations.
  const focusRegion = useMemo(() => {
    for (let i = history.length - 1; i >= 0; i--) {
      const t = history[i];
      if (t.status === "done" && t.citations?.length) {
        return extractFocusRegion(t.citations);
      }
    }
    return null;
  }, [history]);

  function hydrateInputFile(file: File) {
    if (!isPlasmidFile(file.name)) return;
    setAnnotatingNames((prev) => new Set(prev).add(file.name));
    annotatePlasmidFile(file).then((viz) => {
      setAnnotatingNames((prev) => {
        const next = new Set(prev);
        next.delete(file.name);
        return next;
      });
      if (viz) setInputViz((prev) => ({ ...prev, [file.name]: viz }));
    });
  }

  // Hydrate plasmid output files for the embedded viewer as soon as they
  // land — including mid-stream files emitted by find_genomic_record (the
  // retrieved SIRT6.gb / CGAS.gb) and emit_guides_gb (the annotated output)
  // BEFORE the turn finishes. Previously we gated this on status==="done"
  // so the user only saw the viewer light up at the end of the run; now
  // the viewer surfaces the retrieved sequence the moment NCBI returns it
  // and replaces it with the annotated version when emit_guides_gb lands.
  useEffect(() => {
    let cancelled = false;
    (async () => {
      for (const turn of history) {
        if (turn.status === "error" || !turn.outputFiles) continue;
        for (const f of turn.outputFiles) {
          if (!isPlasmidFile(f.fileName) || outputViz[f.fileName]) continue;
          const viz = await annotatePlasmidEnvelope(f);
          if (cancelled) return;
          if (viz) setOutputViz((prev) => ({ ...prev, [f.fileName]: viz }));
        }
      }
    })();
    return () => { cancelled = true; };
  }, [history, outputViz]);

  async function send() {
    const text = message.trim();
    if (!text || sending) return;

    const turnId = `turn_${Date.now()}_${Math.random().toString(36).slice(2, 8)}`;
    const turn: Turn = {
      id: turnId,
      question: text,
      fileNames: files.map((f) => f.name),
      status: "pending",
    };
    setHistory((prev) => [...prev, turn]);
    setSending(true);
    setError(null);

    // Snapshot files for the request; clear the staging input so the
    // user can immediately type the next question while we wait.
    const sentFiles = files.slice();
    setMessage("");
    setFiles([]);

    const form = new FormData();
    form.append("message", text);
    if (sessionId) form.append("session_id", sessionId);
    for (const f of sentFiles) form.append("inventory_files", f, f.name);

    try {
      const res = await fetch("/api/agent_v2/chat-stream", { method: "POST", body: form });
      if (!res.body) throw new Error("No response body from server");

      const reader = res.body.getReader();
      const decoder = new TextDecoder();
      let buf = "";

      while (true) {
        const { done, value } = await reader.read();
        if (done) break;
        buf += decoder.decode(value, { stream: true });
        const events = buf.split("\n\n");
        buf = events.pop() || "";

        for (const ev of events) {
          if (!ev.trim()) continue;
          let eventName = "message";
          let dataStr = "";
          for (const line of ev.split("\n")) {
            if (line.startsWith("event: ")) eventName = line.slice(7).trim();
            else if (line.startsWith("data: ")) dataStr += line.slice(6);
          }
          if (!dataStr) continue;
          let payload: any;
          try { payload = JSON.parse(dataStr); } catch { continue; }

          if (eventName === "shorthand") {
            setHistory((prev) => prev.map((t) =>
              t.id === turnId ? { ...t, shorthand: payload.shorthand || "", intent: payload.intent } : t
            ));
          } else if (eventName === "tool_call") {
            setHistory((prev) => prev.map((t) => {
              if (t.id !== turnId) return t;
              const calls = [...(t.toolCalls || [])];
              const idx = calls.findIndex((c) => c.iteration === payload.iteration);
              if (payload.phase === "start") {
                const newCall: ToolCall = {
                  iteration: payload.iteration,
                  tool: payload.tool,
                  input: payload.input,
                  status: "running",
                  startedAt: Date.now(),
                };
                if (idx >= 0) calls[idx] = { ...calls[idx], ...newCall };
                else calls.push(newCall);
              } else if (payload.phase === "end" && idx >= 0) {
                calls[idx] = {
                  ...calls[idx],
                  status: "done",
                  elapsed_ms: payload.elapsed_ms,
                  n_results: payload.n_results ?? null,
                  ok: payload.ok ?? null,
                };
              }
              // Mid-stream file envelope (e.g. retrieved genomic .gb from
              // find_genomic_record): append to outputFiles so the active
              // viewer renders the .gb as soon as the tool returns,
              // instead of waiting for the full pipeline to finish.
              let outputFiles = t.outputFiles || [];
              if (payload.phase === "end"
                  && payload.file
                  && typeof payload.file === "object"
                  && typeof payload.file.fileName === "string"
                  && typeof payload.file.dataBase64 === "string") {
                const incoming = payload.file as { fileName: string; dataBase64: string };
                const exists = outputFiles.some((f: any) => f.fileName === incoming.fileName);
                if (!exists) outputFiles = [...outputFiles, incoming];
              }
              return { ...t, toolCalls: calls, outputFiles };
            }));
          } else if (eventName === "envelope") {
            if (payload.session_id) setSessionId(payload.session_id);
            setHistory((prev) => prev.map((t) =>
              t.id === turnId ? {
                ...t,
                status: payload.ok === false && payload.error ? "error" : "done",
                reply: payload.reply || "",
                outputFiles: Array.isArray(payload.files) ? payload.files : [],
                nToolCalls: typeof payload.n_tool_calls === "number" ? payload.n_tool_calls : undefined,
                error: payload.ok === false ? payload.error : undefined,
                citations: payload.citations || [],
              } : t
            ));
          } else if (eventName === "error") {
            setHistory((prev) => prev.map((t) =>
              t.id === turnId ? { ...t, status: "error", error: payload.error || "Unknown error" } : t
            ));
          }
        }
      }
    } catch (e) {
      const msg = e instanceof Error ? e.message : String(e);
      setError(msg);
      setHistory((prev) => prev.map((t) =>
        t.id === turnId ? { ...t, status: "error", error: msg } : t
      ));
    } finally {
      setSending(false);
    }
  }

  function newTopic() {
    setSessionId("");
    setMessage("");
    setFiles([]);
    setHistory([]);
    setInputViz({});
    setOutputViz({});
    setError(null);
    Object.values(urlMapRef.current).forEach((u) => URL.revokeObjectURL(u));
    urlMapRef.current = {};
  }

  const hasUploadedPlasmid = Object.keys(inputViz).length > 0;
  const inputViewers = Object.entries(inputViz);

  // Single active viewer surface. Output plasmids from the latest turn
  // win over uploaded inputs so completing a workflow (assembled.gb,
  // guides.gb) replaces the original view in place — Claude-style:
  // one persistent surface that updates as the conversation advances,
  // not a new viewer per turn. The viewer's React key includes
  // activeName so it remounts cleanly when the active plasmid changes.
  let activeName: string | null = null;
  let activeViz: Viz | null = null;
  for (let i = history.length - 1; i >= 0 && activeName === null; i--) {
    const turn = history[i];
    for (let j = (turn.outputFiles || []).length - 1; j >= 0; j--) {
      const f = turn.outputFiles![j];
      if (isPlasmidFile(f.fileName) && outputViz[f.fileName]) {
        activeName = f.fileName;
        activeViz = outputViz[f.fileName];
        break;
      }
    }
  }
  if (activeName === null && inputViewers.length > 0) {
    const [name, viz] = inputViewers[inputViewers.length - 1];
    activeName = name;
    activeViz = viz;
  }

  return (
    <main
      style={{
        minHeight: "100vh",
        background: "linear-gradient(180deg, var(--forest-900) 0%, var(--forest-700) 100%)",
        color: "var(--mint-200)",
        fontFamily: "var(--font-body)",
        padding: "48px 24px 64px",
      }}
    >
      <div style={{ maxWidth: 1080, margin: "0 auto" }}>
        <header style={{ marginBottom: 28 }}>
          <div
            style={{
              fontFamily: "var(--font-display)", fontSize: 12,
              letterSpacing: "0.24em", textTransform: "uppercase",
              color: "var(--brass-400)", marginBottom: 8,
            }}
          >
            Splicify · AI Agent
          </div>
          <h1 style={{ fontFamily: "var(--font-display)", fontSize: 38, margin: 0 }}>
            AI Molecular Biologist
          </h1>
          {history.length === 0 && (
            <p style={{ fontSize: 15, opacity: 0.85, marginTop: 12, maxWidth: 640 }}>
              Describe a plasmid you want to build — or ask a question about an uploaded one
              (<em>"What are the coordinates of Cas9?"</em>, <em>"What promoter drives Cas9?"</em>).
              Drop one or more files (<code>.gb</code>, <code>.fasta</code>, <code>.csv</code>,{" "}
              <code>.json</code>, primers, parts lists…) and the agent returns
              <code> assembled.gb</code>, <code>parts_order.csv</code>,
              <code> protocol.csv</code>, and a workflow trace. Uploaded and assembled plasmids
              render inline so you can inspect features, modules, interaction chords, and
              translation strips.
            </p>
          )}
        </header>

        {/* Chat history — one block per user-assistant turn. Mirrors
            aiplasmiddesign Chat.tsx: a scrolling inner panel (maxHeight
            + overflowY: auto) so the chat doesn't push the rest of the
            page out of view. Flows with the page; no position: sticky. */}
        {history.length > 0 && (
          <section
            style={{
              marginBottom: 16,
              maxHeight: "60vh",
              overflowY: "auto",
              paddingRight: 4,
            }}
          >
            {history.map((turn) => (
              <article key={turn.id} style={{ marginBottom: 18 }}>
                {/* User bubble */}
                <div
                  style={{
                    background: "rgba(8,56,35,0.55)",
                    border: "1px solid rgba(219,239,231,0.1)",
                    borderRadius: 12, padding: "10px 14px",
                    fontSize: 14.5, marginBottom: 8,
                    whiteSpace: "pre-wrap",
                  }}
                >
                  <div style={{ fontSize: 11, opacity: 0.6, marginBottom: 4,
                                 letterSpacing: "0.08em", textTransform: "uppercase" }}>
                    You
                  </div>
                  {turn.question}
                  {turn.fileNames.length > 0 && (
                    <div style={{ marginTop: 8, display: "flex", gap: 6, flexWrap: "wrap" }}>
                      {turn.fileNames.map((n) => (
                        <span
                          key={n}
                          style={{
                            fontSize: 11, padding: "2px 8px",
                            border: "1px solid rgba(219,239,231,0.2)",
                            borderRadius: 999,
                          }}
                        >
                          {n}
                        </span>
                      ))}
                    </div>
                  )}
                </div>

                {/* Shorthand / working indicator */}
                {turn.shorthand && turn.status !== "done" && (
                  <div
                    style={{
                      fontSize: 13, color: "rgba(219,239,231,0.75)", fontStyle: "italic",
                      padding: "8px 14px", borderLeft: "2px solid var(--brass-400)",
                      marginBottom: 8, background: "rgba(216,169,58,0.05)",
                    }}
                  >
                    <strong style={{ color: "var(--brass-400)", fontStyle: "normal" }}>
                      {turn.intent === "REJECT" ? "Out of scope: "
                        : turn.intent === "CRISPR_GUIDE" ? "CRISPR: "
                        : "Working on: "}
                    </strong>
                    {turn.shorthand}
                    {turn.status === "pending" && <span style={{ marginLeft: 6 }}>…</span>}
                  </div>
                )}
                {turn.status === "pending" && !turn.shorthand && (
                  <div
                    style={{
                      fontSize: 13, color: "rgba(219,239,231,0.6)", fontStyle: "italic",
                      padding: "8px 14px", borderLeft: "2px solid rgba(219,239,231,0.2)",
                      marginBottom: 8,
                    }}
                  >
                    Thinking…
                  </div>
                )}

                {/* Live tool-call strip — one row per dispatched tool with
                    a live ⏳/✓ indicator and elapsed-ms timing once the
                    call completes. */}
                {turn.toolCalls && turn.toolCalls.length > 0 && (
                  <div
                    style={{
                      border: "1px solid rgba(219,239,231,0.08)",
                      borderRadius: 10,
                      background: "rgba(8,56,35,0.35)",
                      padding: "6px 10px",
                      marginBottom: 8,
                      fontFamily: "var(--font-mono, monospace)",
                      fontSize: 11.5,
                    }}
                  >
                    <div style={{
                      fontSize: 10, opacity: 0.55, marginBottom: 4,
                      letterSpacing: "0.08em", textTransform: "uppercase",
                      fontFamily: "var(--font-body)",
                    }}>
                      Tool calls
                    </div>
                    {turn.toolCalls.map((c) => (
                      <div
                        key={`tc-${c.iteration}`}
                        style={{
                          display: "grid",
                          gridTemplateColumns: "16px minmax(140px, auto) 1fr minmax(60px, auto) minmax(40px, auto)",
                          gap: 8, alignItems: "center", padding: "2px 0",
                          opacity: c.status === "running" ? 0.85 : 1,
                        }}
                      >
                        <span style={{
                          color: c.status === "running" ? "var(--brass-400)"
                                  : c.ok === false ? "#ff7d7d"
                                  : "var(--mint-200)",
                        }}>
                          {c.status === "running" ? "⏳" : c.ok === false ? "✕" : "✓"}
                        </span>
                        <span style={{ color: "var(--brass-400)" }}>{
                          c.tool === "_phase_explore" ? "Exploring (3 subagents in parallel)" :
                          c.tool === "_phase_plan" ? "Planning" :
                          c.tool === "_phase_main" ? "Designing (main agent)" :
                          c.tool === "_phase_summarize" ? "Summarizing" :
                          c.tool
                        }</span>
                        <span style={{ opacity: 0.6, overflow: "hidden",
                                        textOverflow: "ellipsis", whiteSpace: "nowrap" }}>
                          {briefArgs(c.input)}
                        </span>
                        <span style={{ textAlign: "right", opacity: 0.75 }}>
                          {c.elapsed_ms != null
                            ? `${c.elapsed_ms < 1000
                                ? c.elapsed_ms + " ms"
                                : (c.elapsed_ms / 1000).toFixed(2) + " s"}`
                            : "…"}
                        </span>
                        <span style={{ textAlign: "right", opacity: 0.55 }}>
                          {c.n_results != null ? `${c.n_results} hits` : ""}
                        </span>
                      </div>
                    ))}
                  </div>
                )}

                {/* Assistant reply */}
                {turn.status === "done" && turn.reply && (
                  <div
                    style={{
                      background: "rgba(8,56,35,0.5)",
                      border: "1px solid rgba(219,239,231,0.1)",
                      borderRadius: 12, padding: "12px 16px",
                      whiteSpace: "pre-wrap", fontSize: 14.5, lineHeight: 1.6,
                    }}
                  >
                    <div style={{ fontSize: 11, opacity: 0.6, marginBottom: 6,
                                   letterSpacing: "0.08em", textTransform: "uppercase" }}>
                      Splicify
                    </div>
                    {turn.reply}
                    {typeof turn.nToolCalls === "number" && (
                      <div style={{ marginTop: 10, fontSize: 11, opacity: 0.5, fontStyle: "italic" }}>
                        {turn.nToolCalls} tool call{turn.nToolCalls === 1 ? "" : "s"}
                      </div>
                    )}
                  </div>
                )}

                {turn.status === "error" && (
                  <div
                    style={{
                      padding: "10px 14px", borderRadius: 10,
                      background: "rgba(180,40,40,0.15)",
                      border: "1px solid rgba(180,40,40,0.5)",
                      color: "#ffd7d7", fontSize: 13,
                    }}
                  >
                    <strong>Error:</strong> {turn.error || "unknown"}
                  </div>
                )}

                {/* Output files. Only show ONE plasmid (.gb-class) chip
                    per turn — the LATEST in outputFiles, which is the
                    fully annotated emit_guides_gb / emit_assembled_gb
                    output. Earlier mid-stream .gb files (e.g. the raw
                    retrieved SIRT6.gb from find_genomic_record) feed the
                    viewer but are NOT offered as separate downloads,
                    matching the user's 'provide one annotated .gb file
                    to download' contract. Non-plasmid files (CSV / TXT)
                    are listed individually. */}
                {turn.outputFiles && turn.outputFiles.length > 0 && (() => {
                  const lastPlasmidIdx = (() => {
                    for (let i = turn.outputFiles!.length - 1; i >= 0; i--) {
                      if (isPlasmidFile(turn.outputFiles![i].fileName)) return i;
                    }
                    return -1;
                  })();
                  const chips = turn.outputFiles!.filter((f, i) => (
                    !isPlasmidFile(f.fileName) || i === lastPlasmidIdx
                  ));
                  const showDownloadAll = turn.status === "done" && (chips.length > 0 || (turn.reply && turn.reply.trim().length > 0));
                  const zipName = `splicify_${turn.id}.zip`;
                  return (
                    <div style={{ display: "flex", flexWrap: "wrap", gap: 8, marginTop: 10 }}>
                      {chips.map((f) => (
                        <a
                          key={f.fileName}
                          href={urlFor(f)}
                          download={f.fileName}
                          style={{
                            display: "inline-flex", alignItems: "center", gap: 8,
                            border: "1px solid var(--brass-500)",
                            background: "rgba(216,169,58,0.08)",
                            color: "var(--brass-400)",
                            borderRadius: 10, padding: "6px 12px",
                            fontSize: 12.5, fontFamily: "var(--font-mono, monospace)",
                            textDecoration: "none",
                          }}
                        >
                          ⬇ {f.fileName}
                        </a>
                      ))}
                      {showDownloadAll && (
                        <button
                          type="button"
                          onClick={() => downloadAllAsZip(chips, turn.reply, zipName)}
                          title="Download all files + reply as a .zip"
                          style={{
                            display: "inline-flex", alignItems: "center", gap: 8,
                            border: "1px solid var(--brass-500)",
                            background: "rgba(216,169,58,0.18)",
                            color: "var(--brass-400)",
                            borderRadius: 10, padding: "6px 12px",
                            fontSize: 12.5, fontFamily: "var(--font-mono, monospace)",
                            cursor: "pointer",
                            fontWeight: 600,
                          }}
                        >
                          ⬇ Download All
                        </button>
                      )}
                    </div>
                  );
                })()}

                {/* Per-turn output plasmid viewer removed — the single
                    activeViz viewer at the top of the page surfaces the
                    output instead so the user sees one continuous view
                    that updates as the workflow lands. Pending-annotation
                    hint stays below the download chips. */}
                {turn.outputFiles?.filter((f) => isPlasmidFile(f.fileName) && !outputViz[f.fileName]).map((f) => (
                  <div key={`turn-${turn.id}-out-${f.fileName}-pending`}
                       style={{ fontSize: 12, opacity: 0.7, padding: "10px 14px",
                                border: "1px solid rgba(216,169,58,0.25)",
                                borderRadius: 10, marginTop: 10 }}>
                    Annotating {f.fileName}…
                  </div>
                ))}
              </article>
            ))}
            <div ref={historyEndRef} />
          </section>
        )}

        {/* Pending-input file chips (only visible before the user clicks Send). */}
        {files.some((f) => isPlasmidFile(f.name) && annotatingNames.has(f.name)) && (
          <div style={{ fontSize: 12, opacity: 0.7, marginBottom: 8 }}>
            Annotating {Array.from(annotatingNames).join(", ")}…
          </div>
        )}

        {/* Composer — static, sits directly below the chat history.
            Mirrors aiplasmiddesign Chat.tsx (no sticky / fixed position;
            the page scrolls naturally). The chat history above this
            block has its own maxHeight + overflowY: auto so the chat
            panel itself doesn't grow unboundedly. */}
        <section
          style={{
            background: "rgba(11,74,48,0.55)",
            border: "1px solid rgba(219,239,231,0.12)",
            borderRadius: 14,
            padding: 16,
            marginTop: 8,
          }}
        >
          <textarea
            value={message}
            onChange={(e) => setMessage(e.target.value)}
            placeholder={
              hasUploadedPlasmid
                ? "Ask about the uploaded plasmid, or describe a build…"
                : "e.g. Build me a Gibson assembly with hPGK-GFP into pUC19. Or: What is the application of this plasmid?"
            }
            disabled={sending}
            onKeyDown={(e) => {
              if (e.key === "Enter" && (e.metaKey || e.ctrlKey)) {
                e.preventDefault();
                send();
              }
            }}
            style={{
              width: "100%", minHeight: 88, padding: 12, borderRadius: 10,
              border: "1px solid rgba(219,239,231,0.18)",
              background: "rgba(8,56,35,0.55)",
              color: "var(--mint-200)",
              fontSize: 15, lineHeight: 1.5,
              fontFamily: "var(--font-body)",
              resize: "vertical",
            }}
          />

          <div
            style={{
              display: "flex", alignItems: "center", gap: 10,
              marginTop: 12, flexWrap: "wrap",
            }}
          >
            <label
              style={{
                display: "inline-flex", alignItems: "center", gap: 6,
                fontSize: 12, color: "var(--mint-200)",
                cursor: sending ? "not-allowed" : "pointer",
                opacity: sending ? 0.6 : 1,
              }}
            >
              <span
                style={{
                  border: "1px solid rgba(219,239,231,0.2)",
                  padding: "5px 11px", borderRadius: 999, fontSize: 12.5,
                }}
              >
                + files
              </span>
              <input
                type="file" multiple accept={ACCEPT_TYPES}
                style={{ display: "none" }}
                onChange={(e) => {
                  const fs = Array.from(e.target.files || []);
                  setFiles((prev) => [...prev, ...fs]);
                  for (const f of fs) hydrateInputFile(f);
                  e.currentTarget.value = "";
                }}
                disabled={sending}
              />
            </label>

            {files.map((f, i) => (
              <span
                key={`${f.name}-${i}`}
                style={{
                  display: "inline-flex", alignItems: "center", gap: 6,
                  border: "1px solid rgba(219,239,231,0.2)",
                  borderRadius: 999, padding: "3px 10px", fontSize: 11.5,
                }}
              >
                {f.name}
                <button
                  onClick={() => {
                    setFiles((prev) => prev.filter((_, idx) => idx !== i));
                    setInputViz((prev) => {
                      const next = { ...prev };
                      delete next[f.name];
                      return next;
                    });
                  }}
                  disabled={sending}
                  style={{
                    background: "transparent", border: "none",
                    color: "var(--mint-200)", cursor: "pointer", padding: 0,
                  }}
                  aria-label={`remove ${f.name}`}
                >
                  ✕
                </button>
              </span>
            ))}

            <div style={{ flex: 1 }} />

            {(history.length > 0 || sessionId) && (
              <button
                onClick={newTopic}
                disabled={sending}
                style={{
                  background: "transparent",
                  border: "1px solid rgba(219,239,231,0.18)",
                  color: "var(--mint-200)",
                  padding: "6px 12px", borderRadius: 8, fontSize: 12,
                  cursor: sending ? "not-allowed" : "pointer",
                }}
                title="Clear the conversation and start fresh"
              >
                New topic
              </button>
            )}

            <button
              onClick={send}
              disabled={!message.trim() || sending}
              style={{
                background: !message.trim() || sending ? "rgba(216,169,58,0.35)" : "var(--brass-500)",
                color: "var(--forest-900)", border: "none",
                padding: "8px 18px", borderRadius: 10,
                fontFamily: "var(--font-display)", fontSize: 14, fontWeight: 600,
                cursor: !message.trim() || sending ? "not-allowed" : "pointer",
              }}
            >
              {sending ? "Sending…" : "Send"}
            </button>
          </div>
        </section>

        {error && (
          <div
            style={{
              padding: "10px 14px", borderRadius: 10,
              background: "rgba(180,40,40,0.15)",
              border: "1px solid rgba(180,40,40,0.5)",
              color: "#ffd7d7", fontSize: 13, marginTop: 12,
            }}
          >
            <strong>Error:</strong> {error}
          </div>
        )}

        {/* Single active viewer — present from page load so the user has
            a visual anchor before any plasmid is uploaded. Shows an empty
            placeholder until either (a) the user uploads a .gb, or (b) a
            mid-stream tool result (find_genomic_record / emit_guides_gb /
            emit_assembled_gb) lands a plasmid envelope. The viewer surface
            updates in place — never a separate per-turn viewer — so one
            continuous view tracks the workflow's progress. */}
        <section style={{ marginBottom: 18 }}>
          <div
            style={{
              border: "1px solid rgba(219,239,231,0.12)",
              borderRadius: 12,
              background: "rgba(8,56,35,0.45)",
              padding: 8,
              minHeight: 200,
            }}
          >
            {activeViz ? (
              <>
                <div style={{ fontSize: 12, opacity: 0.75, marginBottom: 6, paddingLeft: 4 }}>
                  {/* Label echoes whether this is an upload or an output. */}
                  {Object.prototype.hasOwnProperty.call(inputViz, activeName!)
                    ? `Uploaded · ${activeName} · ${(activeViz.sequence || "").length} bp`
                    : `Output · ${activeName} · ${(activeViz.sequence || "").length} bp`}
                </div>
                <InteractiveSequenceViewer
                  key={`active-${activeName}`}
                  initialSequence={activeViz.sequence || ""}
                  initialAnnotations={combineVizAnnotations(activeViz) as any}
                  initialInteractions={(activeViz.interactions as any) || []}
                  initialCloningFeatures={(activeViz.cloning_features as any) || null}
                  focusRegion={focusRegion}
                  title={activeName || undefined}
                  plasmidName={(activeName || "").replace(/\.[^.]+$/, "")}
                  viewerMode="both"
                  height={520}
                  circular={activeViz.circular !== false}
                  initialLeftPaneView={
                    (activeViz as any).type === "genomic" || activeViz.circular === false
                      ? "flat"
                      : "circular"
                  }
                  plannotateEndpoint="/api/plannotate"
                  autoAnnotateOnMount={false}
                />
              </>
            ) : (
              <div
                style={{
                  display: "flex", flexDirection: "column",
                  alignItems: "center", justifyContent: "center",
                  minHeight: 200, padding: 24,
                  textAlign: "center", opacity: 0.65,
                  fontSize: 13.5, lineHeight: 1.55,
                  color: "var(--mint-200)",
                }}
              >
                <div style={{ fontSize: 28, marginBottom: 8, opacity: 0.5 }}>◯</div>
                <div style={{ fontWeight: 500, marginBottom: 4 }}>Plasmid viewer</div>
                <div style={{ fontSize: 12.5, opacity: 0.85, maxWidth: 460 }}>
                  Upload a .gb file or ask Splicify for a CRISPR / cloning design — your sequence will render here, then update in place as guides, primers and assemblies land.
                </div>
              </div>
            )}
          </div>
        </section>





        <footer
          style={{
            marginTop: 48, paddingTop: 16,
            borderTop: "1px solid rgba(219,239,231,0.08)",
            fontSize: 11.5, opacity: 0.55,
            display: "flex", justifyContent: "space-between",
            flexWrap: "wrap", gap: 8,
          }}
        >
          <span>Splicify · agent v2 · Claude Sonnet 4.6</span>
          <span>
            Need the full engine + demos?{" "}
            <a href="https://aiplasmiddesign.com" style={{ color: "var(--brass-400)" }}>
              aiplasmiddesign.com
            </a>
          </span>
        </footer>
      </div>
    </main>
  );
}
