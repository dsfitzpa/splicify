import { TopNav } from "../components/splicify/Landing";

const WORKFLOWS: { title: string; blurb: string; beta?: boolean }[] = [
  { title: "Gibson assembly",                blurb: "Build from synthesised or PCR fragments." },
  { title: "Gateway cloning",                blurb: "BP / LR entry and destination vectors." },
  { title: "Golden Gate",                    blurb: "BsaI / BsmBI with scar-free junctions." },
  { title: "sgRNA Golden Gate",              blurb: "CRISPR oligos into lentiCRISPR v2." },
  { title: "Restriction cloning",            blurb: "Insert + vector with site selection." },
  { title: "Site-directed mutagenesis",      blurb: "Point mutations, insertions, deletions." },
  { title: "Plasmid annotation",             blurb: "Six-tier feature identification.", beta: true },
  { title: "Describe-a-plasmid",             blurb: "Natural-language design via semantic retrieval.", beta: true },
];

export default function TutorialPage() {
  return (
    <div style={{
      background: "var(--forest-800)",
      color: "var(--mint-200)",
      minHeight: "100vh",
      fontFamily: "var(--font-body)",
    }}>
      <TopNav variant="dark" active="tutorial" />
      <main style={{ maxWidth: 960, margin: "0 auto", padding: "56px 48px 96px" }}>
        <div className="splicify-mono" style={{
          fontSize: 12, letterSpacing: "0.18em", textTransform: "uppercase",
          color: "var(--brass-400)", marginBottom: 24,
        }}>
          Tutorial
        </div>
        <h1 className="splicify-display" style={{
          fontSize: 72, letterSpacing: "-0.035em", fontWeight: 600,
          margin: "0 0 12px", lineHeight: 0.98,
        }}>
          Prompt Splicify<br />with a design request.
        </h1>
        <p style={{
          fontSize: 18, color: "rgba(219,239,231,0.8)", maxWidth: "65ch",
          marginBottom: 40, lineHeight: 1.65, textWrap: "pretty",
        }}>
          Describe what you want to build in plain language, attach any GenBank files you have,
          and Splicify picks the best workflow. The chat is the front door for designing
          plasmids; the plasmid viewer is where you inspect and edit them. This page walks
          through both.
        </p>

        <h2 className="splicify-display" style={h2}>How a chat request flows</h2>
        <p style={p}>
          Every prompt you send hits a single endpoint, <code>/api/chat</code>, that runs five
          deterministic steps before any cloning workflow executes. You don&rsquo;t need to
          know the internals to use it — but a quick look at the path explains why certain
          phrasings get certain results.
        </p>
        <Card title="The five steps behind every prompt">
          {[
            "Sequence extraction — any DNA pasted in the message is pulled out and replaced with placeholders so the rest of the pipeline can read your intent without choking on long bases.",
            "Intent classification — a deterministic keyword + regex matcher picks one of nine intents (Gibson, Gateway, Golden Gate, sgRNA Golden Gate, restriction, SDM, annotate, describe-a-plasmid, unknown). No LLM is called here.",
            "KB pre-resolution — feature names you mentioned (CMV, eGFP, AmpR, His-tag, …) are matched against the in-house knowledge base, so handlers downstream get pre-resolved sequences instead of guessing.",
            "Unified predesign — for assembly intents, parts are resolved, every part is annotated once (cached by sequence hash), the assembled target is annotated at module-level, and the workflow is scored against the abstract spec of what you asked for.",
            "Per-intent handler — runs the actual designer (Primer3, Gibson, Golden Gate, Gateway operator, …), builds the response, and re-annotates the final plasmid with the full pipeline so the viewer shows every feature.",
          ].map((txt, i) => (
            <div key={i} style={{ display: "flex", gap: 16, alignItems: "flex-start", margin: "10px 0" }}>
              <div className="splicify-mono" style={{
                color: "var(--brass-400)", fontSize: 11, letterSpacing: "0.14em",
                paddingTop: 3, minWidth: 28,
              }}>{String(i + 1).padStart(2, "0")}</div>
              <div style={{ flex: 1, fontSize: 14.5, lineHeight: 1.6, color: "rgba(219,239,231,0.88)" }}>{txt}</div>
            </div>
          ))}
        </Card>

        <h2 className="splicify-display" style={h2}>Phrasing that triggers each workflow</h2>
        <p style={p}>
          The classifier reads your prompt and the files you attach, then picks an intent.
          You don&rsquo;t have to memorise the rules — but if you want a specific workflow,
          phrase the request to match.
        </p>
        <div style={{
          background: "rgba(0,0,0,0.16)",
          border: "1px solid rgba(219,239,231,0.12)",
          borderRadius: 12, padding: "12px 0", margin: "16px 0 28px",
          overflowX: "auto",
        }}>
          <table style={{ width: "100%", borderCollapse: "collapse", fontSize: 13.5 }}>
            <thead>
              <tr style={{ color: "var(--brass-400)" }}>
                <th style={th}>Workflow</th>
                <th style={th}>Triggers</th>
                <th style={th}>What you upload</th>
              </tr>
            </thead>
            <tbody>
              <Row w="Gibson assembly" t="&ldquo;Gibson assembly&rdquo;, &ldquo;assemble these fragments&rdquo;, two or more sequences in the prompt" u="None — fragments are pasted inline" />
              <Row w="Gateway cloning" t="&ldquo;Gateway&rdquo;, &ldquo;BP reaction&rdquo;, &ldquo;LR reaction&rdquo;, &ldquo;attB/attP/attL/attR&rdquo;, &ldquo;pDONR&rdquo;" u="Vector and/or insert .gb files" />
              <Row w="Golden Gate primers" t="&ldquo;Golden Gate primers&rdquo;, &ldquo;BsaI/BsmBI primers&rdquo;, multi-fragment language with &ldquo;+&rdquo; separators" u="Optional — sequences inline" />
              <Row w="sgRNA Golden Gate" t="&ldquo;sgRNA&rdquo; or &ldquo;gRNA&rdquo; + &ldquo;Golden Gate&rdquo; / &ldquo;BsmBI&rdquo; / &ldquo;lentiCRISPR&rdquo; / &ldquo;pX330&rdquo;" u="CRISPR vector .gb file" />
              <Row w="Restriction cloning" t="&ldquo;clone X into Y&rdquo;, &ldquo;restriction enzyme/digest cloning&rdquo;, named Type II enzymes (EcoRI, HindIII, …)" u="Backbone .gb file (optional)" />
              <Row w="Site-directed mutagenesis" t="&ldquo;mutate&rdquo;, &ldquo;delete the X&rdquo;, &ldquo;Y66H&rdquo;, &ldquo;point mutation&rdquo;, &ldquo;remove the NLS&rdquo;" u="Plasmid .gb file (target or inventory)" />
              <Row w="Annotate plasmid" t="&ldquo;annotate&rdquo;, &ldquo;identify features&rdquo;" u="Plasmid .gb file" />
              <Row w="Describe a plasmid" t="&ldquo;make / build / design a vector for X&rdquo;, &ldquo;I need a plasmid that …&rdquo;" u="None" />
            </tbody>
          </table>
        </div>

        <h2 className="splicify-display" style={h2}>Supported workflows</h2>
        <div style={{
          display: "grid", gridTemplateColumns: "repeat(auto-fit, minmax(260px, 1fr))",
          gap: 12, margin: "14px 0 8px",
        }}>
          {WORKFLOWS.map((w) => (
            <div key={w.title} style={{
              background: "rgba(219,239,231,0.04)",
              border: "1px solid rgba(219,239,231,0.1)",
              padding: "14px 16px", borderRadius: 10,
            }}>
              <strong style={{ display: "block", color: "var(--mint-200)", marginBottom: 4, fontSize: 14 }}>
                {w.title}
                {w.beta && <Soon>Beta</Soon>}
              </strong>
              <span style={{ fontSize: 13, color: "rgba(219,239,231,0.65)" }}>{w.blurb}</span>
            </div>
          ))}
        </div>

        <h2 className="splicify-display" style={h2}>Example 1 — Gibson from fragments</h2>
        <p style={p}>
          Paste fragment sequences directly in the prompt. Splicify designs primers with
          homology overhangs at every junction and returns a full assembly + protocol.
        </p>
        <Prompt>
          Design gibson assembly primers for assembling a plasmid with these fragments:
          Frag1: GCCTCCTGCTGGTCCCAAGTTGTGAAATCTTTATCG…,
          Frag2: CTTGCCTTCACCTTTACCCTCGATCGTAAAATCATG…,
          Frag3: TCCACTATTCGAGGCCGTTCGTTAATACTTGTTGCG…,
          frag4: GCCACCATGGTGAGCAAGGGCGAGGAGCTGTTCACC…
        </Prompt>
        <Card title="What you'll get back">
          <ul style={ul}>
            <li>Circular construct visualisation showing fragments, overlaps, and primers</li>
            <li>Primer CSV — sequences, Tm, GC%, secondary structures, quality scores</li>
            <li>Overlap CSV — sequence, length, Tm, quality scores</li>
            <li>Junction table with overlap homology and Tm per junction</li>
            <li>Annotated GenBank of the assembled construct</li>
            <li>Expert explanation — potential issues and optimisation tips</li>
          </ul>
        </Card>

        <h2 className="splicify-display" style={h2}>Example 2 — Golden Gate primers</h2>
        <p style={p}>
          Multi-fragment Golden Gate primer design: Splicify picks orthogonal 4 bp overhangs,
          tails BsaI / BsmBI / BbsI sites onto every primer, and verifies no Type IIs site
          appears inside any fragment.
        </p>
        <Prompt>
          Design Golden Gate primers to assemble EF1a promoter + eGFP + bGH polyA terminator.
        </Prompt>

        <h2 className="splicify-display" style={h2}>Example 3 — Restriction cloning</h2>
        <p style={p}>
          Upload your backbone and ask Splicify to clone an insert into it. The classifier
          recognises &ldquo;clone X into Y&rdquo; phrasing, and naming the enzymes pins the
          digest pair. Without explicit enzymes, Splicify scans the backbone&rsquo;s MCS for
          a unique-cutter pair that is also absent from the insert.
        </p>
        <Prompt>
          Run the test restriction enzyme workflow for cloning GFP into pUC19.
          <br />
          [attached: pUC19_empty.gb]
        </Prompt>
        <Card title="What you'll get back">
          <ul style={ul}>
            <li>Selected enzyme pair, scored on MCS unique-cutter check + insert internal-cut check</li>
            <li>Forward and reverse primers with the chosen sites tailed on</li>
            <li>Synthesis-first manifest when no template exists for the insert (gBlock + ligation, no PCR)</li>
            <li>Diagnostic restriction map of the assembled product</li>
            <li>Annotated GenBank of the assembled plasmid (LLM-annotated, not the upload&rsquo;s original features)</li>
          </ul>
        </Card>

        <h2 className="splicify-display" style={h2}>Example 4 — Site-directed mutagenesis</h2>
        <p style={p}>
          Upload a plasmid and describe the change. SDM understands amino-acid notation
          (Y66H, S65T, D10A), feature deletions (&ldquo;delete the His-tag&rdquo;, &ldquo;remove the NLS&rdquo;),
          terminus insertions (&ldquo;insert FLAG at the N-terminus of eGFP&rdquo;), and position-based
          edits (&ldquo;delete bp 100–150&rdquo;). Before resolving the target, Splicify re-annotates
          the input plasmid with the full pipeline so motifs (His-tag, FLAG, HA, V5, kozak, NLS,
          P2A/T2A/E2A/F2A) and CDS submodules are searchable even when they&rsquo;re not in the
          uploaded file&rsquo;s features.
        </p>
        <Prompt>
          Delete the His-tag from this plasmid.
          <br />
          [attached: pCMV_6xHis_GFP_demo.gb]
        </Prompt>

        <h2 className="splicify-display" style={h2}>Example 5 — sgRNA Golden Gate</h2>
        <p style={p}>
          Annealed-oligo cloning of a 20-bp guide into a CRISPR vector. The handler picks the
          right Type IIs enzyme automatically (BsmBI for lentiCRISPR v2, BbsI for pX330) and
          shows the assembled plasmid annotated with the full pipeline — every feature on the
          backbone is visible, not just the guide cassette.
        </p>
        <Prompt>
          Design oligos to clone gRNA GAGTCCGAGCAGAAGAAGAA (EMX1) into lentiCRISPR v2 using Golden Gate assembly.
          <br />
          [attached: lentiCRISPR_v2_unannotated.gb]
        </Prompt>

        <h2 className="splicify-display" style={h2}>Example 6 — Gateway cloning</h2>
        <p style={p}>
          Splicify scans uploaded plasmids for att sites to decide which is a donor (attP)
          and which carries the insert. The BP/LR product is always emitted as a downloadable
          GenBank, and the response viz shows the recombined plasmid with every feature
          annotated.
        </p>
        <Prompt>
          Design primers for insertion of GFP into pDONR221 using Gateway BP recombination.
          <br />
          [attached: pDONR221.gb]
        </Prompt>

        <h2 className="splicify-display" style={h2}>Example 7 — Describe a plasmid</h2>
        <p style={p}>
          When you describe a plasmid in plain language without naming specific fragments,
          Splicify routes the request to the describe-a-plasmid handler. It builds an
          abstract spec from your description (host, modules, topology), runs a semantic
          search over a 7,256-plasmid corpus to find the closest existing foundation, and
          proposes deterministic edits to make it match your spec — primer-tail edits for
          inserts under 40 bp, synthesis fragments for 40 bp and above.
        </p>
        <Prompt>
          Make a mammalian expression vector for eGFP with PuroR selection.
        </Prompt>
        <Card title="What comes back">
          <ul style={ul}>
            <li><code>&lt;plasmid_id&gt;_foundation.gb</code> — the top-hit plasmid from the corpus, freshly annotated</li>
            <li><code>designed_plasmid.gb</code> — the foundation edited to match your spec, re-annotated</li>
            <li><code>workflow_trace.json</code> — ranked retrieval hits, edit operations with rationale, and any spec gaps deferred to the orchestrator</li>
            <li>Plasmid map of the designed plasmid in the viewer</li>
          </ul>
        </Card>

        <h2 className="splicify-display" style={h2}>The plasmid viewer</h2>
        <p style={p}>
          Every workflow that returns a plasmid renders it in the same viewer panel below the
          chat. The viewer is interactive — drag to select, double-click a feature to inspect
          it, and use the toolbar to add, edit, or import annotations. Selections persist
          across the circular and linear views.
        </p>
        <Card title="Toolbar buttons">
          {[
            ["Upload .gb", "Drop in a GenBank or FASTA file to start a fresh viewer session — no chat round-trip needed."],
            ["Add Annotation (Selected)", "Drag-select a region in either viewer, then click to label it. Annotations are added to the in-memory plasmid; download the updated GenBank with Download GenBank (.gb)."],
            ["Add / Change Sequence", "Insert sequence at the cursor, replace the current selection, or delete it. The viewer re-renders and shifts every annotation past the edit."],
            ["Import Annotation(s)", "Upload a CSV with name + sequence columns. Each sequence is mapped onto the loaded plasmid by exact match (forward + reverse strand, with circular wrap-around). Optional headers — type, location, length, description — show up in the same gene-card popup as KB hits. Max-mismatches knob lets you accept fuzzy hits for noisy guide / primer panels."],
            ["Annotate (LLM + CDS)", "Runs the full annotation pipeline: feature scan across six reference tiers, ORF + CDS submodule resolution, rule-based modules (lentiviral / AAV / Gateway / floxed / Pol III / selection cassettes), Pol II expression cassettes, and 2A polyprotein decomposition."],
            ["Scan Cloning Features", "Adds the cloning-feature layer: every Type II / Type IIs cut site, every Gateway att site, and PCR-warning bars at problematic primer-binding regions. Used by the cloning workflows for unique-cutter selection."],
            ["Analyze Plasmid", "Post-assembly purpose inference — what does this plasmid do? Looks for module conflicts (duplicate nucleases, missing terminators) and produces a plain-English summary."],
            ["Download GenBank (.gb)", "Saves the current sequence + every annotation visible in the viewer — including ones you added in the session and any LLM-annotated features."],
            ["SeqViz / Circular toggle", "Switches between the SeqViz library&rsquo;s detail view and the in-house circular renderer. Both share state; selections survive the toggle."],
          ].map(([title, body], i) => (
            <div key={i} style={{ margin: "12px 0" }}>
              <strong style={{ display: "block", color: "var(--mint-200)", fontSize: 14, marginBottom: 4 }}>
                {title}
              </strong>
              <span style={{ fontSize: 14, color: "rgba(219,239,231,0.78)", lineHeight: 1.55 }}>{body}</span>
            </div>
          ))}
        </Card>

        <Card title="Selecting, editing, navigating">
          <ul style={ul}>
            <li>Drag in either viewer to select a range — the bp count and span are shown next to the toolbar.</li>
            <li>Press <code>Delete</code> with a selection to remove that region from the sequence.</li>
            <li>Double-click any feature to open its gene card: type, location, length, description, KB source, sequence (copy-pasteable).</li>
            <li>Wheel-scroll the linear viewer for fine navigation; click anywhere on the circular to jump.</li>
            <li>Cloning-feature glyphs (Type II / IIs brackets, Gateway att crossovers, PCR warnings) can be filtered by family + cutter count from the Cloning Features menu.</li>
          </ul>
        </Card>

        <h2 className="splicify-display" style={h2}>Tips that change the outcome</h2>
        <ul style={ul}>
          <li>Naming the polymerase (&ldquo;Q5&rdquo;, &ldquo;SuperFi II&rdquo;, &ldquo;KOD One&rdquo;) tunes primer Tm targets and gives you the matching protocol.</li>
          <li>Saying &ldquo;linear&rdquo; vs &ldquo;circular&rdquo; pins Gibson topology; otherwise circular is assumed.</li>
          <li>For SDM, prefer notation that&rsquo;s unambiguous: <code>Y66H in eGFP</code> beats <code>change residue 66 of GFP</code>.</li>
          <li>For restriction cloning, naming both enzymes (&ldquo;EcoRI and HindIII&rdquo;) skips the auto-pick step.</li>
          <li>For describe-a-plasmid, mention the host (&ldquo;mammalian&rdquo;, &ldquo;lentiviral&rdquo;, &ldquo;bacterial&rdquo;) — it biases retrieval toward the right corpus subset.</li>
          <li>Files attached as inventory (multiple .gb) are scored against a target; files attached as target are the thing you&rsquo;re modifying.</li>
        </ul>
      </main>
    </div>
  );
}

const h2: React.CSSProperties = {
  fontSize: 26, letterSpacing: "-0.02em", fontWeight: 600, margin: "48px 0 14px",
};
const p: React.CSSProperties = {
  fontSize: 15.5, lineHeight: 1.65, color: "rgba(219,239,231,0.85)",
  margin: "0 0 14px", maxWidth: "72ch", textWrap: "pretty",
};
const ul: React.CSSProperties = {
  margin: "0 0 14px", paddingLeft: 20, color: "rgba(219,239,231,0.85)",
  lineHeight: 1.55, fontSize: 14.5,
};
const th: React.CSSProperties = {
  textAlign: "left", padding: "8px 16px", fontSize: 11.5, letterSpacing: "0.14em",
  textTransform: "uppercase", borderBottom: "1px solid rgba(219,239,231,0.12)",
};
const td: React.CSSProperties = {
  padding: "10px 16px", borderBottom: "1px solid rgba(219,239,231,0.06)",
  verticalAlign: "top", color: "rgba(219,239,231,0.85)", lineHeight: 1.5,
};

function Row({ w, t, u }: { w: string; t: React.ReactNode; u: string }) {
  return (
    <tr>
      <td style={{ ...td, color: "var(--mint-200)", fontWeight: 500 }}>{w}</td>
      <td style={td}>{t}</td>
      <td style={td}>{u}</td>
    </tr>
  );
}

function Prompt({ children }: { children: React.ReactNode }) {
  return (
    <div className="splicify-mono" style={{
      background: "rgba(0,0,0,0.32)",
      borderLeft: "2px solid var(--brass-500)",
      fontSize: 12.5, lineHeight: 1.55,
      padding: "14px 16px", borderRadius: "0 10px 10px 0",
      margin: "12px 0 20px",
      color: "rgba(219,239,231,0.9)", overflowX: "auto", wordBreak: "break-word",
    }}>{children}</div>
  );
}

function Card({ title, children }: { title: string; children: React.ReactNode }) {
  return (
    <div style={{
      background: "rgba(0,0,0,0.16)",
      border: "1px solid rgba(219,239,231,0.12)",
      borderRadius: 16, padding: "28px 32px", margin: "28px 0",
    }}>
      <h3 className="splicify-display" style={{
        fontSize: 19, fontWeight: 600, margin: "0 0 12px",
        color: "var(--mint-200)", letterSpacing: "-0.015em",
      }}>{title}</h3>
      {children}
    </div>
  );
}

function Soon({ children }: { children: React.ReactNode }) {
  return (
    <span className="splicify-mono" style={{
      display: "inline-block",
      fontSize: 9.5, letterSpacing: "0.15em", textTransform: "uppercase",
      color: "var(--brass-400)",
      border: "1px solid rgba(230,191,85,0.4)",
      padding: "2px 8px", borderRadius: 999, marginLeft: 8,
    }}>{children}</span>
  );
}
