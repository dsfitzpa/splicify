import { TopNav } from "../components/splicify/Landing";

export default function AboutPage() {
  return (
    <div
      style={{
        background: "var(--forest-800)",
        color: "var(--mint-200)",
        minHeight: "100vh",
        fontFamily: "var(--font-body)",
      }}
    >
      <TopNav variant="dark" active="about" />
      <main style={{ maxWidth: 880, margin: "0 auto", padding: "56px 48px 96px" }}>
        <div
          className="splicify-mono"
          style={{
            fontSize: 12,
            letterSpacing: "0.18em",
            textTransform: "uppercase",
            color: "var(--brass-400)",
            marginBottom: 28,
          }}
        >
          About Splicify
        </div>
        <h1
          className="splicify-display"
          style={{
            fontSize: 72,
            letterSpacing: "-0.035em",
            fontWeight: 600,
            margin: "0 0 24px",
            lineHeight: 0.98,
          }}
        >
          A quiet tool
          <br />
          for a precise craft.
        </h1>

        <p style={paragraph}>
          Splicify is a quiet, deterministic design tool for molecular cloning. You describe the plasmid you
          want to build in plain language — or upload the GenBank files of what you already have — and Splicify
          classifies the request, resolves named parts against an in-house knowledge base, scores every cloning
          method that could build the target, and runs the chosen workflow end-to-end. The result is a primer
          set, a protocol, an annotated plasmid map, and a workflow trace that documents every decision.
        </p>
        <p style={paragraph}>
          The intent classifier and predesign pipeline are fully deterministic — keyword and regex rules in
          the front; Primer3, SBOL3, and a clean-room six-tier annotation pipeline in the back. Plain-language
          plasmid descriptions are matched to a corpus of {`>`}7,000 LLM-annotated reference plasmids by
          semantic retrieval, then edited deterministically: insertions under 40 bp ride on primer tails,
          insertions of 40 bp and longer become synthesis fragments. An optional LLM orchestrator slot is
          reserved for the cases where the deterministic edit set leaves gaps; today it ships as a no-op so
          every reply is reproducible.
        </p>
        <p style={paragraph}>
          The primer-design algorithm uses Primer3 to calculate primer characteristics and carefully weighs
          the optimal extensions to maximise the probability of successful PCR and assembly — annealing Tm,
          overlap Tm, mispriming, primer-dimer risk, secondary structures, fragment count, and length.
          The result is a full picture of the factors contributing to experimental success.
        </p>

        <h2 className="splicify-display" style={h2}>
          Open-source acknowledgements
        </h2>
        <p style={paragraph}>
          Splicify stands on the shoulders of the scientific software community. We are grateful to the
          authors and maintainers of every project below for making their work openly available.
        </p>
        <h3 className="splicify-display" style={h3}>Software & libraries</h3>
        <ul style={ul}>
          <li><strong>Primer3</strong> — primer design, thermodynamic calculations, hairpin / homodimer scoring.</li>
          <li><strong>SBOL3</strong> (Synthetic Biology Open Language v3) — standardised export of modules and SBO-typed interactions; round-trip via the <code>pySBOL3</code> reference implementation.</li>
          <li><strong>BioPython</strong> — GenBank / FASTA parsing, sequence record manipulation, feature handling.</li>
          <li><strong>Sentence-Transformers</strong> (UKPLab) and the <strong>all-MiniLM-L6-v2</strong> model — embedding plasmid token streams and natural-language descriptions for semantic retrieval.</li>
          <li><strong>HNSWlib</strong> — approximate nearest-neighbour index used during corpus build (runtime queries are brute-force cosine over a NumPy array).</li>
          <li><strong>SeqViz</strong> (Lattice Automation) — interactive DNA sequence visualisation in the linear viewer.</li>
          <li><strong>BLAST+</strong>, <strong>MMseqs2</strong>, <strong>Infernal</strong> — feature search across the six annotation reference tiers.</li>
          <li><strong>FastAPI</strong>, <strong>Next.js</strong>, <strong>React</strong>, <strong>PyTorch</strong> — the application and ML stack.</li>
        </ul>

        <h3 className="splicify-display" style={h3}>Sequence & feature data</h3>
        <ul style={ul}>
          <li><strong>SnapGene</strong> — 1,767 reference plasmid sequences span nine functional families (basic cloning vectors, CRISPR plasmids, fluorescent-protein vectors, Gateway destination / entry vectors, I.M.A.G.E. Consortium plasmids, insect-cell vectors, luciferase vectors, Lucigen vectors, mammalian expression vectors). Only the DNA sequences from the SnapGene-distributed GenBank files were used — not the SnapGene-authored features, maps, or notes; every annotation rendered in Splicify is generated by our own clean-room annotation pipeline. This corpus is the regression set for the annotation pipeline and part of the retrieval corpus for plain-language plasmid design.</li>
          <li><strong>NCBI RefSeq</strong> and <strong>NCBI engineered plasmids</strong> — 41 RefSeq + 5,414 engineered records contribute to the 7,256-plasmid retrieval corpus, with sequence and metadata fetched via Entrez.</li>
          <li><strong>VectorBuilder</strong> — 34 representative vectors plus 26 shorthand description ↔ token pairs that seeded the description-conditioned generative model.</li>
          <li><strong>GenoLIB</strong> — 1,062 main-tier nucleotide features and 706 GenoLIB CDS translations underpin the clean-room feature reference (post-pLannotate, 2026-04-19).</li>
          <li><strong>FPbase</strong> — 721 fluorescent-protein records; identifies and classifies reporter CDSs.</li>
          <li><strong>UniProt / SwissProt</strong> — 66,221 curated PE-1 and whitelisted protein entries for protein-level feature search.</li>
          <li><strong>Rfam</strong> — 1,737 curated families covering riboswitches, ribozymes, cis-elements, and structured non-coding RNAs.</li>
          <li><strong>Gene Ontology — Sequence Ontology (SO)</strong> and <strong>Systems Biology Ontology (SBO)</strong> — role and interaction URIs that flow through to SBOL3 export.</li>
        </ul>

        <h2 className="splicify-display" style={h2}>Contact</h2>
        <p style={paragraph}>
          Splicify was created by Devon Fitzpatrick, with advice on automation and business development from
          Rishij Mewada and help from many friends in the molecular-biology community.
        </p>

        <ContactRow
          heading="General inquiries, product, customer service"
          name="Devon Fitzpatrick"
          email="devon@splicify.ai"
        />
        <ContactRow
          heading="Business and legal inquiries"
          name="Rishij Mewada"
          email="rishij@splicify.ai"
        />
      </main>
    </div>
  );
}

const paragraph: React.CSSProperties = {
  fontSize: 16,
  lineHeight: 1.65,
  color: "rgba(219,239,231,0.86)",
  margin: "0 0 16px",
  textWrap: "pretty",
  maxWidth: "70ch",
};

const h2: React.CSSProperties = {
  fontSize: 26,
  letterSpacing: "-0.02em",
  fontWeight: 600,
  margin: "48px 0 14px",
};

const h3: React.CSSProperties = {
  fontSize: 17,
  letterSpacing: "-0.01em",
  fontWeight: 600,
  margin: "26px 0 10px",
  color: "var(--mint-200)",
};

const ul: React.CSSProperties = {
  margin: "0 0 16px",
  paddingLeft: 22,
  color: "rgba(219,239,231,0.85)",
  lineHeight: 1.6,
  fontSize: 15,
  maxWidth: "70ch",
};

function ContactRow({ heading, name, email }: { heading: string; name: string; email: string }) {
  return (
    <div style={{ marginBottom: 22 }}>
      <h3
        className="splicify-display"
        style={{ fontSize: 16, fontWeight: 600, margin: "0 0 6px", color: "var(--mint-200)" }}
      >
        {heading}
      </h3>
      <div style={{ color: "rgba(219,239,231,0.8)", margin: "0 0 4px" }}>{name}</div>
      <a
        href={`mailto:${email}`}
        style={{
          color: "var(--brass-400)",
          textDecoration: "none",
          borderBottom: "1px solid rgba(230,191,85,0.4)",
        }}
      >
        {email}
      </a>
    </div>
  );
}
