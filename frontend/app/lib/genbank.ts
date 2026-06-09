// Minimal GenBank generator for Next.js/TypeScript
// Place this file at app/lib/genbank.ts so app/components can import "../lib/genbank".

export type GenbankQualifier = Record<string, string | string[]>;

export type GenbankFeature = {
  type?: string; // e.g., "gene", "CDS", "misc_feature"
  start: number; // 1-based inclusive
  end: number; // 1-based inclusive
  strand?: 1 | -1;
  qualifiers?: GenbankQualifier;
  label?: string;
};

function padRight(s: string, width: number) {
  if (s.length >= width) return s;
  return s + " ".repeat(width - s.length);
}

function formatLocation(start: number, end: number, seqLen: number, strand?: 1 | -1): string {
  let loc: string;
  if (start <= end) {
    loc = `${start}..${end}`;
  } else {
    loc = `join(${start}..${seqLen},1..${end})`;
  }
  // Use complement() for reverse strand instead of /strand qualifier
  if (strand === -1) {
    loc = `complement(${loc})`;
  }
  return loc;
}

function formatFeatureLine(type: string, location: string): string {
  // GenBank spec: 5 spaces indent + 16 char feature key field = 21 chars, location starts at column 22
  const left = "     " + padRight(type, 16);
  return `${left}${location}`;
}

function formatQualifiers(qualifiers: GenbankQualifier | undefined): string[] {
  if (!qualifiers) return [];
  const lines: string[] = [];
  const indent = " ".repeat(21);
  for (const [k, v] of Object.entries(qualifiers)) {
    if (Array.isArray(v)) {
      for (const item of v) {
        lines.push(`${indent}/${k}="${escapeQualifier(item)}"`);
      }
    } else {
      lines.push(`${indent}/${k}="${escapeQualifier(String(v))}"`);
    }
  }
  return lines;
}

function escapeQualifier(value: string) {
  return value.replace(/"/g, "'");
}

export function createGenbankText(options: {
  locusName?: string;
  definition?: string;
  accession?: string;
  version?: string;
  sequence: string;
  features?: GenbankFeature[];
  circular?: boolean;
  date?: Date;
}): string {
  const {
    locusName = "UNKNOWN",
    definition = "",
    accession = "",
    version = "",
    sequence,
    features = [],
    circular = false,
    date = new Date(),
  } = options;

  const seq = sequence.replace(/[^acgtACGTnN]/g, "").toLowerCase();
  const seqLen = seq.length;
  const months = ["JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC"];
  const dd = String(date.getUTCDate()).padStart(2, "0");
  const mon = months[date.getUTCMonth()];
  const yyyy = date.getUTCFullYear();
  const gbDate = `${dd}-${mon}-${yyyy}`;

  const molType = "DNA";
  const topology = circular ? "circular" : "linear";
  const locus = `LOCUS       ${padRight(locusName, 16)} ${String(seqLen).padStart(6, " ")} bp    ${molType}     ${padRight(topology, 10)} ${gbDate}`;

  const definitionLine = `DEFINITION  ${definition || "."}`;
  const accessionLine = accession ? `ACCESSION   ${accession}` : "";
  const versionLine = version ? `VERSION     ${version}` : "";

  const sourceLines = [
    `SOURCE      synthetic`,
    `  ORGANISM  synthetic construct`,
  ];

  const featureHeader = "FEATURES             Location/Qualifiers";
  const featureLines: string[] = [featureHeader];

  const sourceLoc = seqLen > 0 ? `1..${seqLen}` : "1..0";
  featureLines.push(formatFeatureLine("source", sourceLoc));
  featureLines.push(...formatQualifiers({
    organism: "synthetic construct",
    mol_type: "genomic DNA",
  }));

  for (const f of features) {
    const type = f.type || (f.label ? "misc_feature" : "feature");
    const loc = formatLocation(f.start, f.end, seqLen, f.strand);
    featureLines.push(formatFeatureLine(type, loc));
    const q: GenbankQualifier = {};
    if (f.label) q["label"] = f.label;
    if (f.qualifiers) {
      for (const [k, v] of Object.entries(f.qualifiers)) q[k] = v;
    }
    // Note: strand is now encoded in location via complement(), not as a qualifier
    featureLines.push(...formatQualifiers(q));
  }

  const originHeader = "ORIGIN";
  const originLines: string[] = [originHeader];
  const lineWidth = 60;
  for (let i = 0; i < seq.length; i += lineWidth) {
    const chunk = seq.slice(i, i + lineWidth);
    const groups: string[] = [];
    for (let j = 0; j < chunk.length; j += 10) {
      groups.push(chunk.slice(j, j + 10));
    }
    const idx = String(i + 1).padStart(9, " ");
    originLines.push(`${idx} ${groups.join(" ")}`);
  }

  const terminator = "//";

  const parts: string[] = [
    locus,
    definitionLine,
  ];
  if (accessionLine) parts.push(accessionLine);
  if (versionLine) parts.push(versionLine);
  parts.push(...sourceLines);
  parts.push("");
  parts.push(...featureLines);
  parts.push("");
  parts.push(...originLines);
  parts.push(terminator);

  return parts.join("\n");
}

/**
 * Convert seqviz-like feature (0-based start, end-exclusive by default) to GenBank 1-based inclusive coords.
 */
export function convertSeqvizFeatureToGenbank(
  seqvizFeature: { start: number; end: number; name?: string; type?: string; strand?: 1 | -1; qualifiers?: GenbankQualifier },
  seqLen: number,
  options?: { zeroBased?: boolean }
): GenbankFeature {
  const zb = options?.zeroBased ?? true;
  let { start, end } = seqvizFeature;
  if (zb) {
    // Convert 0-based start (inclusive) to 1-based
    start = start + 1;
    // end is typically exclusive; for end-exclusive, the 1-based inclusive end equals end
    // (e.g., 0..100 exclusive -> 1..100 inclusive)
  }
  // Normalize into 1..seqLen
  start = ((start - 1 + seqLen) % seqLen) + 1;
  end = ((end - 1 + seqLen) % seqLen) + 1;
  return {
    start,
    end,
    type: seqvizFeature.type,
    strand: seqvizFeature.strand,
    label: seqvizFeature.name,
    qualifiers: seqvizFeature.qualifiers,
  };
}
