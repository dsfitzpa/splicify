// TypeScript React component - components/GenBankDownloadButton.tsx
// Demonstrates how to wire up a download button that produces a .gb file from seqviz data.
//
// Props:
// - sequence: string of DNA (A/C/G/T/N) -- any non-ACGT characters will be stripped by the generator
// - features: array of seqviz-style features (0-based start, end-exclusive is assumed by default)
// - filename (optional)

import React, { useEffect, useRef, useState } from "react";
import { createGenbankText, convertSeqvizFeatureToGenbank, GenbankFeature } from "../lib/genbank";

type SeqvizFeatureLike = {
  start: number; // seqviz 0-based start
  end: number;   // seqviz 0-based end (exclusive)
  name?: string;
  type?: string;
  strand?: 1 | -1;
  qualifiers?: Record<string, string | string[]>;
};

type Props = {
  sequence: string;
  features?: SeqvizFeatureLike[];
  filename?: string;
  locusName?: string;
  definition?: string;
  circular?: boolean;
  className?: string;
  style?: React.CSSProperties;
};

const stripGbExt = (s: string) => s.replace(/\.(gb|gbk|genbank)$/i, "");
const ensureGbExt = (s: string) => (/\.(gb|gbk|genbank)$/i.test(s) ? s : `${s}.gb`);
const sanitizeFilename = (s: string) => s.replace(/[\/\\\0]/g, "_").trim();

export const GenBankDownloadButton: React.FC<Props> = ({
  sequence,
  features = [],
  filename = "sequence.gb",
  locusName,
  definition,
  circular = true,
  className,
  style,
}) => {
  const [isOpen, setIsOpen] = useState(false);
  const [draftName, setDraftName] = useState(stripGbExt(filename));
  const inputRef = useRef<HTMLInputElement>(null);

  useEffect(() => {
    if (isOpen) {
      setDraftName(stripGbExt(filename));
      // Focus + select the filename so the user can immediately overwrite.
      setTimeout(() => {
        inputRef.current?.focus();
        inputRef.current?.select();
      }, 0);
    }
  }, [isOpen, filename]);

  const writeFile = (finalName: string) => {
    const seq = sequence || "";
    const seqLen = seq.replace(/[^acgtACGTnN]/gi, "").length;

    const gbFeatures: GenbankFeature[] = features.map((f) =>
      convertSeqvizFeatureToGenbank(f, seqLen, { zeroBased: true })
    );

    const gbText = createGenbankText({
      locusName: locusName || "MY_SEQUENCE",
      definition: definition || `Sequence exported from seqviz`,
      sequence: seq,
      features: gbFeatures,
      circular,
      date: new Date(),
    });

    const blob = new Blob([gbText], { type: "text/plain;charset=utf-8" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = finalName;
    document.body.appendChild(a);
    a.click();
    a.remove();
    URL.revokeObjectURL(url);
  };

  const confirmDownload = () => {
    const clean = sanitizeFilename(draftName) || stripGbExt(filename) || "sequence";
    writeFile(ensureGbExt(clean));
    setIsOpen(false);
  };

  return (
    <>
      <button
        type="button"
        onClick={() => setIsOpen(true)}
        className={className}
        style={style}
      >
        Download .gb
      </button>
      {isOpen && (
        <div
          role="dialog"
          aria-modal="true"
          aria-labelledby="gb-download-title"
          onClick={() => setIsOpen(false)}
          style={{
            position: "fixed",
            inset: 0,
            background: "rgba(0,0,0,0.45)",
            display: "flex",
            alignItems: "center",
            justifyContent: "center",
            zIndex: 1000,
          }}
        >
          <div
            onClick={(e) => e.stopPropagation()}
            style={{
              background: "#ffffff",
              color: "#1f2937",
              borderRadius: 12,
              padding: 20,
              width: "min(420px, 92vw)",
              boxShadow: "0 10px 40px rgba(0,0,0,0.25)",
            }}
          >
            <div
              id="gb-download-title"
              style={{ fontWeight: 600, fontSize: 16, marginBottom: 8 }}
            >
              Download GenBank file
            </div>
            <label
              htmlFor="gb-filename-input"
              style={{ display: "block", fontSize: 12, color: "#4b5563", marginBottom: 4 }}
            >
              File name
            </label>
            <div style={{ display: "flex", alignItems: "stretch", gap: 0 }}>
              <input
                id="gb-filename-input"
                ref={inputRef}
                type="text"
                value={draftName}
                onChange={(e) => setDraftName(e.target.value)}
                onKeyDown={(e) => {
                  if (e.key === "Enter") {
                    e.preventDefault();
                    confirmDownload();
                  } else if (e.key === "Escape") {
                    e.preventDefault();
                    setIsOpen(false);
                  }
                }}
                style={{
                  flex: 1,
                  padding: "8px 10px",
                  border: "1px solid #d1d5db",
                  borderRight: "none",
                  borderRadius: "6px 0 0 6px",
                  fontSize: 14,
                  outline: "none",
                }}
              />
              <span
                style={{
                  display: "inline-flex",
                  alignItems: "center",
                  padding: "0 10px",
                  background: "#f3f4f6",
                  border: "1px solid #d1d5db",
                  borderRadius: "0 6px 6px 0",
                  fontSize: 14,
                  color: "#6b7280",
                }}
              >
                .gb
              </span>
            </div>
            <div style={{ display: "flex", justifyContent: "flex-end", gap: 8, marginTop: 16 }}>
              <button
                type="button"
                onClick={() => setIsOpen(false)}
                style={{
                  padding: "6px 12px",
                  background: "transparent",
                  border: "1px solid #d1d5db",
                  borderRadius: 6,
                  cursor: "pointer",
                  fontSize: 14,
                  color: "#374151",
                }}
              >
                Cancel
              </button>
              <button
                type="button"
                onClick={confirmDownload}
                style={{
                  padding: "6px 14px",
                  background: "#105b39",
                  border: "1px solid #105b39",
                  borderRadius: 6,
                  cursor: "pointer",
                  fontSize: 14,
                  color: "#ffffff",
                  fontWeight: 500,
                }}
              >
                Download
              </button>
            </div>
          </div>
        </div>
      )}
    </>
  );
};

export default GenBankDownloadButton;
