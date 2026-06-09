"use client";

import React, { useState } from "react";
import Link from "next/link";
import { useRouter } from "next/navigation";
import { FlockRingLoader, SpinningWheelLoader, WoolHelixLoader } from "./Loaders";

export const FEATURES = [
  {
    n: "01",
    title: "Annotation engine, clean-room",
    body: "MIT-licensed replacement for pLannotate. Six reference tiers searched in parallel — GenoLIB, FPbase, curated SwissProt, curated Rfam, RefSeq, and a motif DB.",
    stat: "~3.8s",
    statLabel: "14.8 kb lentiCRISPR v2",
  },
  {
    n: "02",
    title: "Hierarchical module detection",
    body: "Not features — architecture. LTR-bounded lentiviral payloads. T-DNA LB/RB. Gateway quartets. EBV oriP + EBNA1. OsTIR1+AID degrons.",
    stat: "~80",
    statLabel: "rules across 30 categories",
  },
  {
    n: "03",
    title: "SBO-typed interactions",
    body: "Every plasmid emits a graph of genetic regulation — promoter → CDS [NLS · Cas9 · P2A · PuroR] → Poly A signal. A plasmid has grammar rules and interactions.",
    stat: "graph",
    statLabel: "first-class design artifact",
  },
  {
    n: "04",
    title: "Cloning-feature scanner",
    body: "Type II + Type IIs sites strand-aware. All 22 canonical Gateway att cores classified. PCR feasibility flags for GC extremes, repeats, palindromes.",
    stat: "~0.09s",
    statLabel: "per plasmid",
  },
  {
    n: "05",
    title: "Target-from-Inventory auto-routing",
    body: "Upload target + freezer inventory. We score ten cloning workflows, pick the method with the least hands-on work and highest wet-lab success, and hand off to the right designer.",
    stat: "10",
    statLabel: "workflows scored in parallel",
  },
  {
    n: "06",
    title: "Plasmid language model",
    body: "Build plasmids de novo and create synthetic fragments that align with your description",
    stat: "closed",
    statLabel: "generation ↔ validation loop",
  },
];

type NavActive = "home" | "engine" | "tutorial" | "about";

export function TopNav({ variant = "dark", active = "home" }: { variant?: "light" | "dark"; active?: NavActive }) {
  const dark = variant === "dark";
  const fg = dark ? "var(--mint-200)" : "var(--forest-800)";
  const muted = dark ? "rgba(219,239,231,0.65)" : "var(--ink-500)";
  const accent = dark ? "var(--brass-400)" : "var(--forest-800)";
  const border = dark ? "rgba(219,239,231,0.1)" : "rgba(15,85,54,0.1)";
  const links: { id: NavActive; label: string; href: string }[] = [
    { id: "home",        label: "Home",         href: "/" },
    { id: "engine",      label: "Dashboard",    href: "/engine" },
    { id: "tutorial",    label: "Tutorial",     href: "/tutorial" },
    { id: "about",       label: "About",        href: "/about" },
  ];
  return (
    <nav style={{
      maxWidth: 1400, margin: "0 auto", padding: "22px 56px",
      display: "flex", alignItems: "center", justifyContent: "space-between",
      borderBottom: `1px solid ${border}`,
    }}>
      <Link href="/" style={{ display: "flex", alignItems: "center", gap: 10, textDecoration: "none" }}>
        {/* eslint-disable-next-line @next/next/no-img-element */}
        <img src="/splicify-logo.png" alt="" style={{ width: 34, height: 34, objectFit: "contain" }} />
        <span className="splicify-display" style={{ fontSize: 20, fontWeight: 600, color: fg, letterSpacing: "-0.015em" }}>Splicify</span>
        <span className="splicify-mono" style={{ fontSize: 10, color: "var(--brass-500)", letterSpacing: "0.15em", marginLeft: 8, padding: "2px 7px", border: "1px solid currentColor", borderRadius: 999 }}>BETA</span>
      </Link>
      <div style={{ display: "flex", gap: 28, fontSize: 13.5 }}>
        {links.map((l) => {
          const isActive = active === l.id;
          return (
            <Link key={l.id} href={l.href} style={{
              color: isActive ? accent : muted,
              textDecoration: "none",
              fontWeight: isActive ? 600 : 400,
              borderBottom: isActive ? `1.5px solid ${accent}` : "1.5px solid transparent",
              paddingBottom: 2,
              transition: "color 140ms ease",
            }}>{l.label}</Link>
          );
        })}
      </div>
      <button style={{
        background: dark ? "var(--brass-500)" : "var(--forest-800)",
        color: dark ? "var(--forest-900)" : "var(--mint-200)",
        border: "none", padding: "10px 18px", borderRadius: 10, whiteSpace: "nowrap",
        fontFamily: "var(--font-display)", fontSize: 13.5, fontWeight: 600, cursor: "pointer",
      }}>Sign in</button>
    </nav>
  );
}

function HeroPromptDark() {
  const [v, setV] = useState("");
  const router = useRouter();
  const submit = () => {
    const q = v ? `?q=${encodeURIComponent(v)}` : "";
    router.push(`/engine${q}`);
  };
  return (
    <div style={{ marginTop: 34, maxWidth: 540 }}>
      <div style={{
        display: "flex", alignItems: "center",
        background: "rgba(219,239,231,0.06)",
        border: "1px solid rgba(219,239,231,0.2)",
        borderRadius: 14, padding: "6px 6px 6px 18px",
      }}>
        <span className="splicify-mono" style={{ color: "var(--brass-400)", fontSize: 14, marginRight: 10 }}>›</span>
        <input
          value={v}
          onChange={(e) => setV(e.target.value)}
          onKeyDown={(e) => { if (e.key === "Enter") submit(); }}
          placeholder="Describe the plasmid you want…"
          style={{
            flex: 1, border: "none", outline: "none",
            fontSize: 15.5, color: "var(--mint-200)",
            fontFamily: "var(--font-body)", padding: "14px 4px",
            background: "transparent",
          }}
        />
        <button onClick={submit} style={{
          background: "var(--brass-500)", color: "var(--forest-900)",
          border: "none", padding: "12px 20px", borderRadius: 10, whiteSpace: "nowrap",
          fontFamily: "var(--font-display)", fontSize: 14, fontWeight: 600,
          cursor: "pointer",
        }}>Design it →</button>
      </div>
      <div className="splicify-mono" style={{ fontSize: 10.5, letterSpacing: "0.14em", textTransform: "uppercase", color: "rgba(219,239,231,0.45)", marginTop: 14 }}>
        upload inventory for reference &nbsp;·&nbsp; optionally add a target for a specific end point
      </div>
    </div>
  );
}

function WorkflowScoreboard() {
  const rows: { name: string; success: number; work: number; pick?: boolean }[] = [
    { name: "Gateway",        success: 0.94, work: 0.25, pick: true },
    { name: "Gibson",         success: 0.88, work: 0.45 },
    { name: "Golden Gate",    success: 0.86, work: 0.40 },
    { name: "sgRNA cloning",  success: 0.82, work: 0.55 },
    { name: "SDM",            success: 0.74, work: 0.62 },
    { name: "Synthesis",      success: 1.00, work: 0.92 },
  ];
  return (
    <div style={{
      background: "rgba(8,56,35,0.6)", border: "1px solid rgba(219,239,231,0.12)",
      borderRadius: 16, padding: "22px 24px", fontFamily: "var(--font-mono-splicify)", fontSize: 12,
    }}>
      <div style={{
        display: "grid", gridTemplateColumns: "1fr 80px 80px 28px",
        gap: 10, paddingBottom: 12, borderBottom: "1px solid rgba(219,239,231,0.14)",
        color: "rgba(219,239,231,0.55)", letterSpacing: "0.12em", textTransform: "uppercase", fontSize: 10,
      }}>
        <span>workflow</span><span style={{ textAlign: "right" }}>success</span><span style={{ textAlign: "right" }}>work</span><span></span>
      </div>
      {rows.map((r) => {
        const score = r.success / (1 + r.work);
        return (
          <div key={r.name} style={{
            display: "grid", gridTemplateColumns: "1fr 80px 80px 28px",
            gap: 10, padding: "9px 0",
            color: r.pick ? "var(--mint-200)" : "rgba(219,239,231,0.72)",
            borderBottom: "1px dashed rgba(219,239,231,0.08)",
            alignItems: "center",
          }}>
            <span style={{ display: "flex", alignItems: "center", gap: 8 }}>
              {r.pick && <span style={{ width: 6, height: 6, borderRadius: "50%", background: "var(--brass-400)" }} />}
              <span style={{ fontWeight: r.pick ? 600 : 400, color: r.pick ? "var(--brass-400)" : "inherit" }}>{r.name}</span>
            </span>
            <span style={{ textAlign: "right", fontVariantNumeric: "tabular-nums" }}>{r.success.toFixed(2)}</span>
            <span style={{ textAlign: "right", fontVariantNumeric: "tabular-nums" }}>{r.work.toFixed(2)}</span>
            <span style={{ textAlign: "right", color: r.pick ? "var(--brass-400)" : "rgba(219,239,231,0.35)", fontSize: 10 }}>
              {r.pick ? "PICK" : score.toFixed(2)}
            </span>
          </div>
        );
      })}
    </div>
  );
}

function GridBackdrop() {
  return (
    <svg style={{ position: "absolute", inset: 0, width: "100%", height: "100%", opacity: 0.05, pointerEvents: "none" }} aria-hidden>
      <defs>
        <pattern id="bg-grid" width="48" height="48" patternUnits="userSpaceOnUse">
          <path d="M 48 0 L 0 0 0 48" fill="none" stroke="var(--mint-200)" strokeWidth="0.5" />
        </pattern>
      </defs>
      <rect width="100%" height="100%" fill="url(#bg-grid)" />
    </svg>
  );
}

export function LandingScientific({
  loaderKind = "flock",
  pace = 2.4,
  showDolly = true,
}: {
  loaderKind?: "flock" | "wheel" | "wool";
  pace?: number;
  showDolly?: boolean;
}) {
  const Loader = loaderKind === "flock" ? FlockRingLoader : loaderKind === "wheel" ? SpinningWheelLoader : WoolHelixLoader;
  return (
    <div style={{
      background: "var(--forest-800)",
      color: "var(--mint-200)",
      fontFamily: "var(--font-body)",
      width: "100%", minHeight: "100%", position: "relative", overflow: "hidden",
    }}>
      <GridBackdrop />
      <TopNav variant="dark" active="home" />

      <section style={{
        maxWidth: 1240, margin: "0 auto", padding: "48px 48px 72px",
        display: "grid", gridTemplateColumns: "1fr 1fr", gap: 48,
      }}>
        <div style={{ display: "flex", flexDirection: "column", justifyContent: "center" }}>
          <div className="splicify-mono" style={{ fontSize: 11, letterSpacing: "0.22em", textTransform: "uppercase", color: "var(--brass-400)", marginBottom: 24 }}>
            ◆ &nbsp;v0.1 &nbsp;/&nbsp; annotation + cloning
          </div>
          <h1 className="splicify-display" style={{
            fontSize: 72, lineHeight: 1, letterSpacing: "-0.035em", margin: 0,
            color: "var(--mint-200)", fontWeight: 600, textWrap: "balance",
          }}>
            A plasmid is<br />
            a <span style={{ color: "var(--brass-400)" }}>validated</span><br />
            functional circuit.
          </h1>
          <p style={{ fontSize: 17, lineHeight: 1.55, color: "rgba(219,239,231,0.78)", maxWidth: 480, marginTop: 28, textWrap: "pretty" }}>
            For 15 years plasmid tools gave us flat feature lists.
            Splicify&nbsp;·&nbsp;Clone emits an interactive module graph
            representation, chooses your parts and workflow, clones your
            plasmid in silico, and validates your product&rsquo;s internal
            grammar.
          </p>
          <HeroPromptDark />
        </div>

        <div style={{ display: "flex", flexDirection: "column", gap: 16, alignItems: "stretch", position: "relative" }}>
          {showDolly && (
            // eslint-disable-next-line @next/next/no-img-element
            <img src="/splicify-logo.png" alt=""
              style={{ position: "absolute", width: 96, height: 96, objectFit: "contain",
                       right: -18, top: -42, zIndex: 2, opacity: 0.95,
                       filter: "drop-shadow(0 14px 28px rgba(0,0,0,0.35))",
                       pointerEvents: "none" }} />
          )}
          <Loader variant="dark" pace={pace} />
          <div style={{
            display: "grid", gridTemplateColumns: "repeat(3, 1fr)", gap: 1,
            background: "rgba(219,239,231,0.1)",
            borderRadius: 14, overflow: "hidden",
            border: "1px solid rgba(219,239,231,0.1)",
          }}>
            {[
              ["14.8 kb", "lentiCRISPR v2 annotated"],
              ["0.09 s", "cloning-feature scan"],
              ["10 / 10", "workflows ranked"],
            ].map(([k, v]) => (
              <div key={k} style={{ background: "var(--forest-800)", padding: "16px 18px" }}>
                <div className="splicify-display" style={{ fontSize: 22, fontWeight: 600, color: "var(--brass-400)", letterSpacing: "-0.015em" }}>{k}</div>
                <div className="splicify-mono" style={{ fontSize: 10.5, color: "rgba(219,239,231,0.7)", letterSpacing: "0.08em", textTransform: "uppercase", marginTop: 4 }}>{v}</div>
              </div>
            ))}
          </div>
        </div>
      </section>

      <section style={{ maxWidth: 1240, margin: "0 auto", padding: "56px 48px 56px" }}>
        <div style={{ display: "flex", alignItems: "baseline", justifyContent: "space-between", marginBottom: 36, borderBottom: "1px solid rgba(219,239,231,0.12)", paddingBottom: 18 }}>
          <div className="splicify-mono" style={{ fontSize: 11, letterSpacing: "0.22em", textTransform: "uppercase", color: "var(--brass-400)" }}>
            / What shipped
          </div>
          <div className="splicify-mono" style={{ fontSize: 11, letterSpacing: "0.18em", textTransform: "uppercase", color: "rgba(219,239,231,0.55)" }}>
            6 modules · all first-class
          </div>
        </div>

        <div style={{ display: "grid", gridTemplateColumns: "repeat(3, 1fr)", gap: 1, background: "rgba(219,239,231,0.12)" }}>
          {FEATURES.map((f) => (
            <article key={f.n} style={{
              background: "var(--forest-800)", padding: "28px 26px 26px",
              display: "flex", flexDirection: "column", gap: 14, minHeight: 260,
              position: "relative",
            }}>
              <div style={{ display: "flex", justifyContent: "space-between", alignItems: "baseline" }}>
                <span className="splicify-mono" style={{ fontSize: 11, color: "var(--brass-400)", letterSpacing: "0.15em" }}>
                  {f.n}
                </span>
                <span className="splicify-mono" style={{ fontSize: 10, color: "rgba(219,239,231,0.45)", letterSpacing: "0.1em" }}>
                  /splicify/core
                </span>
              </div>
              <h3 className="splicify-display" style={{ fontSize: 20, letterSpacing: "-0.015em", margin: 0, color: "var(--mint-200)", fontWeight: 600, lineHeight: 1.15 }}>
                {f.title}
              </h3>
              <p style={{ fontSize: 13.5, lineHeight: 1.55, color: "rgba(219,239,231,0.72)", flex: 1, textWrap: "pretty" }}>
                {f.body}
              </p>
              <div style={{ display: "flex", alignItems: "baseline", gap: 10, paddingTop: 12, borderTop: "1px dashed rgba(219,239,231,0.18)" }}>
                <span className="splicify-display" style={{ fontSize: 22, color: "var(--brass-400)", fontWeight: 600, letterSpacing: "-0.02em" }}>{f.stat}</span>
                <span className="splicify-mono" style={{ fontSize: 10, color: "rgba(219,239,231,0.55)", letterSpacing: "0.08em", textTransform: "uppercase" }}>{f.statLabel}</span>
              </div>
            </article>
          ))}
        </div>
      </section>

      <section style={{ maxWidth: 1240, margin: "0 auto", padding: "56px 48px 96px" }}>
        <div style={{ display: "grid", gridTemplateColumns: "1fr 1.2fr", gap: 56, alignItems: "center" }}>
          <div>
            <div className="splicify-mono" style={{ fontSize: 11, letterSpacing: "0.22em", textTransform: "uppercase", color: "var(--brass-400)", marginBottom: 18 }}>
              / Target-from-inventory
            </div>
            <h2 className="splicify-display" style={{ fontSize: 40, letterSpacing: "-0.025em", lineHeight: 1.05, margin: 0, color: "var(--mint-200)", fontWeight: 500 }}>
              Drop in a target. Drop in your freezer. We pick the path.
            </h2>
            <p style={{ fontSize: 15, lineHeight: 1.6, color: "rgba(219,239,231,0.75)", marginTop: 20, maxWidth: 480, textWrap: "pretty" }}>
              Six workflows are scored in parallel — Gateway (BP/LR and
              MultiSite), Gibson Assembly, Golden Gate (fragments and gRNA
              cloning), restriction cloning, SDM, and synthesis. We rank
              by <span className="splicify-mono" style={{ color: "var(--brass-400)" }}>success_estimate / (1 + work_estimate)</span>,
              hand off to the right designer with pre-resolved arguments,
              and validate the in-silico product against your target
              to confirm it preserves function.
            </p>
          </div>
          <WorkflowScoreboard />
        </div>
      </section>

      <section style={{ borderTop: "1px solid rgba(219,239,231,0.1)", padding: "28px 48px", maxWidth: 1240, margin: "0 auto", display: "flex", justifyContent: "space-between", fontSize: 12, color: "rgba(219,239,231,0.55)" }}>
        <span className="splicify-mono" style={{ letterSpacing: "0.12em" }}>© SPLICIFY · DEVON@SPLICIFY.AI</span>
        <span className="splicify-mono" style={{ letterSpacing: "0.12em" }}>MIT-LICENSED · BUILT FOR WET-LAB SCIENTISTS</span>
      </section>
    </div>
  );
}
