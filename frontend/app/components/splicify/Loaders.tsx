"use client";

import React, { useEffect, useMemo, useState } from "react";

export const STAGES = [
  { key: "parse",    label: "Parsing sequence",     hint: "Reading GenBank · detecting topology" },
  { key: "annotate", label: "Annotating features",  hint: "6 reference tiers · feature type and function lookup" },
  { key: "score",    label: "Scoring workflows",    hint: "6 cloning strategies · feasibility ranked" },
  { key: "assemble", label: "Assembling in silico", hint: "Picking least-work, highest-success route" },
  { key: "validate", label: "Validating product",   hint: "Sequence + module + interaction invariants" },
];

type Variant = "light" | "dark";

export function usePipelineProgress(stageSeconds = 2.2, loop = true) {
  const [t, setT] = useState(0);
  useEffect(() => {
    let raf = 0;
    let last = performance.now();
    const total = stageSeconds * STAGES.length;
    const tick = (now: number) => {
      const dt = (now - last) / 1000;
      last = now;
      setT((prev) => {
        const next = prev + dt;
        if (next >= total) return loop ? 0 : total;
        return next;
      });
      raf = requestAnimationFrame(tick);
    };
    raf = requestAnimationFrame(tick);
    return () => cancelAnimationFrame(raf);
  }, [stageSeconds, loop]);
  const total = stageSeconds * STAGES.length;
  const overall = Math.min(t / total, 1);
  const stageIdx = Math.min(Math.floor(t / stageSeconds), STAGES.length - 1);
  const stageFrac = (t - stageIdx * stageSeconds) / stageSeconds;
  return { t, overall, stageIdx, stageFrac };
}

export function StageLabels({
  stageIdx,
  stageFrac,
  overall,
  variant = "light",
}: {
  stageIdx: number;
  stageFrac: number;
  overall: number;
  variant?: Variant;
}) {
  const dark = variant === "dark";
  const fg = dark ? "var(--mint-200)" : "var(--forest-800)";
  const muted = dark ? "rgba(219,239,231,0.55)" : "var(--ink-500)";
  const rule = dark ? "rgba(219,239,231,0.15)" : "rgba(15,85,54,0.14)";
  const accent = "var(--brass-500)";
  return (
    <div style={{ width: "100%", display: "flex", flexDirection: "column", gap: 14 }}>
      <div style={{ display: "flex", alignItems: "baseline", justifyContent: "space-between" }}>
        <div className="splicify-mono" style={{ fontSize: 11, letterSpacing: "0.14em", color: muted, textTransform: "uppercase" }}>
          Step {String(stageIdx + 1).padStart(2, "0")} / {String(STAGES.length).padStart(2, "0")}
        </div>
        <div className="splicify-mono" style={{ fontSize: 11, letterSpacing: "0.14em", color: accent, textTransform: "uppercase", fontVariantNumeric: "tabular-nums" }}>
          {String(Math.round(overall * 100)).padStart(2, "0")}%
        </div>
      </div>
      <div style={{ display: "flex", flexDirection: "column", gap: 8 }}>
        {STAGES.map((s, i) => {
          const state = i < stageIdx ? "done" : i === stageIdx ? "active" : "pending";
          return (
            <div key={s.key} style={{
              display: "grid", gridTemplateColumns: "18px 1fr auto", alignItems: "center",
              gap: 12, padding: "6px 0", borderTop: i === 0 ? "none" : `1px solid ${rule}`,
              opacity: state === "pending" ? 0.42 : 1,
              transition: "opacity 400ms ease",
            }}>
              <StageDot state={state} frac={state === "active" ? stageFrac : 0} accent={accent} fg={fg} muted={muted} />
              <div style={{ display: "flex", flexDirection: "column" }}>
                <div className="splicify-display" style={{
                  fontSize: 15, fontWeight: state === "active" ? 600 : 500, color: fg,
                  letterSpacing: "-0.01em",
                }}>{s.label}</div>
                <div className="splicify-mono" style={{ fontSize: 10.5, color: muted, letterSpacing: "0.02em", marginTop: 2 }}>{s.hint}</div>
              </div>
              <div className="splicify-mono" style={{ fontSize: 10.5, color: state === "done" ? accent : muted, minWidth: 28, textAlign: "right" }}>
                {state === "done" ? "✓" : state === "active" ? "···" : ""}
              </div>
            </div>
          );
        })}
      </div>
    </div>
  );
}

function StageDot({ state, frac, accent, fg, muted }: { state: string; frac: number; accent: string; fg: string; muted: string }) {
  const size = 12;
  if (state === "done") {
    return (
      <div style={{ width: size, height: size, borderRadius: "50%", background: accent, boxShadow: `0 0 0 3px color-mix(in oklch, ${accent} 20%, transparent)` }} />
    );
  }
  if (state === "active") {
    return (
      <svg width={18} height={18} viewBox="0 0 18 18">
        <circle cx={9} cy={9} r={6} fill="none" stroke="color-mix(in oklch, currentColor 18%, transparent)" strokeWidth={1.5} style={{ color: fg }} />
        <circle cx={9} cy={9} r={6}
          fill="none" stroke={accent} strokeWidth={1.8}
          strokeLinecap="round"
          strokeDasharray={2 * Math.PI * 6}
          strokeDashoffset={(1 - frac) * 2 * Math.PI * 6}
          transform="rotate(-90 9 9)" />
      </svg>
    );
  }
  return <div style={{ width: size, height: size, borderRadius: "50%", border: `1.5px dashed ${muted}` }} />;
}

function LoaderShell({ children, bg, fg, variant }: { children: React.ReactNode; bg: string; fg: string; variant: Variant }) {
  const dark = variant === "dark";
  const border = dark ? "rgba(219,239,231,0.1)" : "rgba(15,85,54,0.1)";
  return (
    <div style={{
      background: bg, color: fg,
      borderRadius: 18,
      padding: "22px 24px 20px",
      boxShadow: dark ? "0 30px 60px -30px rgba(0,0,0,0.6)" : "var(--shadow-md)",
      border: `1px solid ${border}`,
      display: "flex", flexDirection: "column", gap: 18,
      minHeight: 460, width: "100%",
    }}>
      {children}
    </div>
  );
}

export function WoolHelixLoader({ variant = "light", pace = 2.4 }: { variant?: Variant; pace?: number }) {
  const { t, overall, stageIdx, stageFrac } = usePipelineProgress(pace);
  const dark = variant === "dark";
  const bg = dark ? "var(--forest-800)" : "var(--cream-50)";
  const fg = dark ? "var(--mint-200)" : "var(--forest-800)";
  const wool = dark ? "var(--mint-200)" : "var(--forest-700)";
  const brass = "var(--brass-500)";

  const W = 360, H = 260;
  const cy = H / 2;
  const turns = 3.2;
  const amp = 42;
  const bases = 90;

  const fill = overall;
  const revealedBases = Math.floor(bases * fill);

  const strandA: [number, number][] = [];
  const strandB: [number, number][] = [];
  for (let i = 0; i <= bases; i++) {
    const x = 40 + i * ((W - 80) / bases);
    const phase = (i / bases) * Math.PI * 2 * turns + t * 1.2;
    const yA = cy + Math.sin(phase) * amp;
    const yB = cy + Math.sin(phase + Math.PI) * amp;
    strandA.push([x, yA]);
    strandB.push([x, yB]);
  }
  const toPath = (pts: [number, number][]) => pts.map((p, i) => `${i === 0 ? "M" : "L"}${p[0].toFixed(1)},${p[1].toFixed(1)}`).join(" ");

  const rungs: [[number, number], [number, number], number][] = [];
  for (let i = 0; i <= revealedBases; i += 2) {
    rungs.push([strandA[i], strandB[i], i]);
  }

  return (
    <LoaderShell bg={bg} fg={fg} variant={variant}>
      <div style={{ position: "relative", flex: 1, minHeight: 260, display: "flex", alignItems: "center", justifyContent: "center" }}>
        <svg width={W} height={H} viewBox={`0 0 ${W} ${H}`} style={{ maxWidth: "100%", height: "auto" }}>
          <defs>
            <radialGradient id="woolGrad" cx="0.5" cy="0.5" r="0.6">
              <stop offset="0" stopColor={dark ? "#f5fbf8" : "#ffffff"} />
              <stop offset="0.6" stopColor={dark ? "var(--mint-200)" : "var(--mint-100)"} />
              <stop offset="1" stopColor={dark ? "var(--mint-200)" : "var(--mint-200)"} stopOpacity="0.3" />
            </radialGradient>
            <linearGradient id="strandFade" x1="0" x2="1">
              <stop offset="0" stopColor={wool} stopOpacity="0.2" />
              <stop offset={Math.min(0.98, fill + 0.05)} stopColor={wool} stopOpacity="1" />
              <stop offset={Math.min(1, fill + 0.08)} stopColor={wool} stopOpacity="0.0" />
            </linearGradient>
          </defs>

          <g transform={`translate(${28 + fill * 10}, ${cy})`}>
            <circle r={28 - fill * 20} fill="url(#woolGrad)" />
            {Array.from({ length: 8 }).map((_, i) => {
              const a = (i / 8) * Math.PI * 2 + t * 0.3;
              const r = 22 - fill * 16;
              return <circle key={i} cx={Math.cos(a) * r} cy={Math.sin(a) * r} r={8 - fill * 6} fill="url(#woolGrad)" opacity={0.8} />;
            })}
            <path d={`M ${22 - fill * 14} 0 Q ${30} ${-10 + Math.sin(t * 2) * 4} ${40 - (28 + fill * 10) + 40} ${strandA[0][1] - cy}`}
              stroke={wool} strokeWidth={1.2} fill="none" opacity={0.6} />
          </g>

          <path d={toPath(strandA.slice(0, revealedBases + 1))} stroke="url(#strandFade)" strokeWidth={2.2} fill="none" strokeLinecap="round" />
          <path d={toPath(strandB.slice(0, revealedBases + 1))} stroke="url(#strandFade)" strokeWidth={2.2} fill="none" strokeLinecap="round" />

          {rungs.map(([a, b, i]) => {
            const phase = (i / bases) * Math.PI * 2 * turns + t * 1.2;
            const inFront = Math.cos(phase) > 0;
            return (
              <line key={i} x1={a[0]} y1={a[1]} x2={b[0]} y2={b[1]}
                stroke={inFront ? brass : wool}
                strokeWidth={inFront ? 1.4 : 0.9}
                strokeOpacity={inFront ? 0.85 : 0.35}
                strokeLinecap="round" />
            );
          })}

          {revealedBases < bases && (
            <circle cx={strandA[revealedBases][0]} cy={cy} r={4}
              fill={brass} opacity={0.9 - 0.4 * Math.sin(t * 6)} />
          )}
        </svg>
      </div>
      <StageLabels stageIdx={stageIdx} stageFrac={stageFrac} overall={overall} variant={variant} />
    </LoaderShell>
  );
}

function easeProgress(overall: number, delay: number, span: number) {
  const local = (overall - delay) / span;
  const x = Math.max(0, Math.min(1, local));
  return 1 - Math.pow(1 - x, 3);
}

function SheepGlyph({ x, y, arrived, color, accent }: { x: number; y: number; arrived: boolean; color: string; accent: string }) {
  return (
    <g transform={`translate(${x} ${y})`}>
      <circle r="5" fill="white" stroke={color} strokeWidth="0.8" />
      <circle cx="-3.5" cy="-2" r="3" fill="white" stroke={color} strokeWidth="0.8" />
      <circle cx="3.5" cy="-2" r="3" fill="white" stroke={color} strokeWidth="0.8" />
      <circle cx="0" cy="2" r="3.2" fill="white" stroke={color} strokeWidth="0.8" />
      <circle cx="5.5" cy="0" r="1.6" fill={color} />
      {arrived && <circle r="8" fill="none" stroke={accent} strokeWidth="0.8" opacity="0.35" />}
    </g>
  );
}

export function FlockRingLoader({ variant = "light", pace = 2.4 }: { variant?: Variant; pace?: number }) {
  const { t, overall, stageIdx, stageFrac } = usePipelineProgress(pace);
  const dark = variant === "dark";
  const bg = dark ? "var(--forest-800)" : "var(--cream-50)";
  const fg = dark ? "var(--mint-200)" : "var(--forest-800)";
  const forest = dark ? "var(--mint-200)" : "var(--forest-700)";
  const brass = "var(--brass-500)";

  const W = 360, H = 260;
  const cx = W / 2, cy = H / 2 + 4;
  const R = 82;
  const N = 14;

  const sheep = useMemo(() => {
    const rand = (i: number) => {
      const x = Math.sin(i * 9301 + 49297) * 233280;
      return x - Math.floor(x);
    };
    return Array.from({ length: N }, (_, i) => ({
      i,
      angle: (i / N) * Math.PI * 2 - Math.PI / 2,
      startX: 40 + rand(i) * (W - 80),
      startY: 30 + rand(i + 100) * (H - 60),
      delay: rand(i + 200) * 0.5,
    }));
  }, []);

  const featureLabels = ["CMV", "EGFP", "P2A", "Cas9", "bGH", "AmpR", "ori", "attR1", "ccdB", "CmR", "attR2", "LTR", "LTR", "WPRE"];

  return (
    <LoaderShell bg={bg} fg={fg} variant={variant}>
      <div style={{ position: "relative", flex: 1, minHeight: 260, display: "flex", alignItems: "center", justifyContent: "center" }}>
        <svg width={W} height={H} viewBox={`0 0 ${W} ${H}`} style={{ maxWidth: "100%", height: "auto" }}>
          <defs>
            <filter id="softGlow" x="-40%" y="-40%" width="180%" height="180%">
              <feGaussianBlur stdDeviation="2.5" />
            </filter>
          </defs>

          <circle cx={cx} cy={cy} r={R} fill="none"
            stroke={forest} strokeOpacity={0.18} strokeWidth={1.2} strokeDasharray="2 4" />

          {sheep.map((s, i) => {
            const nextI = (i + 1) % N;
            const p1 = easeProgress(overall, s.delay, 0.7);
            const p2 = easeProgress(overall, sheep[nextI].delay, 0.7);
            if (p1 < 1 || p2 < 1) return null;
            const a1 = s.angle, a2 = sheep[nextI].angle;
            const x1 = cx + Math.cos(a1) * R, y1 = cy + Math.sin(a1) * R;
            const x2 = cx + Math.cos(a2) * R, y2 = cy + Math.sin(a2) * R;
            return <path key={`seg-${i}`}
              d={`M ${x1} ${y1} A ${R} ${R} 0 0 1 ${x2} ${y2}`}
              stroke={forest} strokeOpacity={0.55} strokeWidth={1.4} fill="none" />;
          })}

          {sheep.map((s, i) => {
            const p = easeProgress(overall, s.delay, 0.7);
            if (p < 1) return null;
            const tickX = cx + Math.cos(s.angle) * (R + 14);
            const tickY = cy + Math.sin(s.angle) * (R + 14);
            const labelX = cx + Math.cos(s.angle) * (R + 22);
            const labelY = cy + Math.sin(s.angle) * (R + 22);
            const anchor: "start" | "end" | "middle" = Math.cos(s.angle) > 0.25 ? "start" : Math.cos(s.angle) < -0.25 ? "end" : "middle";
            return (
              <g key={`lbl-${i}`} opacity={Math.min(1, (p - 1 + 0.2) * 5)}>
                <line x1={cx + Math.cos(s.angle) * R} y1={cy + Math.sin(s.angle) * R}
                  x2={tickX} y2={tickY} stroke={brass} strokeWidth={1} />
                <text x={labelX} y={labelY}
                  fontSize="8.5" textAnchor={anchor} dominantBaseline="middle"
                  fill={forest} style={{ fontFamily: "var(--font-mono-splicify)", letterSpacing: "0.05em" }}>
                  {featureLabels[i]}
                </text>
              </g>
            );
          })}

          {sheep.map((s) => {
            const p = easeProgress(overall, s.delay, 0.7);
            const targetX = cx + Math.cos(s.angle) * R;
            const targetY = cy + Math.sin(s.angle) * R;
            const x = s.startX + (targetX - s.startX) * p;
            const y = s.startY + (targetY - s.startY) * p;
            const arrived = p >= 1;
            const floatY = arrived ? Math.sin(t * 1.5 + s.i) * 0.8 : 0;
            return <SheepGlyph key={s.i} x={x} y={y + floatY} arrived={arrived} color={forest} accent={brass} />;
          })}

          <text x={cx} y={cy - 6} fontSize="10" textAnchor="middle" fill={forest} opacity={0.7}
            style={{ fontFamily: "var(--font-mono-splicify)", letterSpacing: "0.1em" }}>
            pSPLICIFY.ai
          </text>
          <text x={cx} y={cy + 10} fontSize="14" textAnchor="middle" fill={forest} fontWeight="600"
            style={{ fontFamily: "var(--font-display)", letterSpacing: "-0.01em" }}>
            {Math.round(overall * 14800).toLocaleString()} bp
          </text>
        </svg>
      </div>
      <StageLabels stageIdx={stageIdx} stageFrac={stageFrac} overall={overall} variant={variant} />
    </LoaderShell>
  );
}

export function SpinningWheelLoader({ variant = "light", pace = 2.4 }: { variant?: Variant; pace?: number }) {
  const { t, overall, stageIdx, stageFrac } = usePipelineProgress(pace);
  const dark = variant === "dark";
  const bg = dark ? "var(--forest-800)" : "var(--cream-50)";
  const fg = dark ? "var(--mint-200)" : "var(--forest-800)";
  const forest = dark ? "var(--mint-200)" : "var(--forest-700)";
  const brass = "var(--brass-500)";

  const W = 360, H = 260;

  const wx = 92, wy = 150, wR = 62;
  const rot = t * 2.2;

  const sx = 268, sy = 150;
  const spoolR = 14 + overall * 26;

  const feedX = wx + Math.cos(-Math.PI / 4) * wR;
  const feedY = wy + Math.sin(-Math.PI / 4) * wR;
  const midX = (feedX + sx) / 2;
  const sagBase = (sy + feedY) / 2 + 22;
  const sag = sagBase + Math.sin(t * 2) * 2;

  const rings = Math.floor(overall * 18);

  return (
    <LoaderShell bg={bg} fg={fg} variant={variant}>
      <div style={{ position: "relative", flex: 1, minHeight: 260, display: "flex", alignItems: "center", justifyContent: "center" }}>
        <svg width={W} height={H} viewBox={`0 0 ${W} ${H}`} style={{ maxWidth: "100%", height: "auto" }}>
          <defs>
            <linearGradient id="spoolGrad" x1="0" x2="1" y1="0" y2="1">
              <stop offset="0" stopColor={brass} stopOpacity="0.95" />
              <stop offset="1" stopColor="var(--brass-400)" stopOpacity="0.85" />
            </linearGradient>
          </defs>

          <line x1="30" y1="210" x2={W - 30} y2="210" stroke={forest} strokeOpacity="0.14" strokeWidth="1" />

          <g transform={`translate(${wx} ${wy})`}>
            <circle r={wR} fill="none" stroke={forest} strokeOpacity="0.75" strokeWidth="1.6" />
            <circle r={wR - 4} fill="none" stroke={forest} strokeOpacity="0.22" strokeWidth="0.8" />
            <g transform={`rotate(${rot * 180 / Math.PI})`}>
              {Array.from({ length: 10 }).map((_, i) => {
                const a = (i / 10) * Math.PI * 2;
                return <line key={i} x1="0" y1="0"
                  x2={Math.cos(a) * (wR - 6)} y2={Math.sin(a) * (wR - 6)}
                  stroke={forest} strokeOpacity="0.55" strokeWidth="1" />;
              })}
              <circle r="3.5" fill={brass} />
            </g>
            <line x1="0" y1={wR} x2="0" y2="58" stroke={forest} strokeWidth="2" />
          </g>

          <g transform={`translate(${wx - wR - 28} ${wy})`}>
            <ellipse rx={24 - overall * 14} ry={18 - overall * 10} fill={dark ? "var(--mint-200)" : "white"} stroke={forest} strokeWidth="0.8" />
            {Array.from({ length: 5 }).map((_, i) => {
              const a = (i / 5) * Math.PI * 2 + t * 0.3;
              const rr = (20 - overall * 12);
              return <circle key={i}
                cx={Math.cos(a) * rr * 0.6} cy={Math.sin(a) * rr * 0.5}
                r={5 - overall * 2}
                fill={dark ? "var(--mint-200)" : "white"}
                stroke={forest} strokeWidth="0.6" opacity="0.85" />;
            })}
          </g>

          <path
            d={`M ${feedX} ${feedY} Q ${midX} ${sag} ${sx - spoolR} ${sy}`}
            stroke={forest} strokeWidth="1.2" fill="none" opacity="0.55" />
          <path
            d={`M ${feedX} ${feedY} Q ${midX} ${sag} ${sx - spoolR} ${sy}`}
            stroke={brass} strokeWidth="1.8" fill="none"
            strokeDasharray="5 9"
            strokeDashoffset={-t * 40} />

          <g transform={`translate(${sx} ${sy})`}>
            <circle r={spoolR + 6} fill="none" stroke={forest} strokeOpacity="0.2" strokeWidth="1" />
            {Array.from({ length: rings }).map((_, i) => {
              const rr = 14 + i * 1.6;
              return <circle key={i} r={rr}
                fill="none"
                stroke={i % 3 === 0 ? brass : forest}
                strokeOpacity={i % 3 === 0 ? 0.8 : 0.55}
                strokeWidth={i % 3 === 0 ? 1.2 : 0.8} />;
            })}
            <circle r="12" fill="url(#spoolGrad)" />
            <circle r="6" fill={dark ? "var(--forest-800)" : "var(--cream-50)"} />
          </g>

          <text x={wx} y={230} fontSize="9" textAnchor="middle" fill={forest} opacity="0.6"
            style={{ fontFamily: "var(--font-mono-splicify)", letterSpacing: "0.14em" }}>RAW SEQUENCE</text>
          <text x={sx} y={230} fontSize="9" textAnchor="middle" fill={forest} opacity="0.6"
            style={{ fontFamily: "var(--font-mono-splicify)", letterSpacing: "0.14em" }}>PLASMID</text>
        </svg>
      </div>
      <StageLabels stageIdx={stageIdx} stageFrac={stageFrac} overall={overall} variant={variant} />
    </LoaderShell>
  );
}
