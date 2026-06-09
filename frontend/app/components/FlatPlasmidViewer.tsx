"use client";

import React, {
  useCallback,
  useEffect,
  useMemo,
  useRef,
  useState,
} from "react";

import type { SeqVizAnnotation } from "./InteractiveSequenceViewer";

// Single-track horizontal whole-sequence viewer. Whole plasmid lays out
// left-to-right at zoom=1; the user zooms in with the mouse wheel and pans
// by dragging on the empty track. Click a feature -> onAnnotationClick.
// Drag the empty track with the mouse -> onSelectionChange.

type Props = {
  sequence: string;
  annotations: SeqVizAnnotation[];
  height?: number;
  selection?: { start: number; end: number } | null;
  onSelectionChange?: (sel: { start: number; end: number } | null) => void;
  onAnnotationClick?: (ann: SeqVizAnnotation, idx: number) => void;
  /** Fired on single-click of a feature: parent should select + center the linear viewer on it. */
  onAnnotationSelect?: (sel: { start: number; end: number }) => void;
};

// Vivid palette for features the design pipeline injected (pegRNAs,
// ngRNAs, primers, edit sites). Cycles by index so each design-added
// feature gets its own hue and pops against the muted preserved
// annotations. Matches CircularPlasmidViewer + LinearSequenceViewer.
const DESIGN_PALETTE = [
  "#ec4899", // magenta
  "#f59e0b", // amber
  "#06b6d4", // cyan
  "#8b5cf6", // violet
  "#22c55e", // lime
  "#ef4444", // red
  "#14b8a6", // teal
  "#f97316", // orange
];

const PRESERVED_GRAY_PALETTE = ["#9ca3af", "#6b7280", "#cbd5e1"];

const FEATURE_ALT_PALETTE = [
  "#4a7c59",
  "#6b5b95",
  "#3b82f6",
  "#f97316",
  "#22c55e",
  "#ef4444",
];

function featureColor(
  ann: SeqVizAnnotation,
  index: number = 0,
  hasDesignAdditions: boolean = false,
): string {
  if ((ann as any).added_by_design) {
    return DESIGN_PALETTE[index % DESIGN_PALETTE.length];
  }
  if (hasDesignAdditions) {
    return PRESERVED_GRAY_PALETTE[index % PRESERVED_GRAY_PALETTE.length];
  }
  if (ann.color) return ann.color;
  const layer = (ann as any).layer as string | undefined;
  if (layer === "module") return "#4a7c59";
  if (layer === "motif") return "#6b5b95";
  if (layer === "gap") return "#5a6c7d";
  if (layer === "cloning_feature") return "#C2185B";
  return FEATURE_ALT_PALETTE[index % FEATURE_ALT_PALETTE.length];
}

function fmtBp(n: number): string {
  if (n >= 1_000_000) return `${(n / 1_000_000).toFixed(1)}M`;
  if (n >= 1_000) return `${(n / 1_000).toFixed(n >= 10_000 ? 0 : 1)}k`;
  return String(n);
}

function spanOf(ann: SeqVizAnnotation, L: number): { s: number; e: number; spans: boolean } {
  // start <= end : normal segment [start, end)
  // start  > end : origin-spanning -> [start, L) + [0, end)
  if (ann.start <= ann.end) return { s: ann.start, e: ann.end, spans: false };
  return { s: ann.start, e: ann.end, spans: true };
}

// Pack features into lanes (row index) using a simple greedy interval
// algorithm. We expand origin-spanning features in advance.
type Seg = {
  origIdx: number;
  ann: SeqVizAnnotation;
  s: number;
  e: number;     // exclusive
};

function packLanes(segs: Seg[]): Seg[][] {
  const sorted = [...segs].sort((a, b) => a.s - b.s);
  const lanes: Seg[][] = [];
  for (const seg of sorted) {
    let placed = false;
    for (const lane of lanes) {
      const last = lane[lane.length - 1];
      if (seg.s >= last.e) {
        lane.push(seg);
        placed = true;
        break;
      }
    }
    if (!placed) lanes.push([seg]);
  }
  return lanes;
}

const MIN_ZOOM = 1;
const MIN_VISIBLE_BP = 10;  // max-zoom floor: never less than 10 bp in view
const BASE_HEIGHT = 28;
const HUD_HEIGHT = 32;
const RULER_HEIGHT = 18;
const FEATURE_HEIGHT = 16;
const FEATURE_GAP = 4;
const SEQUENCE_THRESHOLD_PX = 9;   // base width at which we start to draw letters

const COLOR_TRACK = "#ffffff";
const COLOR_TEXT_DARK = "#105b39";
const COLOR_TICK = "rgba(16,91,57,0.45)";
const COLOR_TICK_LABEL = "#105b39";
const COLOR_SELECTION = "rgba(59, 130, 246, 0.20)";
const COLOR_SELECTION_BORDER = "#3b82f6";
const COLOR_LABEL = "#105b39";
const COLOR_FEATURE_LABEL = "#ffffff";
const SEQ_FONT_FAMILY = "'Droid Sans Mono', 'Courier New', Courier, monospace";
const SANS_FONT_FAMILY = "'Source Sans Pro', sans-serif";
const FEATURE_OPACITY = 0.5;

export default function FlatPlasmidViewer({
  sequence,
  annotations,
  height = 460,
  selection,
  onSelectionChange,
  onAnnotationClick,
  onAnnotationSelect,
}: Props) {
  const wrapRef = useRef<HTMLDivElement | null>(null);
  const [containerWidth, setContainerWidth] = useState<number>(800);
  const [zoom, setZoom] = useState<number>(1);
  const [pan, setPan] = useState<number>(0);
  const panRef = useRef(0);
  const zoomRef = useRef(1);
  const containerWidthRef = useRef(800);
  const sequenceLenRef = useRef(0);

  const L = sequence ? sequence.length : 0;

  // Whether the .gb on screen is a design-output file (contains at least
  // one /added_by feature). When true, featureColor gray-mutes preserved
  // annotations so the vivid design overlay reads at a glance.
  const hasDesignAdditions = useMemo(
    () => annotations.some((a) => (a as any).added_by_design === true),
    [annotations],
  );

  // Watch container size to fit the full sequence at zoom=1.
  useEffect(() => {
    const el = wrapRef.current;
    if (!el) return;
    const obs = new ResizeObserver((entries) => {
      for (const e of entries) {
        const w = Math.max(200, e.contentRect.width);
        setContainerWidth(w);
      }
    });
    obs.observe(el);
    setContainerWidth(Math.max(200, el.getBoundingClientRect().width));
    return () => obs.disconnect();
  }, []);

  // Native (non-passive) wheel listener so preventDefault() actually
  // stops the page from scrolling underneath the zoom.
  useEffect(() => {
    const el = wrapRef.current;
    if (!el) return;
    const onWheel = (ev: WheelEvent) => {
      ev.preventDefault();
      const rect = el.getBoundingClientRect();
      const cursorViewX = ev.clientX - rect.left;
      const cursorContentX = panRef.current + cursorViewX;
      const cw = containerWidthRef.current;
      const lenRef = sequenceLenRef.current;
      const zoomCur = zoomRef.current;
      const contentW = cw * zoomCur;
      const cursorBp = lenRef > 0
        ? Math.max(0, Math.min(lenRef, (cursorContentX / contentW) * lenRef))
        : 0;
      const dir = -Math.sign(ev.deltaY);
      const factor = dir > 0 ? 1.2 : 1 / 1.2;
      const maxZoomDyn = lenRef > 0 ? Math.max(MIN_ZOOM, lenRef / MIN_VISIBLE_BP) : MIN_ZOOM;
      const newZoom = Math.max(MIN_ZOOM, Math.min(maxZoomDyn, zoomCur * factor));
      const newContentWidth = cw * newZoom;
      const newCursorContentX = lenRef > 0 ? (cursorBp / lenRef) * newContentWidth : 0;
      const newPan = Math.max(
        0,
        Math.min(Math.max(0, newContentWidth - cw), newCursorContentX - cursorViewX)
      );
      setZoom(newZoom);
      setPan(newPan);
    };
    el.addEventListener("wheel", onWheel, { passive: false });
    return () => el.removeEventListener("wheel", onWheel);
  }, []);

  useEffect(() => { panRef.current = pan; }, [pan]);
  useEffect(() => { zoomRef.current = zoom; }, [zoom]);
  useEffect(() => { containerWidthRef.current = containerWidth; }, [containerWidth]);
  useEffect(() => { sequenceLenRef.current = sequence ? sequence.length : 0; }, [sequence]);

  // Geometry: at zoom = z, total content width = containerWidth * z.
  // pan is content-space x offset of the leftmost visible pixel.
  const contentWidth = containerWidth * zoom;
  const bpPerPx = L > 0 ? L / contentWidth : 1;
  const pxPerBp = L > 0 ? contentWidth / L : 1;

  const visibleStartBp = Math.max(0, Math.floor(pan * bpPerPx));
  const visibleEndBp = Math.min(
    L,
    Math.ceil((pan + containerWidth) * bpPerPx) + 1
  );

  const clampPan = useCallback(
    (p: number) => Math.max(0, Math.min(contentWidth - containerWidth, p)),
    [contentWidth, containerWidth]
  );

  useEffect(() => {
    // Re-clamp pan when zoom or width changes (e.g. window resize zooms out).
    setPan((p) => Math.max(0, Math.min(Math.max(0, contentWidth - containerWidth), p)));
  }, [contentWidth, containerWidth]);

  // Build expanded segments (origin-spanning features split in two).
  const segs = useMemo<Seg[]>(() => {
    if (L === 0) return [];
    const out: Seg[] = [];
    annotations.forEach((ann, origIdx) => {
      const { s, e, spans } = spanOf(ann, L);
      if (!spans) {
        out.push({ origIdx, ann, s, e });
      } else {
        out.push({ origIdx, ann, s: s, e: L });
        out.push({ origIdx, ann, s: 0, e: e });
      }
    });
    return out;
  }, [annotations, L]);

  // Pack into lanes for stacking.
  const lanes = useMemo(() => packLanes(segs), [segs]);
  const featureLanes = Math.max(1, lanes.length);
  const trackHeight =
    HUD_HEIGHT +
    RULER_HEIGHT +
    BASE_HEIGHT +                                 // hud + ruler tick row + sequence row
    featureLanes * (FEATURE_HEIGHT + FEATURE_GAP);

  const svgHeight = Math.max(height, trackHeight + 20);

  // Bp <-> px helpers (content-space).
  const bpToContentX = useCallback(
    (bp: number) => (L > 0 ? (bp / L) * contentWidth : 0),
    [L, contentWidth]
  );
  const contentXToBp = useCallback(
    (x: number) => Math.max(0, Math.min(L, Math.round((x / contentWidth) * L))),
    [L, contentWidth]
  );

  // Mouse wheel: zoom around cursor.
  const handleWheel = useCallback(
    (e: React.WheelEvent<HTMLDivElement>) => {
      e.preventDefault();
      const rect = wrapRef.current?.getBoundingClientRect();
      if (!rect) return;
      const cursorViewX = e.clientX - rect.left;
      const cursorContentX = pan + cursorViewX;
      const cursorBp = contentXToBp(cursorContentX);
      const dir = -Math.sign(e.deltaY);
      const factor = dir > 0 ? 1.2 : 1 / 1.2;
      const maxZoomDyn = L > 0 ? Math.max(MIN_ZOOM, L / MIN_VISIBLE_BP) : MIN_ZOOM;
      const newZoom = Math.max(MIN_ZOOM, Math.min(maxZoomDyn, zoom * factor));
      const newContentWidth = containerWidth * newZoom;
      const newCursorContentX = (cursorBp / Math.max(1, L)) * newContentWidth;
      const newPan = Math.max(
        0,
        Math.min(
          Math.max(0, newContentWidth - containerWidth),
          newCursorContentX - cursorViewX
        )
      );
      setZoom(newZoom);
      setPan(newPan);
    },
    [zoom, pan, contentXToBp, containerWidth, L]
  );

  // Pointer interactions: shift = pan, plain drag = select.
  const dragRef = useRef<{
    mode: "pan" | "select" | null;
    startViewX: number;
    startPan: number;
    startBp: number;
  }>({ mode: null, startViewX: 0, startPan: 0, startBp: 0 });

  const localViewX = useCallback(
    (clientX: number) => {
      const rect = wrapRef.current?.getBoundingClientRect();
      if (!rect) return 0;
      return clientX - rect.left;
    },
    []
  );

  const handlePointerDown = useCallback(
    (e: React.PointerEvent<HTMLDivElement>) => {
      const tgt = e.target as Element | null;
      // Buttons + interactive elements own their own click events. Also
      // bail out for anything inside the zoom-controls overlay so a click
      // on the ±/Fit/range pill never starts a selection drag.
      if (
        tgt && (
          tgt.closest("button") ||
          tgt.closest("input") ||
          tgt.closest("select") ||
          tgt.closest('[data-ui="zoom-controls"]')
        )
      ) return;
      const vx = localViewX(e.clientX);
      const bp = contentXToBp(pan + vx);
      const isFeatureClick = tgt?.getAttribute?.("data-feat") === "1";
      if (isFeatureClick) return;
      const mode: "pan" | "select" = e.shiftKey || e.button === 1 ? "pan" : "select";
      dragRef.current = {
        mode,
        startViewX: vx,
        startPan: pan,
        startBp: bp,
      };
      (e.currentTarget as Element).setPointerCapture(e.pointerId);
      if (mode === "select") {
        onSelectionChange?.({ start: bp, end: bp });
      }
      e.preventDefault();
    },
    [pan, contentXToBp, localViewX, onSelectionChange]
  );

  const handlePointerMove = useCallback(
    (e: React.PointerEvent<HTMLDivElement>) => {
      const d = dragRef.current;
      if (!d.mode) return;
      const vx = localViewX(e.clientX);
      if (d.mode === "pan") {
        const dx = vx - d.startViewX;
        setPan(clampPan(d.startPan - dx));
      } else {
        const bp = contentXToBp(pan + vx);
        onSelectionChange?.({ start: d.startBp, end: bp });
      }
    },
    [pan, contentXToBp, localViewX, onSelectionChange, clampPan]
  );

  const handlePointerUp = useCallback(
    (e: React.PointerEvent<HTMLDivElement>) => {
      try {
        (e.currentTarget as Element).releasePointerCapture(e.pointerId);
      } catch {}
      dragRef.current.mode = null;
    },
    []
  );

  // Ruler ticks: pick a tick step that yields ~6-10 labels in the visible window.
  const visBp = Math.max(1, visibleEndBp - visibleStartBp);
  const targetTicks = Math.max(4, Math.min(12, Math.round(containerWidth / 110)));
  let step = Math.pow(10, Math.floor(Math.log10(visBp / targetTicks)));
  // Snap step to {1,2,5} * 10^k
  for (const m of [1, 2, 5, 10]) {
    if (m * step >= visBp / targetTicks) {
      step = m * step;
      break;
    }
  }
  const firstTick = Math.ceil(visibleStartBp / step) * step;
  const ticks: number[] = [];
  for (let t = firstTick; t <= visibleEndBp; t += step) ticks.push(t);

  // Selection rendering.
  const sel = selection || null;
  const selA = sel ? Math.min(sel.start, sel.end) : null;
  const selB = sel ? Math.max(sel.start, sel.end) : null;
  const selVisible =
    selA !== null && selB !== null && selB > selA &&
    selB >= visibleStartBp && selA <= visibleEndBp;

  // Choose lane positioning helpers.
  const featuresTopY = HUD_HEIGHT + RULER_HEIGHT + BASE_HEIGHT;

  // For each lane, render features that are visible.
  const renderLanes = lanes.map((laneSegs, laneIdx) => {
    const items: React.ReactNode[] = [];
    laneSegs.forEach((seg) => {
      // Visibility: any overlap with [visibleStartBp, visibleEndBp]
      if (seg.e < visibleStartBp || seg.s > visibleEndBp) return;
      const xContentStart = bpToContentX(seg.s);
      const xContentEnd = bpToContentX(seg.e);
      const xViewStart = xContentStart - pan;
      const xViewEnd = xContentEnd - pan;
      const widthRaw = Math.max(1, xViewEnd - xViewStart);
      const ann = seg.ann;
      const dir = ann.direction === -1 ? -1 : 1;
      const color = featureColor(ann, seg.origIdx, hasDesignAdditions);
      const y = featuresTopY + laneIdx * (FEATURE_HEIGHT + FEATURE_GAP);
      const showLabel = widthRaw > 32;
      // Draw a notched arrow when wide enough; otherwise just a rectangle.
      const arrowSize = 6;
      const w = widthRaw;
      let path: string;
      if (w > arrowSize * 2) {
        if (dir === 1) {
          path = `M ${xViewStart} ${y}
                  L ${xViewStart + w - arrowSize} ${y}
                  L ${xViewStart + w} ${y + FEATURE_HEIGHT / 2}
                  L ${xViewStart + w - arrowSize} ${y + FEATURE_HEIGHT}
                  L ${xViewStart} ${y + FEATURE_HEIGHT} Z`;
        } else {
          path = `M ${xViewStart + w} ${y}
                  L ${xViewStart + arrowSize} ${y}
                  L ${xViewStart} ${y + FEATURE_HEIGHT / 2}
                  L ${xViewStart + arrowSize} ${y + FEATURE_HEIGHT}
                  L ${xViewStart + w} ${y + FEATURE_HEIGHT} Z`;
        }
      } else {
        path = `M ${xViewStart} ${y}
                L ${xViewStart + w} ${y}
                L ${xViewStart + w} ${y + FEATURE_HEIGHT}
                L ${xViewStart} ${y + FEATURE_HEIGHT} Z`;
      }
      items.push(
        <g
          key={`${seg.origIdx}-${seg.s}-${seg.e}-${laneIdx}`}
          onClick={(ev) => {
            ev.stopPropagation();
            // Single click matches CircularPlasmidViewer's behavior: select +
            // center the linear viewer on the feature. Double click opens
            // the edit-annotation modal.
            onAnnotationSelect?.({ start: seg.ann.start, end: seg.ann.end });
          }}
          onDoubleClick={(ev) => {
            ev.stopPropagation();
            onAnnotationClick?.(seg.ann, seg.origIdx);
          }}
        >
          <path
            d={path}
            data-feat="1"
            fill={color}
            opacity={FEATURE_OPACITY}
            stroke={color}
            strokeOpacity={0.85}
            strokeWidth={0.75}
            style={{ cursor: "pointer" }}
          />
          {showLabel && (() => {
            const visStart = Math.max(0, xViewStart);
            const visEnd = Math.min(containerWidth, xViewStart + w);
            const labelMid = (visStart + visEnd) / 2;
            const labelMaxChars = Math.max(2, Math.floor((visEnd - visStart - 8) / 6));
            const labelText =
              seg.ann.name && seg.ann.name.length > labelMaxChars
                ? seg.ann.name.slice(0, Math.max(1, labelMaxChars - 1)) + "\u2026"
                : seg.ann.name;
            return (
              <text
                x={labelMid}
                y={y + FEATURE_HEIGHT - 4}
                fontSize={10}
                textAnchor="middle"
                fill={COLOR_FEATURE_LABEL}
                fontFamily={SANS_FONT_FAMILY}
                style={{ pointerEvents: "none", userSelect: "none" }}
              >
                {labelText}
              </text>
            );
          })()}
        </g>
      );
    });
    return items;
  });

  // Sequence letters when zoomed in enough.
  const drawSequenceLetters = pxPerBp >= SEQUENCE_THRESHOLD_PX && L > 0;
  let letters: React.ReactNode = null;
  if (drawSequenceLetters) {
    const startBp = Math.max(0, visibleStartBp - 1);
    const endBp = Math.min(L, visibleEndBp + 1);
    const arr: React.ReactNode[] = [];
    for (let i = startBp; i < endBp; i++) {
      const xMid = bpToContentX(i + 0.5) - pan;
      arr.push(
        <text
          key={`l-${i}`}
          x={xMid}
          y={HUD_HEIGHT + RULER_HEIGHT + BASE_HEIGHT - 6}
          fontSize={11}
          textAnchor="middle"
          fill={COLOR_TEXT_DARK}
          fontFamily={SEQ_FONT_FAMILY}
          style={{ pointerEvents: "none", userSelect: "none" }}
        >
          {sequence[i] || ""}
        </text>
      );
    }
    letters = arr;
  }

  // Render.
  // Zoom by a multiplicative factor (> 1 zooms in, < 1 zooms out).
  // Anchor:
  //   - zoom-in with an active selection: midpoint of the selection
  //   - otherwise: current visible-window midpoint
  // Pan is recomputed so the anchor bp lands at the horizontal centre of
  // the viewer, which keeps the user's focus stable across +/- clicks.
  const handleZoomBy = (factor: number) => {
    if (L <= 0 || containerWidth <= 0) return;
    const maxZoomDyn = Math.max(MIN_ZOOM, L / MIN_VISIBLE_BP);
    const newZoom = Math.max(MIN_ZOOM, Math.min(maxZoomDyn, zoom * factor));
    const visStartBp = (pan / contentWidth) * L;
    const visEndBp = ((pan + containerWidth) / contentWidth) * L;
    const viewCentreBp = (visStartBp + visEndBp) / 2;
    let anchorBp = viewCentreBp;
    if (factor > 1 && selection && selection.start !== selection.end) {
      const a = Math.min(selection.start, selection.end);
      const b = Math.max(selection.start, selection.end);
      anchorBp = (a + b) / 2;
    }
    const newContentWidth = containerWidth * newZoom;
    const anchorContentX = (anchorBp / L) * newContentWidth;
    const newPan = Math.max(
      0,
      Math.min(
        Math.max(0, newContentWidth - containerWidth),
        anchorContentX - containerWidth / 2
      )
    );
    setZoom(newZoom);
    setPan(newPan);
  };

  const handleResetView = () => {
    setZoom(1);
    setPan(0);
  };

  return (
    <div
      ref={wrapRef}
      style={{
        position: "relative",
        width: "100%",
        height,
        overflow: "hidden",
        background: COLOR_TRACK,
        borderRadius: 8,
      }}
      onPointerDown={handlePointerDown}
      onPointerMove={handlePointerMove}
      onPointerUp={handlePointerUp}
      onPointerCancel={handlePointerUp}
    >
      {/* Top bar: zoom controls + visible range.
          Stops pointer events from bubbling to the SVG wrapper so the
          ±/Fit buttons and the surrounding pill never start a selection
          drag on the track underneath — fixes double-click on +/− being
          interpreted as a feature selection. */}
      <div
        data-ui="zoom-controls"
        onPointerDown={(e) => e.stopPropagation()}
        onPointerUp={(e) => e.stopPropagation()}
        onPointerMove={(e) => e.stopPropagation()}
        onClick={(e) => e.stopPropagation()}
        onDoubleClick={(e) => e.stopPropagation()}
        style={{
          position: "absolute",
          top: 6,
          right: 8,
          zIndex: 5,
          display: "flex",
          gap: 4,
          alignItems: "center",
          fontSize: 11,
          color: "#dbefe7",
          opacity: 0.95,
          background: "#46896c",
          padding: "4px 8px",
          borderRadius: 6,
          fontFamily: SANS_FONT_FAMILY,
        }}
      >
        <span>
          {fmtBp(visibleStartBp + 1)}–{fmtBp(visibleEndBp)} / {fmtBp(L)} bp
        </span>
        <button
          onClick={() => handleZoomBy(1 / 1.5)}
          style={{
            marginLeft: 4,
            padding: "1px 6px",
            background: "#46896c",
            color: "#dbefe7",
            border: "none",
            borderRadius: 4,
            cursor: "pointer",
            fontSize: 11,
          }}
          title="Zoom out"
        >
          −
        </button>
        <button
          onClick={() => handleZoomBy(1.5)}
          style={{
            padding: "1px 6px",
            background: "#46896c",
            color: "#dbefe7",
            border: "none",
            borderRadius: 4,
            cursor: "pointer",
            fontSize: 11,
          }}
          title="Zoom in"
        >
          +
        </button>
        <button
          onClick={handleResetView}
          style={{
            padding: "1px 6px",
            background: "#2d4a3e",
            color: "#dbefe7",
            border: "none",
            borderRadius: 4,
            cursor: "pointer",
            fontSize: 11,
          }}
          title="Reset zoom"
        >
          Fit
        </button>
      </div>

      <svg
        width="100%"
        height={svgHeight}
        style={{ display: "block", userSelect: "none" }}
      >
        {/* Alternating shading background for the visible track */}
        <rect x={0} y={0} width={containerWidth} height={svgHeight} fill={COLOR_TRACK} />

        {/* Selection highlight (renders behind features) */}
        {selVisible && (
          <rect
            x={bpToContentX(selA as number) - pan}
            y={HUD_HEIGHT + RULER_HEIGHT}
            width={Math.max(1, bpToContentX(selB as number) - bpToContentX(selA as number))}
            height={svgHeight - HUD_HEIGHT - RULER_HEIGHT}
            fill={COLOR_SELECTION}
            stroke={COLOR_SELECTION_BORDER}
            strokeDasharray="3 3"
          />
        )}

        {/* Ruler tick marks + labels */}
        <line
          x1={0}
          y1={HUD_HEIGHT + RULER_HEIGHT + BASE_HEIGHT - 1}
          x2={containerWidth}
          y2={HUD_HEIGHT + RULER_HEIGHT + BASE_HEIGHT - 1}
          stroke={COLOR_TICK}
        />
        {ticks.map((t) => {
          const x = bpToContentX(t) - pan;
          if (x < -40 || x > containerWidth + 40) return null;
          return (
            <g key={`tk-${t}`}>
              <line
                x1={x}
                y1={HUD_HEIGHT + RULER_HEIGHT + BASE_HEIGHT - 5}
                x2={x}
                y2={HUD_HEIGHT + RULER_HEIGHT + BASE_HEIGHT}
                stroke={COLOR_TICK}
              />
              <text
                x={x}
                y={HUD_HEIGHT + RULER_HEIGHT - 4}
                fontSize={10}
                textAnchor="middle"
                fill={COLOR_TICK_LABEL}
                fontFamily={SANS_FONT_FAMILY}
              >
                {fmtBp(t)}
              </text>
            </g>
          );
        })}

        {/* Sequence letters when zoomed in enough */}
        {letters}

        {/* Features by lane */}
        {renderLanes}
      </svg>
    </div>
  );
}
