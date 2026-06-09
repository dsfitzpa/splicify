"use client";

import React, { useCallback, useEffect, useLayoutEffect, useMemo, useRef, useState } from "react";
import type { SeqVizAnnotation } from "./InteractiveSequenceViewer";

// ─────────────────────────────────────────────────────────────────────────────
// Constants — fixed bp-per-pixel, Droid Sans Mono 14px
// ─────────────────────────────────────────────────────────────────────────────
// Sizes scaled +25 % (2026-04-29) to match the modular (circular) viewer's
// glyph weight: 14 px → 17.5 px nucleotides, etc.
const FONT_FAMILY = '"Droid Sans Mono", "DM Mono", "Courier New", monospace';
const BASE_FONT_SIZE = 17.5;
const BP_WIDTH = 14;            // px per base (fits 17.5 px glyph with air)
const TICK_LANE_HEIGHT = 18;
const STRAND_HEIGHT = 22;
const FEATURE_STRIP_HEIGHT = 15;
const FEATURE_STRIP_GAP = 3;
const CUT_TICK_HEIGHT = 8;
const LABEL_LANE_HEIGHT = 14; // lane below features for narrow-feature names + leader lines
const NARROW_LABEL_FONT = 11;
const ROW_V_GAP = 15;           // gap between rows
const GUTTER_W = 64;
const TICK_FONT_SIZE = 11;
const FEATURE_LABEL_FONT_SIZE = 12.5;
const ARROW_HEAD_PX = 8;        // tip-of-arrow horizontal length

const COLOR_TOP_STRAND = "#2d3135";
const COLOR_BOTTOM_STRAND = "#6b7075";
const COLOR_TICK = "#9aa0a6";
const COLOR_TICK_LABEL = "#72777e";
const COLOR_SELECTION = "#ffe066";
const COLOR_CUT_RE = "#C2185B";
const COLOR_CUT_GATEWAY = "#6A1B9A";

type Selection = { start: number; end: number } | null;

type Props = {
  sequence: string;
  annotations: SeqVizAnnotation[];
  height?: number;
  centerOnPosition?: number | null;
  topOnPosition?: number | null;
  selection?: Selection;
  onSelectionChange?: (sel: Selection) => void;
  onAnnotationClick?: (ann: SeqVizAnnotation, idx: number) => void;
};

type FeatureSegment = {
  ann: SeqVizAnnotation;
  idx: number;
  bpStart: number; // virtual strip coords
  bpEnd: number;   // inclusive
};
type PlacedSegment = FeatureSegment & { layer: number };

type CutGlyph = {
  ann: SeqVizAnnotation;
  topX: number;      // boundary x (left edge of cut_top base)
  bottomX: number;   // boundary x
  color: string;
  blunt: boolean;
  onRowTop: boolean; // cut_top position falls on this row
  onRowBottom: boolean; // cut_bottom falls on this row
};

type RenderedRow = {
  rowIdx: number;
  firstBp: number;
  lastBp: number;
  segments: PlacedSegment[];
  featureLayers: number;
  yOffset: number;
  height: number;
  rowBpStart: number; // virtual strip bp
  cuts: CutGlyph[];
  labelLaneHeight: number;
  translationGlyphs: TranslationGlyph[];
  aaLaneHeight: number;
};

const COMPLEMENT: Record<string, string> = {
  A: "T", T: "A", G: "C", C: "G", N: "N",
  a: "t", t: "a", g: "c", c: "g", n: "n",
};

function complement(base: string): string {
  return COMPLEMENT[base] ?? base;
}

// ─────────────────────────────────────────────────────────────────────────────
// Translation annotation rendering — inline AA strip with per-AA hit boxes
// ─────────────────────────────────────────────────────────────────────────────
const AA_LANE_GAP = 4;
const AA_LANE_HEIGHT = 20;
const AA_FONT_SIZE = 12.5;

const AA_FULL_NAME: Record<string, string> = {
  A: "Alanine", R: "Arginine", N: "Asparagine", D: "Aspartic acid",
  C: "Cysteine", E: "Glutamic acid", Q: "Glutamine", G: "Glycine",
  H: "Histidine", I: "Isoleucine", L: "Leucine", K: "Lysine",
  M: "Methionine", F: "Phenylalanine", P: "Proline", S: "Serine",
  T: "Threonine", W: "Tryptophan", Y: "Tyrosine", V: "Valine",
  "*": "Stop", X: "Unknown",
};

type FeatureRegion = { name: string; aa_start: number; aa_end: number; feature_type?: string };

type TranslationGlyph = {
  ann: SeqVizAnnotation;
  annIdx: number;
  aaIdx: number;     // 1-based AA position within the ORF
  letter: string;
  x: number;         // pixel x within the row svg
  ntStart: number;   // 0-based bp index where this AA starts
};

type AaTooltipState = {
  x: number;
  y: number;
  orfName: string;
  orfLen: number;
  aaIdx: number;
  letter: string;
  region: FeatureRegion | null;
};

function isTranslationAnn(a: SeqVizAnnotation): boolean {
  return (a as any).layer === "translation";
}

// Palette + selection logic mirrors CircularPlasmidViewer.getAnnotationColor
// so the same annotation shows the same color across the two viewers.
const FEATURE_ALT_PALETTE = ["#4a7c59", "#6b5b95", "#3b82f6", "#f97316", "#22c55e", "#ef4444"];
const DESIGN_PALETTE = [
  "#ec4899", "#f59e0b", "#06b6d4", "#8b5cf6",
  "#22c55e", "#ef4444", "#14b8a6", "#f97316",
];
const PRESERVED_GRAY_PALETTE = ["#9ca3af", "#6b7280", "#cbd5e1"];

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

function featureSpan(ann: SeqVizAnnotation, totalLength: number): { s: number; e: number } {
  if (ann.start <= ann.end) return { s: ann.start, e: ann.end - 1 };
  return { s: ann.start, e: ann.end + totalLength - 1 };
}

// Packs segments into layers using greedy first-fit.
function packLayers(segs: FeatureSegment[]): PlacedSegment[] {
  const sorted = [...segs].sort((a, b) => a.bpStart - b.bpStart || a.bpEnd - b.bpEnd);
  const layerEnds: number[] = [];
  const placed: PlacedSegment[] = [];
  for (const seg of sorted) {
    let layer = -1;
    for (let li = 0; li < layerEnds.length; li++) {
      if (layerEnds[li] < seg.bpStart) {
        layer = li;
        layerEnds[li] = seg.bpEnd;
        break;
      }
    }
    if (layer < 0) {
      layer = layerEnds.length;
      layerEnds.push(seg.bpEnd);
    }
    placed.push({ ...seg, layer });
  }
  return placed;
}

export default function LinearSequenceViewer({
  sequence,
  annotations,
  height = 500,
  centerOnPosition = null,
  topOnPosition = null,
  selection = null,
  onSelectionChange,
  onAnnotationClick,
}: Props) {
  const totalLength = sequence.length;
  const containerRef = useRef<HTMLDivElement | null>(null);
  const [containerWidth, setContainerWidth] = useState(600);

  // If any annotation carries added_by_design, this is a design-output .gb;
  // gray out the preserved annotations so the design overlay reads at a glance.
  const hasDesignAdditions = useMemo(
    () => annotations.some((a) => (a as any).added_by_design === true),
    [annotations],
  );

  useEffect(() => {
    const el = containerRef.current;
    if (!el) return;
    const ro = new ResizeObserver((entries) => {
      for (const e of entries) setContainerWidth(e.contentRect.width);
    });
    ro.observe(el);
    setContainerWidth(el.clientWidth);
    return () => ro.disconnect();
  }, []);

  const basesPerRow = useMemo(() => {
    const usable = Math.max(60, containerWidth - GUTTER_W * 2 - 16);
    return Math.max(10, Math.floor(usable / BP_WIDTH));
  }, [containerWidth]);

  // Native-scroll model (matches SeqViz): the strip is laid out at its full
  // natural height inside an overflow-y: auto scroller. The browser handles
  // every wheel / trackpad / scrollbar interaction natively. We mirror
  // scrollerRef.current.scrollTop into React state (rAF-throttled) so the
  // virtualisation memo can pick the visible row window.
  const scrollerRef = useRef<HTMLDivElement | null>(null);
  const [scrollTop, setScrollTop] = useState<number>(0);

  // Viewport (minus 8px top + 8px bottom padding); used by everything below.
  const availableHeight = Math.max(40, height - 16);

  // External re-center / re-top requests from the find panel / circular
  // viewer. Both flow through useLayoutEffects below whose deps include
  // the trigger value so each new request actually scrolls (the old
  // pendingRef pattern only re-fired the resolver when allRows changed,
  // which never happened for plain selection / search updates).
  const lastCenterRef = useRef<number | null>(null);
  const lastTopRef = useRef<number | null>(null);

  const realBp = useCallback(
    (offset: number) => Math.max(0, Math.min(totalLength - 1, offset)),
    [totalLength]
  );

  // Pre-compute every row in the sequence (height + segments + cuts) once
  // per (annotations, basesPerRow) change, with cumulative yOffsets. The
  // strip renders all rows positioned at their absolute yOffset; only the
  // ones in the visible scroll window are actually mounted.
  const allRows = useMemo<RenderedRow[]>(() => {
    if (totalLength === 0 || basesPerRow === 0) return [];

    // Translation annotations get their own per-AA strip lane below the
    // bottom strand — do not pack them into the regular feature lanes.
    const allSegs: FeatureSegment[] = [];
    const translationAnns: { ann: SeqVizAnnotation; idx: number }[] = [];
    annotations.forEach((ann, idx) => {
      if (isTranslationAnn(ann)) {
        translationAnns.push({ ann, idx });
        return;
      }
      const { s, e } = featureSpan(ann, totalLength);
      // Origin-spanning features (start > end on circular plasmids) split
      // into [s..totalLength-1] and [0..e-totalLength] in linear projection.
      if (e < totalLength) {
        allSegs.push({ ann, idx, bpStart: s, bpEnd: e });
      } else {
        allSegs.push({ ann, idx, bpStart: s, bpEnd: totalLength - 1 });
        allSegs.push({ ann, idx, bpStart: 0, bpEnd: e - totalLength });
      }
    });
    allSegs.sort((a, b) => a.bpStart - b.bpStart);

    const rowCount = Math.ceil(totalLength / basesPerRow);
    const out: RenderedRow[] = [];
    let y = 0;
    let segCursor = 0;

    for (let rowIdx = 0; rowIdx < rowCount; rowIdx++) {
      const rowBpStart = rowIdx * basesPerRow;
      const rowBpEnd = Math.min(totalLength - 1, rowBpStart + basesPerRow - 1);

      while (segCursor < allSegs.length && allSegs[segCursor].bpEnd < rowBpStart) {
        segCursor++;
      }
      const rowSegs: FeatureSegment[] = [];
      for (let i = segCursor; i < allSegs.length; i++) {
        const seg = allSegs[i];
        if (seg.bpStart > rowBpEnd) break;
        rowSegs.push({
          ann: seg.ann,
          idx: seg.idx,
          bpStart: Math.max(seg.bpStart, rowBpStart),
          bpEnd: Math.min(seg.bpEnd, rowBpEnd),
        });
      }
      const placed = packLayers(rowSegs);
      const featureLayers = placed.reduce((m, s) => Math.max(m, s.layer + 1), 0);

      const cuts: CutGlyph[] = [];
      annotations.forEach((ann) => {
        const cp = (ann as any).cut_profile;
        if (!cp || typeof cp.cut_top !== "number" || typeof cp.cut_bottom !== "number") return;
        if ((ann as any).feature_family === "primer_design_warning") return;
        const vTop = cp.cut_top;
        const vBot = cp.cut_bottom;
        const onRowTop = vTop >= rowBpStart && vTop <= rowBpEnd + 1;
        const onRowBottom = vBot >= rowBpStart && vBot <= rowBpEnd + 1;
        if (!onRowTop && !onRowBottom) return;
        const topX = (vTop - rowBpStart) * BP_WIDTH;
        const bottomX = (vBot - rowBpStart) * BP_WIDTH;
        const color = (ann as any).feature_family === "gateway_att" ? COLOR_CUT_GATEWAY : COLOR_CUT_RE;
        cuts.push({
          ann,
          topX,
          bottomX,
          color,
          blunt: cp.cut_top === cp.cut_bottom,
          onRowTop,
          onRowBottom,
        });
      });

      const hasNarrowFeatures = placed.some((seg) => {
        const w = (seg.bpEnd - seg.bpStart + 1) * BP_WIDTH;
        return seg.ann.name && w <= 36;
      });
      const labelLaneHeight = hasNarrowFeatures ? LABEL_LANE_HEIGHT : 0;

      // Per-row AA glyphs: for each translation annotation overlapping this
      // row, emit one glyph per amino acid whose first codon nt falls into
      // the row's bp window. nt → AA position uses (i)*3 from the ORF start
      // (forward) or (orfEnd - 3*(i+1)) (reverse).
      const translationGlyphs: TranslationGlyph[] = [];
      translationAnns.forEach(({ ann, idx: annIdx }) => {
        const dir = (ann.direction ?? ann.strand ?? 1) === -1 ? -1 : 1;
        const orfStart = ann.start;
        const orfEnd = ann.end;
        const aaSeq = (ann.metadata?.aa_sequence as string) || "";
        if (!aaSeq) return;
        for (let i = 0; i < aaSeq.length; i++) {
          const ntStart = dir === 1 ? orfStart + 3 * i : orfEnd - 3 * (i + 1);
          if (ntStart < rowBpStart || ntStart > rowBpEnd) continue;
          translationGlyphs.push({
            ann,
            annIdx,
            aaIdx: i + 1,
            letter: aaSeq[i],
            x: (ntStart - rowBpStart) * BP_WIDTH,
            ntStart,
          });
        }
      });
      const aaLaneHeight = translationGlyphs.length > 0
        ? AA_LANE_HEIGHT + AA_LANE_GAP
        : 0;

      const rowHeight =
        TICK_LANE_HEIGHT +
        STRAND_HEIGHT +
        Math.max(0, featureLayers * (FEATURE_STRIP_HEIGHT + FEATURE_STRIP_GAP)) +
        labelLaneHeight +
        STRAND_HEIGHT +
        aaLaneHeight +
        CUT_TICK_HEIGHT +
        ROW_V_GAP;

      out.push({
        rowIdx,
        firstBp: rowBpStart,
        lastBp: rowBpEnd,
        segments: placed,
        featureLayers,
        yOffset: y,
        height: rowHeight,
        rowBpStart,
        cuts,
        labelLaneHeight,
        translationGlyphs,
        aaLaneHeight,
      });
      y += rowHeight;
    }
    return out;
  }, [annotations, basesPerRow, totalLength]);

  const totalContentHeight = useMemo(
    () => allRows.reduce((s, r) => s + r.height, 0),
    [allRows]
  );

  // Visible row range — binary search on cumulative yOffset, with a 2-row
  // overscan above and below so a fast scroll never shows empty gaps.
  const visibleRange = useMemo(() => {
    if (allRows.length === 0) return { first: 0, last: -1 };
    const top = scrollTop;
    const bottom = scrollTop + availableHeight;
    let lo = 0;
    let hi = allRows.length - 1;
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if (allRows[mid].yOffset + allRows[mid].height > top) hi = mid;
      else lo = mid + 1;
    }
    const first = lo;
    let lo2 = first;
    let hi2 = allRows.length - 1;
    while (lo2 < hi2) {
      const mid = (lo2 + hi2 + 1) >> 1;
      if (allRows[mid].yOffset < bottom) lo2 = mid;
      else hi2 = mid - 1;
    }
    const last = lo2;
    return {
      first: Math.max(0, first - 2),
      last: Math.min(allRows.length - 1, last + 2),
    };
  }, [allRows, scrollTop, availableHeight]);

  const rows = useMemo<RenderedRow[]>(
    () => allRows.slice(visibleRange.first, visibleRange.last + 1),
    [allRows, visibleRange.first, visibleRange.last]
  );

  const visibleBases = rows.length * basesPerRow;

  // Resolve centerOnPosition into a scrollTop once allRows is available.
  // Re-fires whenever centerOnPosition changes OR allRows recomputes —
  // either condition can be the one that makes the scroll possible.
  useLayoutEffect(() => {
    if (centerOnPosition == null || allRows.length === 0 || totalLength === 0) return;
    lastCenterRef.current = centerOnPosition;
    const targetIdx = Math.min(
      allRows.length - 1,
      Math.max(0, Math.floor(centerOnPosition / basesPerRow))
    );
    const r = allRows[targetIdx];
    const desired = r.yOffset + r.height / 2 - availableHeight / 2;
    const max = Math.max(0, totalContentHeight - availableHeight);
    const clamped = Math.max(0, Math.min(max, desired));
    if (scrollerRef.current) scrollerRef.current.scrollTop = clamped;
    setScrollTop(clamped);
  }, [centerOnPosition, allRows, basesPerRow, availableHeight, totalContentHeight, totalLength]);

  // Top-aligned counterpart: row containing topOnPosition becomes the
  // first visible row. Used by the find panel + circular-viewer
  // selections so the start of the match / feature lands at the top.
  useLayoutEffect(() => {
    if (topOnPosition == null || allRows.length === 0 || totalLength === 0) return;
    if (lastTopRef.current === topOnPosition) return;
    lastTopRef.current = topOnPosition;
    const targetIdx = Math.min(
      allRows.length - 1,
      Math.max(0, Math.floor(topOnPosition / basesPerRow))
    );
    const r = allRows[targetIdx];
    const desired = r.yOffset;
    const max = Math.max(0, totalContentHeight - availableHeight);
    const clamped = Math.max(0, Math.min(max, desired));
    if (scrollerRef.current) scrollerRef.current.scrollTop = clamped;
    setScrollTop(clamped);
  }, [topOnPosition, allRows, basesPerRow, availableHeight, totalContentHeight, totalLength]);

  // Reset to the top of the strip when the sequence (totalLength) flips.
  useEffect(() => {
    if (scrollerRef.current) scrollerRef.current.scrollTop = 0;
    setScrollTop(0);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [totalLength]);
  // `firstVisibleBp` / `lastVisibleBp` were only used by the 5'/3'
  // banners; both removed 2026-04-29.

  // ───────────────────────────────────────────────────────────────────────────
  // Mouse → virtual bp offset
  // ───────────────────────────────────────────────────────────────────────────
  const hitTest = useCallback(
    (clientX: number, clientY: number): number | null => {
      const scroller = scrollerRef.current;
      if (!scroller || allRows.length === 0) return null;
      const rect = scroller.getBoundingClientRect();
      const y = Math.max(0, clientY - rect.top + scroller.scrollTop);
      let lo = 0;
      let hi = allRows.length - 1;
      while (lo < hi) {
        const mid = (lo + hi) >> 1;
        if (allRows[mid].yOffset + allRows[mid].height > y) hi = mid;
        else lo = mid + 1;
      }
      const row = allRows[lo];
      const xInRow = clientX - rect.left - GUTTER_W;
      const col = Math.max(0, Math.min(basesPerRow - 1, Math.floor(xInRow / BP_WIDTH)));
      return Math.max(0, Math.min(totalLength - 1, row.rowIdx * basesPerRow + col));
    },
    [allRows, basesPerRow, totalLength]
  );

  // ───────────────────────────────────────────────────────────────────────────
  // Native scroll mirror — onScroll → setScrollTop, throttled by rAF
  // ───────────────────────────────────────────────────────────────────────────
  const scrollRafRef = useRef<number | null>(null);
  const onScroll = useCallback(() => {
    if (scrollRafRef.current != null) return;
    scrollRafRef.current = requestAnimationFrame(() => {
      scrollRafRef.current = null;
      const el = scrollerRef.current;
      if (el) setScrollTop(el.scrollTop);
    });
  }, []);
  useEffect(() => {
    return () => {
      if (scrollRafRef.current != null) {
        cancelAnimationFrame(scrollRafRef.current);
        scrollRafRef.current = null;
      }
    };
  }, []);

  // ───────────────────────────────────────────────────────────────────────────
  // Drag selection + click-to-select-feature + double-click
  // ───────────────────────────────────────────────────────────────────────────
  const dragStartOffsetRef = useRef<number | null>(null);
  const [dragRange, setDragRange] = useState<{ a: number; b: number } | null>(null);
  const clickTimer = useRef<number | null>(null);

  const emitSelection = useCallback(
    (aOff: number, bOff: number) => {
      const lo = Math.min(aOff, bOff);
      const hi = Math.max(aOff, bOff);
      const realStart = realBp(lo);
      const realEndInclusive = realBp(hi);
      // Half-open end
      const end = (realEndInclusive + 1) % totalLength;
      if (lo === hi) {
        onSelectionChange?.(null);
      } else {
        onSelectionChange?.({ start: realStart, end });
      }
    },
    [realBp, totalLength, onSelectionChange]
  );

  // In native-scroll mode, "offset" is the absolute base index — same as
  // the real bp.
  const realBpToOffset = useCallback(
    (realPos: number): number => Math.max(0, Math.min(totalLength - 1, realPos)),
    [totalLength]
  );

  const onMouseDownBg = useCallback(
    (e: React.MouseEvent) => {
      if ((e.target as Element).closest("[data-annotation]")) return;
      if ((e.target as Element).closest("[data-cut]")) return;
      const off = hitTest(e.clientX, e.clientY);
      if (off == null) return;
      // Shift-click: extend the existing selection. We preserve BOTH of
      // the existing endpoints so a click before the selection extends
      // backward without dropping the far end (and vice-versa).
      if (e.shiftKey && selection && selection.start !== selection.end) {
        const exStart = selection.start;
        // selection.end is half-open in this codebase, so the far inclusive
        // endpoint is end-1 (mod totalLength).
        const exEndIncl = (selection.end - 1 + totalLength) % totalLength;
        const clickedReal = realBp(off);
        const lo = Math.min(exStart, exEndIncl, clickedReal);
        const hi = Math.max(exStart, exEndIncl, clickedReal);
        onSelectionChange?.({ start: lo, end: (hi + 1) % totalLength });
        // Begin a drag from the anchor opposite to the click so further
        // mouse-move keeps extending in the direction the user is dragging.
        const anchorReal = clickedReal <= exStart ? exEndIncl : exStart;
        const anchorOff = realBpToOffset(anchorReal);
        dragStartOffsetRef.current = anchorOff;
        setDragRange({ a: anchorOff, b: off });
        e.preventDefault();
        return;
      }
      dragStartOffsetRef.current = off;
      setDragRange({ a: off, b: off });
      e.preventDefault();
    },
    [hitTest, selection, realBp, totalLength, onSelectionChange, realBpToOffset]
  );

  // Live cursor pos for the auto-scroll RAF loop.
  const lastMouseRef = useRef<{ x: number; y: number } | null>(null);
  const autoScrollRafRef = useRef<number | null>(null);

  useEffect(() => {
    if (dragStartOffsetRef.current == null) return;

    const tickAutoScroll = () => {
      const el = containerRef.current;
      const pos = lastMouseRef.current;
      if (!el || !pos) {
        autoScrollRafRef.current = requestAnimationFrame(tickAutoScroll);
        return;
      }
      const rect = el.getBoundingClientRect();
      const EDGE = 24; // px hot-zone above/below the viewer
      let direction = 0;
      if (pos.y < rect.top + EDGE) direction = -1;
      else if (pos.y > rect.bottom - EDGE) direction = 1;
      if (direction !== 0 && scrollerRef.current) {
        const overshoot =
          direction < 0 ? rect.top + EDGE - pos.y : pos.y - (rect.bottom - EDGE);
        // ~30 px / frame at the deepest edge; matches SeqViz's drag pull.
        const pxPerFrame = Math.max(4, Math.min(40, Math.ceil(overshoot * 1.2)));
        scrollerRef.current.scrollTop += direction * pxPerFrame;
        const off = hitTest(pos.x, pos.y);
        if (off != null) setDragRange((prev) => (prev ? { a: prev.a, b: off } : null));
      }
      autoScrollRafRef.current = requestAnimationFrame(tickAutoScroll);
    };
    autoScrollRafRef.current = requestAnimationFrame(tickAutoScroll);

    const move = (ev: MouseEvent) => {
      lastMouseRef.current = { x: ev.clientX, y: ev.clientY };
      const off = hitTest(ev.clientX, ev.clientY);
      if (off == null) return;
      setDragRange((prev) => (prev ? { a: prev.a, b: off } : null));
    };
    const up = (ev: MouseEvent) => {
      const start = dragStartOffsetRef.current;
      dragStartOffsetRef.current = null;
      lastMouseRef.current = null;
      if (autoScrollRafRef.current != null) {
        cancelAnimationFrame(autoScrollRafRef.current);
        autoScrollRafRef.current = null;
      }
      if (start == null) return;
      const off = hitTest(ev.clientX, ev.clientY) ?? start;
      setDragRange(null);
      emitSelection(start, off);
    };
    window.addEventListener("mousemove", move);
    window.addEventListener("mouseup", up);
    return () => {
      window.removeEventListener("mousemove", move);
      window.removeEventListener("mouseup", up);
      if (autoScrollRafRef.current != null) {
        cancelAnimationFrame(autoScrollRafRef.current);
        autoScrollRafRef.current = null;
      }
    };
  }, [dragRange, hitTest, emitSelection, totalLength, basesPerRow]);

  const onFeatureClick = useCallback(
    (seg: PlacedSegment, ev?: React.MouseEvent) => {
      // Shift-click on a feature: preserve BOTH existing endpoints so the
      // selection grows in either direction depending on where the
      // feature sits relative to the prior selection.
      if (ev?.shiftKey && selection && selection.start !== selection.end) {
        const exStart = selection.start;
        const exEndIncl = (selection.end - 1 + totalLength) % totalLength;
        const fStart = seg.ann.start;
        const fEnd = seg.ann.end;
        const points = [exStart, exEndIncl, fStart, fEnd];
        const lo = Math.min(...points);
        const hi = Math.max(...points);
        onSelectionChange?.({ start: lo, end: hi });
        return;
      }
      if (clickTimer.current != null) {
        window.clearTimeout(clickTimer.current);
        clickTimer.current = null;
        onAnnotationClick?.(seg.ann, seg.idx);
        return;
      }
      clickTimer.current = window.setTimeout(() => {
        clickTimer.current = null;
        onSelectionChange?.({ start: seg.ann.start, end: seg.ann.end });
      }, 240);
    },
    [onAnnotationClick, onSelectionChange, selection]
  );

  // In-selection check (half-open, handles origin wrap)
  const inSelection = useCallback(
    (bp: number): boolean => {
      if (!selection) return false;
      const { start, end } = selection;
      if (start === end) return false;
      if (start < end) return bp >= start && bp < end;
      return bp >= start || bp < end;
    },
    [selection]
  );

  const inDragRange = useCallback(
    (offset: number): boolean => {
      if (!dragRange) return false;
      const lo = Math.min(dragRange.a, dragRange.b);
      const hi = Math.max(dragRange.a, dragRange.b);
      return offset >= lo && offset <= hi;
    },
    [dragRange]
  );

  const rowSvgWidth = Math.max(40, containerWidth - GUTTER_W * 2);

  const [aaTooltip, setAaTooltip] = useState<AaTooltipState | null>(null);

  // Close tooltip on outside click / ESC.
  useEffect(() => {
    if (!aaTooltip) return;
    const onDoc = (ev: MouseEvent) => {
      const tgt = ev.target as HTMLElement | null;
      if (tgt && tgt.closest("[data-aa-glyph], [data-aa-tooltip]")) return;
      setAaTooltip(null);
    };
    const onKey = (ev: KeyboardEvent) => {
      if (ev.key === "Escape") setAaTooltip(null);
    };
    window.addEventListener("mousedown", onDoc);
    window.addEventListener("keydown", onKey);
    return () => {
      window.removeEventListener("mousedown", onDoc);
      window.removeEventListener("keydown", onKey);
    };
  }, [aaTooltip]);

  const onAaClick = useCallback(
    (g: TranslationGlyph, ev: React.MouseEvent) => {
      ev.stopPropagation();
      const meta = (g.ann.metadata || {}) as { feature_regions?: FeatureRegion[]; aa_length?: number };
      const regions = meta.feature_regions || [];
      const region = regions.find((r) => g.aaIdx >= r.aa_start && g.aaIdx <= r.aa_end) || null;
      // Tooltip is positioned in viewer-local coordinates anchored to the
      // mouse click. The wrapper has position:relative so we use clientX/Y
      // minus the container's offset.
      const containerEl = containerRef.current;
      const rect = containerEl ? containerEl.getBoundingClientRect() : null;
      const x = rect ? ev.clientX - rect.left + 6 : ev.clientX;
      const y = rect ? ev.clientY - rect.top + 6 : ev.clientY;
      setAaTooltip({
        x,
        y,
        orfName: g.ann.name || "Translation",
        orfLen: meta.aa_length || 0,
        aaIdx: g.aaIdx,
        letter: g.letter,
        region,
      });
    },
    []
  );

  return (
    <div
      ref={containerRef}
      onMouseDown={onMouseDownBg}
      style={{
        height,
        width: "100%",
        backgroundColor: "#ffffff",
        padding: 8,
        boxSizing: "border-box",
        overflow: "hidden",
        fontFamily: FONT_FAMILY,
        position: "relative",
        userSelect: "none",
        cursor: "text",
      }}
    >

      <div
        ref={scrollerRef}
        onScroll={onScroll}
        style={{
          position: "relative",
          width: "100%",
          height: availableHeight,
          overflowY: "auto",
          overflowX: "hidden",
        }}
      >
        <div style={{ position: "relative", width: "100%", height: totalContentHeight }}>
        {rows.map((row) => {
          const tickY = TICK_LANE_HEIGHT - 2;
          const topStrandY = TICK_LANE_HEIGHT + STRAND_HEIGHT / 2;
          const featureLaneY = TICK_LANE_HEIGHT + STRAND_HEIGHT;
          const featureLaneHeight = row.featureLayers * (FEATURE_STRIP_HEIGHT + FEATURE_STRIP_GAP);
          const labelLaneY = featureLaneY + featureLaneHeight; // lane below features
          const bottomStrandY =
            TICK_LANE_HEIGHT + STRAND_HEIGHT + featureLaneHeight + row.labelLaneHeight + STRAND_HEIGHT / 2;
          const aaLaneY =
            TICK_LANE_HEIGHT + STRAND_HEIGHT + featureLaneHeight + row.labelLaneHeight + STRAND_HEIGHT + AA_LANE_GAP;
          const cutLaneY = TICK_LANE_HEIGHT + STRAND_HEIGHT + featureLaneHeight + row.labelLaneHeight + STRAND_HEIGHT + row.aaLaneHeight;

          return (
            <div
              key={row.rowIdx}
              style={{
                position: "absolute",
                top: row.yOffset,
                left: 0,
                right: 0,
                height: row.height,
                display: "flex",
                alignItems: "flex-start",
              }}
            >
              {/* Left gutter — kept for layout symmetry; position label
                  removed (2026-04-29) to keep only the inline 10-bp ticks. */}
              <div style={{ width: GUTTER_W }} />

              <svg width={rowSvgWidth} height={row.height} style={{ display: "block" }}>
                {/* Selection + live drag highlight */}
                {Array.from({ length: basesPerRow }).map((_, col) => {
                  const offset = row.rowIdx * basesPerRow + col;
                  if (offset >= totalLength) return null;
                  const bp = offset;
                  const hit = inSelection(bp) || inDragRange(offset);
                  if (!hit) return null;
                  return (
                    <rect
                      key={`sel-${col}`}
                      x={col * BP_WIDTH}
                      y={topStrandY - STRAND_HEIGHT / 2}
                      width={BP_WIDTH}
                      height={
                        STRAND_HEIGHT + featureLaneHeight + row.labelLaneHeight + STRAND_HEIGHT
                      }
                      fill={COLOR_SELECTION}
                      opacity={0.55}
                    />
                  );
                })}

                {/* Tick marks every 10 bases */}
                {Array.from({ length: basesPerRow }).map((_, col) => {
                  const bp = row.rowIdx * basesPerRow + col;
                  if (bp >= totalLength) return null;
                  const oneIndexed = bp + 1;
                  if (oneIndexed % 10 !== 0) return null;
                  return (
                    <g key={`tick-${col}`}>
                      <line
                        x1={col * BP_WIDTH + BP_WIDTH / 2}
                        y1={tickY - 4}
                        x2={col * BP_WIDTH + BP_WIDTH / 2}
                        y2={tickY + 3}
                        stroke={COLOR_TICK}
                        strokeWidth={1.2}
                      />
                      <text
                        x={col * BP_WIDTH + BP_WIDTH / 2}
                        y={tickY - 6}
                        fontSize={TICK_FONT_SIZE}
                        fontFamily={FONT_FAMILY}
                        fill={COLOR_TICK_LABEL}
                        textAnchor="middle"
                      >
                        {oneIndexed}
                      </text>
                    </g>
                  );
                })}

                {/* Top strand */}
                {Array.from({ length: basesPerRow }).map((_, col) => {
                  const bp = row.rowIdx * basesPerRow + col;
                  if (bp >= totalLength) return null;
                  const base = sequence[bp] ?? "";
                  return (
                    <text
                      key={`top-${col}`}
                      x={col * BP_WIDTH + BP_WIDTH / 2}
                      y={topStrandY}
                      fontSize={BASE_FONT_SIZE}
                      fontFamily={FONT_FAMILY}
                      fill={COLOR_TOP_STRAND}
                      textAnchor="middle"
                      dominantBaseline="middle"
                    >
                      {base}
                    </text>
                  );
                })}

                {/* Feature segments — arrow-shaped to indicate direction
                    (matches the modular/circular viewer). The forward arrow
                    has a chevron tip at the right end at the segment's last
                    base; reverse points left. Segments that are clipped at
                    the row edge keep a flat (clipped) end and only render
                    the chevron when the actual feature endpoint is visible
                    within this row. */}
                {row.segments.map((seg, i) => {
                  const s = seg.bpStart;
                  const e = seg.bpEnd;
                  const x = (s - row.rowBpStart) * BP_WIDTH;
                  const w = (e - s + 1) * BP_WIDTH;
                  const y = featureLaneY + seg.layer * (FEATURE_STRIP_HEIGHT + FEATURE_STRIP_GAP);
                  const color = featureColor(seg.ann, seg.idx, hasDesignAdditions);
                  const dir = (seg.ann.direction ?? seg.ann.strand ?? 0);
                  const fullSpan = featureSpan(seg.ann, totalLength);
                  // Tip of the arrow corresponds to the feature's directional
                  // endpoint. Only show a chevron when that endpoint actually
                  // sits inside this row (i.e. not clipped at the row edge).
                  const featStart = fullSpan.s;
                  const featEnd = fullSpan.e;
                  const showRightChevron = dir === 1 && (
                    e === featEnd ||
                    e === featEnd + totalLength ||
                    e === featEnd - totalLength
                  );
                  const showLeftChevron = dir === -1 && (
                    s === featStart ||
                    s === featStart + totalLength ||
                    s === featStart - totalLength
                  );
                  const headW = Math.min(ARROW_HEAD_PX, Math.max(0, w * 0.45));
                  const yMid = y + FEATURE_STRIP_HEIGHT / 2;
                  let path = "";
                  if (showRightChevron && headW > 0) {
                    const shaftRight = x + w - headW;
                    path = `M ${x} ${y}
                            L ${shaftRight} ${y}
                            L ${x + w} ${yMid}
                            L ${shaftRight} ${y + FEATURE_STRIP_HEIGHT}
                            L ${x} ${y + FEATURE_STRIP_HEIGHT} Z`;
                  } else if (showLeftChevron && headW > 0) {
                    const shaftLeft = x + headW;
                    path = `M ${x + w} ${y}
                            L ${shaftLeft} ${y}
                            L ${x} ${yMid}
                            L ${shaftLeft} ${y + FEATURE_STRIP_HEIGHT}
                            L ${x + w} ${y + FEATURE_STRIP_HEIGHT} Z`;
                  }
                  return (
                    <g
                      key={`seg-${row.rowIdx}-${i}`}
                      data-annotation
                      style={{ cursor: "pointer" }}
                      onMouseDown={(e) => e.stopPropagation()}
                      onClick={(e) => {
                        e.stopPropagation();
                        onFeatureClick(seg, e);
                      }}
                    >
                      {path ? (
                        <path d={path} fill={color} opacity={0.8} />
                      ) : (
                        <rect
                          x={x}
                          y={y}
                          width={w}
                          height={FEATURE_STRIP_HEIGHT}
                          fill={color}
                          opacity={0.8}
                          rx={2}
                        />
                      )}
                      {w > 40 ? (
                        <text
                          x={x + w / 2}
                          y={yMid + 1}
                          fontSize={FEATURE_LABEL_FONT_SIZE}
                          fontFamily={FONT_FAMILY}
                          fill="#ffffff"
                          textAnchor="middle"
                          dominantBaseline="middle"
                          pointerEvents="none"
                        >
                          {seg.ann.name.length * 7 > w - 6
                            ? seg.ann.name.slice(0, Math.max(1, Math.floor((w - 6) / 7)))
                            : seg.ann.name}
                        </text>
                      ) : (
                        seg.ann.name && (
                          <g pointerEvents="none">
                            <line
                              x1={x + w / 2}
                              y1={y + FEATURE_STRIP_HEIGHT}
                              x2={x + w / 2}
                              y2={labelLaneY + LABEL_LANE_HEIGHT - 2}
                              stroke={color}
                              strokeWidth={0.7}
                              opacity={0.65}
                            />
                            <text
                              x={x + w / 2}
                              y={labelLaneY + LABEL_LANE_HEIGHT - 1}
                              fontSize={NARROW_LABEL_FONT}
                              fontFamily={FONT_FAMILY}
                              fill={color}
                              textAnchor="middle"
                              dominantBaseline="alphabetic"
                            >
                              {seg.ann.name.length > 18 ? seg.ann.name.slice(0, 17) + "…" : seg.ann.name}
                            </text>
                          </g>
                        )
                      )}
                    </g>
                  );
                })}

                {/* Bottom strand — complement */}
                {Array.from({ length: basesPerRow }).map((_, col) => {
                  const bp = row.rowIdx * basesPerRow + col;
                  if (bp >= totalLength) return null;
                  const base = complement(sequence[bp] ?? "");
                  return (
                    <text
                      key={`bot-${col}`}
                      x={col * BP_WIDTH + BP_WIDTH / 2}
                      y={bottomStrandY}
                      fontSize={BASE_FONT_SIZE}
                      fontFamily={FONT_FAMILY}
                      fill={COLOR_BOTTOM_STRAND}
                      textAnchor="middle"
                      dominantBaseline="middle"
                    >
                      {base}
                    </text>
                  );
                })}

                {/* Translation AA strip — one clickable glyph per residue,
                    centered over the residue's first nucleotide (3 bp wide). */}
                {row.translationGlyphs.map((g, i) => {
                  const glyphX = g.x;
                  const glyphW = 3 * BP_WIDTH;
                  const dir = (g.ann.direction ?? g.ann.strand ?? 1) === -1 ? -1 : 1;
                  const fill = dir === -1 ? "#5e35b1" : "#673AB7";
                  const isHovered = !!aaTooltip
                    && aaTooltip.orfName === (g.ann.name || "Translation")
                    && aaTooltip.aaIdx === g.aaIdx;
                  return (
                    <g
                      key={`aa-${row.rowIdx}-${i}`}
                      data-aa-glyph
                      style={{ cursor: "pointer" }}
                      onMouseDown={(e) => e.stopPropagation()}
                      onClick={(e) => onAaClick(g, e)}
                    >
                      <rect
                        x={glyphX}
                        y={aaLaneY}
                        width={glyphW}
                        height={AA_LANE_HEIGHT}
                        fill={isHovered ? fill : "transparent"}
                        opacity={isHovered ? 0.18 : 1}
                        rx={2}
                      />
                      <rect
                        x={glyphX + 0.5}
                        y={aaLaneY + 0.5}
                        width={glyphW - 1}
                        height={AA_LANE_HEIGHT - 1}
                        fill="none"
                        stroke={fill}
                        strokeOpacity={0.25}
                        strokeWidth={0.7}
                        rx={2}
                      />
                      <text
                        x={glyphX + glyphW / 2}
                        y={aaLaneY + AA_LANE_HEIGHT / 2 + 1}
                        fontSize={AA_FONT_SIZE}
                        fontFamily={FONT_FAMILY}
                        fill={fill}
                        textAnchor="middle"
                        dominantBaseline="middle"
                        pointerEvents="none"
                      >
                        {g.letter}
                      </text>
                    </g>
                  );
                })}

                {/* Cut-profile glyphs — staggered bracket. Forward tick on top
                    strand edge, reverse tick on bottom strand edge, horizontal
                    connector between them at the mid-strand y. */}
                {row.cuts.map((cut, i) => {
                  const topY1 = topStrandY - STRAND_HEIGHT / 2;
                  const topY2 = topY1 - CUT_TICK_HEIGHT;
                  const botY1 = bottomStrandY + STRAND_HEIGHT / 2;
                  const botY2 = botY1 + CUT_TICK_HEIGHT;
                  const midY = featureLaneY + featureLaneHeight / 2;
                  return (
                    <g key={`cut-${i}`} data-cut style={{ pointerEvents: "none" }}>
                      {cut.onRowTop && (
                        <line
                          x1={cut.topX}
                          y1={topY1}
                          x2={cut.topX}
                          y2={topY2}
                          stroke={cut.color}
                          strokeWidth={1.2}
                          strokeLinecap="round"
                        />
                      )}
                      {cut.onRowBottom && (
                        <line
                          x1={cut.bottomX}
                          y1={botY1}
                          x2={cut.bottomX}
                          y2={botY2}
                          stroke={cut.color}
                          strokeWidth={1.2}
                          strokeLinecap="round"
                        />
                      )}
                      {!cut.blunt && cut.onRowTop && cut.onRowBottom && (
                        <>
                          <line
                            x1={cut.topX}
                            y1={topY1}
                            x2={cut.topX}
                            y2={midY}
                            stroke={cut.color}
                            strokeWidth={1.2}
                            strokeLinecap="round"
                          />
                          <line
                            x1={cut.topX}
                            y1={midY}
                            x2={cut.bottomX}
                            y2={midY}
                            stroke={cut.color}
                            strokeWidth={1.2}
                            strokeLinecap="round"
                          />
                          <line
                            x1={cut.bottomX}
                            y1={midY}
                            x2={cut.bottomX}
                            y2={botY1}
                            stroke={cut.color}
                            strokeWidth={1.2}
                            strokeLinecap="round"
                          />
                        </>
                      )}
                      {cut.blunt && cut.onRowTop && cut.onRowBottom && (
                        <line
                          x1={cut.topX}
                          y1={topY1}
                          x2={cut.bottomX}
                          y2={botY1}
                          stroke={cut.color}
                          strokeWidth={1.2}
                          strokeLinecap="round"
                        />
                      )}
                    </g>
                  );
                })}
              </svg>

              {/* Right gutter — kept for layout symmetry; position label
                  removed (2026-04-29) to keep only the inline 10-bp ticks. */}
              <div style={{ width: GUTTER_W }} />
            </div>
          );
        })}
        </div>
      </div>

      {aaTooltip && (
        <div
          data-aa-tooltip
          style={{
            position: "absolute",
            left: aaTooltip.x,
            top: aaTooltip.y,
            zIndex: 10,
            background: "#ffffff",
            border: "1px solid #d0d4d8",
            borderRadius: 4,
            padding: "8px 10px",
            boxShadow: "0 2px 8px rgba(0,0,0,0.08)",
            fontFamily: '"Source Sans Pro", sans-serif',
            fontSize: 12.5,
            color: "#2d3135",
            minWidth: 180,
            lineHeight: 1.45,
            pointerEvents: "auto",
          }}
        >
          <div>
            <strong>CDS ORF</strong>{" "}
            ({aaTooltip.orfLen} aa): {aaTooltip.letter}{aaTooltip.aaIdx}
          </div>
          {aaTooltip.region && (
            <div>
              <strong>{aaTooltip.region.name}</strong>{" "}
              ({aaTooltip.region.aa_end - aaTooltip.region.aa_start + 1} aa):{" "}
              {aaTooltip.letter}{aaTooltip.aaIdx - aaTooltip.region.aa_start + 1}
            </div>
          )}
          <div style={{ color: "#5a6c7d" }}>
            {AA_FULL_NAME[aaTooltip.letter] || "Unknown"}
          </div>
        </div>
      )}
    </div>
  );
}
