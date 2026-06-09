"use client";

import React, {
  useRef,
  useState,
  useEffect,
  useMemo,
  useCallback,
} from "react";
import type { SeqVizAnnotation } from "./InteractiveSequenceViewer";

// ─────────────────────────────────────────────────────────────────────────────
// Types
// ─────────────────────────────────────────────────────────────────────────────

export interface InteractionParticipant {
  name: string;
  start: number | null;
  end: number | null;
  strand?: number | null;
  role?: string;
  sbo_role?: string;
  so_role?: string;
  external?: boolean;
}

export interface PlasmidInteraction {
  interaction_id: string;
  interaction_type: string;
  sbo_term?: string;
  rule_id?: string;
  confidence?: number;
  notes?: string;
  participants: InteractionParticipant[];
}

interface CircularPlasmidViewerProps {
  sequence: string;
  annotations: SeqVizAnnotation[];
  title?: string;
  height?: number;
  onAnnotationClick?: (annotation: SeqVizAnnotation, index: number) => void;
  onSelectionChange?: (selection: { start: number; end: number } | null) => void;
  selection?: { start: number; end: number } | null;
  interactions?: PlasmidInteraction[];
  showInteractions?: boolean;
  highlightedInteractionId?: string | null;
}

interface LabelPlacement {
  index: number;
  midAngle: number;
  isInside: boolean;
  labelRadius: number;
  yOffset: number;
  isReverse: boolean;
}

// ─────────────────────────────────────────────────────────────────────────────
// Constants
// ─────────────────────────────────────────────────────────────────────────────

const MIN_ZOOM = 0.5;
const MAX_ZOOM = 200;
const BASE_RADIUS = 150;
const BACKBONE_STROKE = 1; // Thin black line

// DNA sequence display gap - space reserved for nucleotide letters
// This creates a clear separation between features and the backbone/DNA area
const DNA_STRAND_GAP = 6; // Gap for each DNA strand (forward above, reverse below backbone)

// Feature sizing: at 1x zoom, features are 2.5x base height
// As zoom increases, height decreases proportionally (at 10x = 1/10th of 1x height)
const FEATURE_BASE_STROKE = 10.5; // Base stroke width
const FEATURE_MULTIPLIER_AT_1X = 2.5; // Features are 2.5x larger at 1x zoom

// Font sizing: 50% larger at 1x, scales down with zoom
const BASE_FONT_SIZE = 10;
const FONT_MULTIPLIER_AT_1X = 1.5;

const EXTERNAL_LABEL_OFFSET_FORWARD = 55;
const EXTERNAL_LABEL_OFFSET_REVERSE = -55;

// Zoom threshold for sequence display
const SEQUENCE_ZOOM_THRESHOLD = 4;

// Opacity values
const FEATURE_OPACITY = 0.5;
const BACKBONE_OPACITY = 0.3;

const CHAR_WIDTH_ESTIMATE = 6;

const COLORS = {
  backbone: "#000000", // Black backbone line
  text: "#dbefe7",
  textDark: "#105b39",
  feature: "#7C3AED",
  featureAlt: ["#4a7c59", "#6b5b95", "#3b82f6", "#f97316", "#22c55e", "#ef4444"],
  selection: "rgba(59, 130, 246, 0.4)",
  selectionStroke: "#3b82f6",
  labelLine: "#888888",
};

// Base pair colors
const BASE_COLORS: Record<string, string> = {
  A: "#22c55e", // green
  T: "#ef4444", // red
  G: "#f59e0b", // amber
  C: "#3b82f6", // blue
};

// Complementary base mapping
const COMPLEMENT: Record<string, string> = {
  A: "T",
  T: "A",
  G: "C",
  C: "G",
};

// ─────────────────────────────────────────────────────────────────────────────
// Geometry Utilities
// ─────────────────────────────────────────────────────────────────────────────

function bpToAngle(bp: number, totalLength: number): number {
  return (bp / totalLength) * 2 * Math.PI - Math.PI / 2;
}

function angleToBp(angle: number, totalLength: number): number {
  const TAU = 2 * Math.PI;
  const normalizedAngle = ((angle + Math.PI / 2) % TAU + TAU) % TAU;
  return Math.round((normalizedAngle / TAU) * totalLength) % totalLength;
}

function polarToCartesian(angle: number, radius: number): { x: number; y: number } {
  return {
    x: radius * Math.cos(angle),
    y: radius * Math.sin(angle),
  };
}

function cartesianToAngle(x: number, y: number): number {
  return Math.atan2(y, x);
}

function normalizeAngle(angle: number): number {
  const TAU = 2 * Math.PI;
  return ((angle % TAU) + TAU) % TAU;
}

function calculateArcSpan(startBp: number, endBp: number, totalLength: number): number {
  if (startBp <= endBp) {
    return (endBp - startBp) / totalLength;
  } else {
    return (totalLength - startBp + endBp) / totalLength;
  }
}

function calculateArcLength(startBp: number, endBp: number, totalLength: number, radius: number): number {
  const arcSpan = calculateArcSpan(startBp, endBp, totalLength);
  return arcSpan * 2 * Math.PI * radius;
}

function createFilledArcPath(
  startBp: number,
  endBp: number,
  totalLength: number,
  innerRadius: number,
  outerRadius: number
): string {
  const startAngle = bpToAngle(startBp, totalLength);
  const arcSpan = calculateArcSpan(startBp, endBp, totalLength);
  const endAngle = startAngle + arcSpan * 2 * Math.PI;
  const largeArcFlag = arcSpan > 0.5 ? 1 : 0;

  const outerStart = polarToCartesian(startAngle, outerRadius);
  const outerEnd = polarToCartesian(endAngle, outerRadius);
  const innerStart = polarToCartesian(startAngle, innerRadius);
  const innerEnd = polarToCartesian(endAngle, innerRadius);

  return `
    M ${outerStart.x.toFixed(3)} ${outerStart.y.toFixed(3)}
    A ${outerRadius} ${outerRadius} 0 ${largeArcFlag} 1 ${outerEnd.x.toFixed(3)} ${outerEnd.y.toFixed(3)}
    L ${innerEnd.x.toFixed(3)} ${innerEnd.y.toFixed(3)}
    A ${innerRadius} ${innerRadius} 0 ${largeArcFlag} 0 ${innerStart.x.toFixed(3)} ${innerStart.y.toFixed(3)}
    Z
  `;
}

// Create arrow-shaped (pentagon) path for directional features
function createArrowArcPath(
  startBp: number,
  endBp: number,
  totalLength: number,
  innerRadius: number,
  outerRadius: number,
  direction: number // 1 = forward (arrow at end), -1 = reverse (arrow at start)
): string {
  const startAngle = bpToAngle(startBp, totalLength);
  const arcSpan = calculateArcSpan(startBp, endBp, totalLength);
  const endAngle = startAngle + arcSpan * 2 * Math.PI;
  const largeArcFlag = arcSpan > 0.5 ? 1 : 0;

  const midRadius = (innerRadius + outerRadius) / 2;
  // Arrow length scales with feature size but capped for large features
  // Smaller features get proportionally larger arrows, large features get subtle arrows
  const arrowLength = Math.min(arcSpan * 0.15, 0.015) * 2 * Math.PI; // Arrow takes up to 15% of arc or 1.5% of circle

  if (direction === 1) {
    // Forward: arrow points at the end
    const arrowAngle = endAngle - arrowLength;
    const outerStart = polarToCartesian(startAngle, outerRadius);
    const outerArrowBase = polarToCartesian(arrowAngle, outerRadius);
    const arrowTip = polarToCartesian(endAngle, midRadius);
    const innerArrowBase = polarToCartesian(arrowAngle, innerRadius);
    const innerStart = polarToCartesian(startAngle, innerRadius);

    const bodyLargeArc = (arrowAngle - startAngle) / (2 * Math.PI) > 0.5 ? 1 : 0;

    return `
      M ${outerStart.x.toFixed(3)} ${outerStart.y.toFixed(3)}
      A ${outerRadius} ${outerRadius} 0 ${bodyLargeArc} 1 ${outerArrowBase.x.toFixed(3)} ${outerArrowBase.y.toFixed(3)}
      L ${arrowTip.x.toFixed(3)} ${arrowTip.y.toFixed(3)}
      L ${innerArrowBase.x.toFixed(3)} ${innerArrowBase.y.toFixed(3)}
      A ${innerRadius} ${innerRadius} 0 ${bodyLargeArc} 0 ${innerStart.x.toFixed(3)} ${innerStart.y.toFixed(3)}
      Z
    `;
  } else {
    // Reverse: arrow points at the start
    const arrowAngle = startAngle + arrowLength;
    const outerEnd = polarToCartesian(endAngle, outerRadius);
    const outerArrowBase = polarToCartesian(arrowAngle, outerRadius);
    const arrowTip = polarToCartesian(startAngle, midRadius);
    const innerArrowBase = polarToCartesian(arrowAngle, innerRadius);
    const innerEnd = polarToCartesian(endAngle, innerRadius);

    const bodySpan = (endAngle - arrowAngle) / (2 * Math.PI);
    const bodyLargeArc = bodySpan > 0.5 ? 1 : 0;

    return `
      M ${outerArrowBase.x.toFixed(3)} ${outerArrowBase.y.toFixed(3)}
      A ${outerRadius} ${outerRadius} 0 ${bodyLargeArc} 1 ${outerEnd.x.toFixed(3)} ${outerEnd.y.toFixed(3)}
      L ${innerEnd.x.toFixed(3)} ${innerEnd.y.toFixed(3)}
      A ${innerRadius} ${innerRadius} 0 ${bodyLargeArc} 0 ${innerArrowBase.x.toFixed(3)} ${innerArrowBase.y.toFixed(3)}
      L ${arrowTip.x.toFixed(3)} ${arrowTip.y.toFixed(3)}
      Z
    `;
  }
}

function createCirclePath(radius: number): string {
  return `M 0 ${-radius} A ${radius} ${radius} 0 1 1 0 ${radius} A ${radius} ${radius} 0 1 1 0 ${-radius}`;
}

// Create an arc path for curved text (textPath)
function createTextArcPath(
  startBp: number,
  endBp: number,
  totalLength: number,
  radius: number
): string {
  const startAngle = bpToAngle(startBp, totalLength);
  const arcSpan = calculateArcSpan(startBp, endBp, totalLength);
  const endAngle = startAngle + arcSpan * 2 * Math.PI;
  const largeArcFlag = arcSpan > 0.5 ? 1 : 0;

  const start = polarToCartesian(startAngle, radius);
  const end = polarToCartesian(endAngle, radius);

  return `M ${start.x.toFixed(3)} ${start.y.toFixed(3)} A ${radius} ${radius} 0 ${largeArcFlag} 1 ${end.x.toFixed(3)} ${end.y.toFixed(3)}`;
}

function getAnnotationMidAngle(startBp: number, endBp: number, totalLength: number): number {
  const startAngle = bpToAngle(startBp, totalLength);
  const arcSpan = calculateArcSpan(startBp, endBp, totalLength);
  return startAngle + (arcSpan * Math.PI);
}

// Vivid palette for CRISPR/cloning-pipeline-added features (pegRNAs,
// ngRNAs, RT templates, primers, edit sites). Cycles by annotation index
// so each design-added feature gets a distinct hue. Matches the same
// palette in FlatPlasmidViewer + LinearSequenceViewer.
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

// Muted gray shades for the original .gb annotations preserved alongside
// design-added features — keeps the genomic backdrop legible without
// competing with the vivid design overlay.
const PRESERVED_GRAY_PALETTE = ["#9ca3af", "#6b7280", "#cbd5e1"];

function getAnnotationColor(
  annotation: SeqVizAnnotation,
  index: number,
  hasDesignAdditions: boolean,
): string {
  if ((annotation as any).added_by_design) {
    return DESIGN_PALETTE[index % DESIGN_PALETTE.length];
  }
  // When viewing a design-output .gb, mute the pre-existing annotations to
  // gray so the design additions read at a glance. Plain plasmid views
  // (no added_by_design features) keep the original palette.
  if (hasDesignAdditions) {
    return PRESERVED_GRAY_PALETTE[index % PRESERVED_GRAY_PALETTE.length];
  }
  if (annotation.color) return annotation.color;
  if (annotation.layer === "module") return "#4a7c59";
  if (annotation.layer === "motif") return "#6b5b95";
  return COLORS.featureAlt[index % COLORS.featureAlt.length];
}

function getAnnotationLength(start: number, end: number, totalLength: number): number {
  if (start <= end) {
    return end - start;
  } else {
    return (totalLength - start) + end;
  }
}

function annotationsOverlap(
  a: { start: number; end: number },
  b: { start: number; end: number },
  totalLength: number
): boolean {
  // Convert to normalized ranges for overlap detection
  const aSpansOrigin = a.start > a.end;
  const bSpansOrigin = b.start > b.end;

  if (!aSpansOrigin && !bSpansOrigin) {
    // Neither spans origin - simple overlap check
    return a.start < b.end && b.start < a.end;
  } else if (aSpansOrigin && bSpansOrigin) {
    // Both span origin - they always overlap
    return true;
  } else {
    // One spans origin, one doesn't
    const spanning = aSpansOrigin ? a : b;
    const normal = aSpansOrigin ? b : a;
    // Spanning annotation covers [spanning.start, totalLength) and [0, spanning.end)
    return normal.start < spanning.end || normal.end > spanning.start;
  }
}

// ─────────────────────────────────────────────────────────────────────────────
// Gesture Handlers
// ─────────────────────────────────────────────────────────────────────────────

function getTouchDistance(touches: TouchList): number {
  if (touches.length < 2) return 0;
  const dx = touches[1].clientX - touches[0].clientX;
  const dy = touches[1].clientY - touches[0].clientY;
  return Math.sqrt(dx * dx + dy * dy);
}

function clamp(value: number, min: number, max: number): number {
  return Math.max(min, Math.min(max, value));
}

// ─────────────────────────────────────────────────────────────────────────────
// Component
// ─────────────────────────────────────────────────────────────────────────────

export default function CircularPlasmidViewer({
  sequence,
  annotations,
  title = "Plasmid",
  height = 500,
  onAnnotationClick,
  onSelectionChange,
  interactions = [],
  showInteractions = true,
  highlightedInteractionId = null,
  selection = null,
}: CircularPlasmidViewerProps) {
  const svgRef = useRef<SVGSVGElement>(null);
  const containerRef = useRef<HTMLDivElement>(null);

  const [rotation, setRotation] = useState(0);
  const [zoom, setZoom] = useState(1);
  const [hoveredAnnotation, setHoveredAnnotation] = useState<number | null>(null);

  const [isDragging, setIsDragging] = useState(false);
  const [selectionStart, setSelectionStart] = useState<number | null>(null);
  const [selectionEnd, setSelectionEnd] = useState<number | null>(null);
  const [selectionDirection, setSelectionDirection] = useState<'forward' | 'reverse'>('forward');

  // Sync external selection (linear-side selection / find match) into
  // internal state so the circular viewer paints the same highlight the
  // linear viewer is showing.
  useEffect(() => {
    if (selection && typeof selection.start === 'number' && typeof selection.end === 'number') {
      if (selection.start === selectionStart && selection.end === selectionEnd) return;
      setSelectionStart(selection.start);
      setSelectionEnd(selection.end);
    } else if (!selection && (selectionStart !== null || selectionEnd !== null)) {
      setSelectionStart(null);
      setSelectionEnd(null);
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [selection?.start, selection?.end]);
  const [isShiftHeld, setIsShiftHeld] = useState(false);
  const dragStartBp = useRef<number | null>(null);
  const lastDragBp = useRef<number | null>(null);
  const cumulativeDragAngle = useRef<number>(0);

  const lastTouchDist = useRef<number>(0);
  const isTouching = useRef(false);

  // Track shift key state
  useEffect(() => {
    const handleKeyDown = (e: KeyboardEvent) => {
      if (e.key === "Shift") setIsShiftHeld(true);
    };
    const handleKeyUp = (e: KeyboardEvent) => {
      if (e.key === "Shift") setIsShiftHeld(false);
    };
    window.addEventListener("keydown", handleKeyDown);
    window.addEventListener("keyup", handleKeyUp);
    return () => {
      window.removeEventListener("keydown", handleKeyDown);
      window.removeEventListener("keyup", handleKeyUp);
    };
  }, []);

  const totalLength = sequence.length;

  // Feature stroke width: 2.5x at zoom 1, scales down inversely with zoom
  // At zoom 1: 2.5 * 14 = 35
  // At zoom 10: 35 / 10 = 3.5
  const featureStrokeWidth = useMemo(() => {
    return (FEATURE_MULTIPLIER_AT_1X * FEATURE_BASE_STROKE) / zoom;
  }, [zoom]);

  // Label font size: 50% larger at 1x, scales down with zoom
  const labelFontSize = useMemo(() => {
    const size = (FONT_MULTIPLIER_AT_1X * BASE_FONT_SIZE) / zoom;
    return Math.max(1.5, Math.min(size, 20)); // Floor matches feature-label floor so ticks scale with zoom
  }, [zoom]);

  // Whether to show DNA sequence (at 8x+ zoom)
  const shouldShowSequence = useMemo(() => {
    return zoom >= SEQUENCE_ZOOM_THRESHOLD;
  }, [zoom]);

  // ───────────────────────────────────────────────────────────────────────────
  // Mouse Position to BP Conversion
  // ───────────────────────────────────────────────────────────────────────────

  const getMouseBp = useCallback((e: MouseEvent | React.MouseEvent): number | null => {
    const svg = svgRef.current;
    if (!svg) return null;

    // Use SVG's native coordinate transformation for accurate conversion
    const pt = svg.createSVGPoint();
    pt.x = e.clientX;
    pt.y = e.clientY;

    // Get the inverse of the screen transformation matrix
    const screenCTM = svg.getScreenCTM();
    if (!screenCTM) return null;

    // Transform screen coordinates to SVG coordinates
    const svgPt = pt.matrixTransform(screenCTM.inverse());

    // Now we have coordinates in the SVG viewBox space
    // But we still need to undo the rotation transform applied to the <g> element
    const cosR = Math.cos(-rotation);
    const sinR = Math.sin(-rotation);
    const unrotatedX = svgPt.x * cosR - svgPt.y * sinR;
    const unrotatedY = svgPt.x * sinR + svgPt.y * cosR;

    // Check if click is near the plasmid backbone (with generous tolerance)
    const distFromCenter = Math.sqrt(unrotatedX * unrotatedX + unrotatedY * unrotatedY);
    const tolerance = Math.max(60, featureStrokeWidth * 2);
    if (distFromCenter < BASE_RADIUS - tolerance || distFromCenter > BASE_RADIUS + tolerance) {
      return null;
    }

    const angle = cartesianToAngle(unrotatedX, unrotatedY);
    return angleToBp(angle, totalLength);
  }, [rotation, totalLength, featureStrokeWidth]);

  // Unbounded variant used for resolving feature clicks off the backbone.
  const getMouseBpUnbounded = useCallback((e: MouseEvent | React.MouseEvent): number | null => {
    const svg = svgRef.current;
    if (!svg) return null;
    const pt = svg.createSVGPoint();
    pt.x = e.clientX;
    pt.y = e.clientY;
    const screenCTM = svg.getScreenCTM();
    if (!screenCTM) return null;
    const svgPt = pt.matrixTransform(screenCTM.inverse());
    const cosR = Math.cos(-rotation);
    const sinR = Math.sin(-rotation);
    const unrotatedX = svgPt.x * cosR - svgPt.y * sinR;
    const unrotatedY = svgPt.x * sinR + svgPt.y * cosR;
    const angle = cartesianToAngle(unrotatedX, unrotatedY);
    return angleToBp(angle, totalLength);
  }, [rotation, totalLength]);

  // ───────────────────────────────────────────────────────────────────────────
  // Selection Handlers
  // ───────────────────────────────────────────────────────────────────────────

  const handleMouseDown = useCallback((e: React.MouseEvent) => {
    if ((e.target as Element).closest('[data-annotation]')) return;

    const bp = getMouseBp(e);
    if (bp !== null) {
      if (isShiftHeld && selectionStart !== null) {
        // Shift+click: extend selection from current start to this position
        setSelectionEnd(bp);
        onSelectionChange?.({ start: selectionStart, end: bp });
      } else {
        // Normal click: start new selection
        setIsDragging(true);
        dragStartBp.current = bp;
        lastDragBp.current = bp;
        cumulativeDragAngle.current = 0;
        setSelectionStart(bp);
        setSelectionEnd(bp);
        setSelectionDirection('forward');
      }
    }
  }, [getMouseBp, isShiftHeld, selectionStart, onSelectionChange]);

  // Handle clicking on an annotation to select its region
  const handleAnnotationSelect = useCallback((annotation: SeqVizAnnotation, e: React.MouseEvent) => {
    if (isShiftHeld && selectionStart !== null) {
      // Shift+click on feature: extend selection to include this feature
      // Extend to the furthest point of the feature from current selection start
      const distToStart = getAnnotationLength(selectionStart, annotation.start, totalLength);
      const distToEnd = getAnnotationLength(selectionStart, annotation.end, totalLength);
      const newEnd = distToEnd > distToStart ? annotation.end : annotation.start;
      setSelectionEnd(newEnd);
      onSelectionChange?.({ start: selectionStart, end: newEnd });
    } else {
      // Normal click on feature: select the feature's region
      setSelectionStart(annotation.start);
      setSelectionEnd(annotation.end);
      onSelectionChange?.({ start: annotation.start, end: annotation.end });
    }
  }, [isShiftHeld, selectionStart, totalLength, onSelectionChange]);

  const handleMouseMove = useCallback((e: MouseEvent) => {
    if (!isDragging || dragStartBp.current === null || lastDragBp.current === null) return;

    const bp = getMouseBp(e);
    if (bp !== null) {
      // Calculate angle change from the last position
      const lastAngle = bpToAngle(lastDragBp.current, totalLength);
      const currentAngle = bpToAngle(bp, totalLength);

      // Calculate the shortest angle difference
      let angleDiff = currentAngle - lastAngle;
      if (angleDiff > Math.PI) angleDiff -= 2 * Math.PI;
      if (angleDiff < -Math.PI) angleDiff += 2 * Math.PI;

      cumulativeDragAngle.current += angleDiff;
      lastDragBp.current = bp;

      // Calculate both possible selection spans
      const startBp = dragStartBp.current;
      const forwardSpan = bp >= startBp ? bp - startBp : totalLength - startBp + bp;
      const reverseSpan = startBp >= bp ? startBp - bp : totalLength - bp + startBp;

      // Determine direction based on cumulative angle AND which span makes more sense
      // Use a hysteresis threshold to prevent flipping back and forth
      const hysteresis = 0.1; // radians
      let newDirection: 'forward' | 'reverse';

      if (Math.abs(cumulativeDragAngle.current) < hysteresis) {
        // Very small movement - use the shorter span
        newDirection = forwardSpan <= reverseSpan ? 'forward' : 'reverse';
      } else {
        // Use the cumulative angle direction
        // Positive = clockwise = forward, Negative = counter-clockwise = reverse
        newDirection = cumulativeDragAngle.current >= 0 ? 'forward' : 'reverse';
      }

      setSelectionDirection(newDirection);

      // Set selection based on direction
      // For forward: start at dragStart, end at current bp (arc goes clockwise)
      // For reverse: start at current bp, end at dragStart (arc goes clockwise from current to dragStart)
      // This ensures the arc always covers the bases in the direction the user dragged
      if (newDirection === 'forward') {
        setSelectionStart(startBp);
        setSelectionEnd(bp);
      } else {
        // Reverse: we want to select from current position going forward to dragStart
        // This means setting start=bp, end=dragStart
        setSelectionStart(bp);
        setSelectionEnd(startBp);
      }
    }
  }, [isDragging, getMouseBp, totalLength]);

  const handleMouseUp = useCallback(() => {
    if (isDragging && selectionStart !== null && selectionEnd !== null) {
      if (selectionStart !== selectionEnd) {
        onSelectionChange?.({ start: selectionStart, end: selectionEnd });
      } else {
        setSelectionStart(null);
        setSelectionEnd(null);
        onSelectionChange?.(null);
      }
    }
    setIsDragging(false);
    dragStartBp.current = null;
  }, [isDragging, selectionStart, selectionEnd, onSelectionChange]);

  // ───────────────────────────────────────────────────────────────────────────
  // Gesture Event Handlers
  // ───────────────────────────────────────────────────────────────────────────

  const handleWheel = useCallback((e: WheelEvent) => {
    e.preventDefault();

    if (e.ctrlKey || e.metaKey) {
      setZoom((z) => clamp(z - e.deltaY * 0.04, MIN_ZOOM, MAX_ZOOM));
    } else {
      // Vertical scroll rotates the plasmid.
      setRotation((r) => r + e.deltaY * 0.002);
    }
  }, []);

  const handleTouchStart = useCallback((e: TouchEvent) => {
    if (e.touches.length === 2) {
      isTouching.current = true;
      lastTouchDist.current = getTouchDistance(e.touches);
    }
  }, []);

  const handleTouchMove = useCallback((e: TouchEvent) => {
    if (e.touches.length === 2 && isTouching.current) {
      e.preventDefault();
      const newDist = getTouchDistance(e.touches);
      if (lastTouchDist.current > 0) {
        const scale = lastTouchDist.current / newDist;
        setZoom((z) => clamp(z * scale, MIN_ZOOM, MAX_ZOOM));
      }
      lastTouchDist.current = newDist;
    }
  }, []);

  const handleTouchEnd = useCallback(() => {
    isTouching.current = false;
    lastTouchDist.current = 0;
  }, []);

  // ───────────────────────────────────────────────────────────────────────────
  // Event Listener Setup
  // ───────────────────────────────────────────────────────────────────────────

  useEffect(() => {
    const container = containerRef.current;
    if (!container) return;

    container.addEventListener("wheel", handleWheel, { passive: false });
    container.addEventListener("touchstart", handleTouchStart, { passive: true });
    container.addEventListener("touchmove", handleTouchMove, { passive: false });
    container.addEventListener("touchend", handleTouchEnd, { passive: true });

    return () => {
      container.removeEventListener("wheel", handleWheel);
      container.removeEventListener("touchstart", handleTouchStart);
      container.removeEventListener("touchmove", handleTouchMove);
      container.removeEventListener("touchend", handleTouchEnd);
    };
  }, [handleWheel, handleTouchStart, handleTouchMove, handleTouchEnd]);

  useEffect(() => {
    if (isDragging) {
      document.addEventListener("mousemove", handleMouseMove);
      document.addEventListener("mouseup", handleMouseUp);
      return () => {
        document.removeEventListener("mousemove", handleMouseMove);
        document.removeEventListener("mouseup", handleMouseUp);
      };
    }
  }, [isDragging, handleMouseMove, handleMouseUp]);

  // ───────────────────────────────────────────────────────────────────────────
  // Computed Values
  // ───────────────────────────────────────────────────────────────────────────

  const viewBoxSize = useMemo(() => {
    const baseSize = (BASE_RADIUS + EXTERNAL_LABEL_OFFSET_FORWARD + 50) * 2;
    return baseSize / zoom;
  }, [zoom]);

  const viewBoxOffset = useMemo(() => {
    if (zoom <= 1) return { x: 0, y: 0 };
    // Focal target sits ~10% outside the backbone so the DNA sequence (which
    // straddles the backbone and extends outward into the feature area) is
    // more vertically centered in the zoomed viewport.
    const targetY = -BASE_RADIUS * 1.07;
    const offsetY = targetY * (1 - 1 / zoom);
    return { x: 0, y: offsetY };
  }, [zoom]);

  const visibleBases = useMemo(() => {
    if (!shouldShowSequence || !sequence || totalLength === 0) return [];

    const visibleFraction = 1 / zoom;
    const basesVisible = Math.ceil(totalLength * visibleFraction * 1.5);

    // Ensure rotationBp stays within bounds (Math.round can produce totalLength)
    const rotationBp = Math.round(
      (normalizeAngle(-rotation) / (2 * Math.PI)) * totalLength
    ) % totalLength;

    const halfVisible = Math.floor(basesVisible / 2);
    const startBp = (rotationBp - halfVisible + totalLength) % totalLength;

    const bases: { bp: number; base: string; complement: string; angle: number }[] = [];
    for (let i = 0; i < Math.min(basesVisible, totalLength); i++) {
      const bp = (startBp + i) % totalLength;
      // Ensure bp is valid before accessing sequence
      if (bp < 0 || bp >= sequence.length) continue;
      const base = (sequence[bp] || "N").toUpperCase();
      const complement = COMPLEMENT[base] || "N";
      bases.push({
        bp,
        base,
        complement,
        angle: bpToAngle(bp, totalLength),
      });
    }

    return bases;
  }, [shouldShowSequence, zoom, rotation, totalLength, sequence]);

  // ───────────────────────────────────────────────────────────────────────────
  // Annotation Data and Label Placement
  // ───────────────────────────────────────────────────────────────────────────

  const backbonePath = useMemo(() => createCirclePath(BASE_RADIUS), []);

  // Tick marks - computed outside JSX to avoid hook rules violation.
  // Interval is zoom-adaptive: 1x=1000, 2x=500, 5x=200, and finer as zoom increases,
  // picked from a 1/2/5 "nice number" ladder and scaled by plasmid size.
  const tickMarks = useMemo(() => {
    const baseInterval =
      totalLength > 50000 ? 5000 :
      totalLength > 20000 ? 2000 :
      totalLength > 8000  ? 1000 :
      totalLength > 2000  ? 500  :
                            100;
    const niceSteps = [10000, 5000, 2000, 1000, 500, 200, 100, 50];
    const target = baseInterval / zoom;
    let interval = niceSteps[niceSteps.length - 1];
    for (const s of niceSteps) {
      if (s <= target * 1.5) { interval = s; break; }
    }

    const ticks: { bp: number; angle: number; inner: { x: number; y: number }; outer: { x: number; y: number }; labelPos: { x: number; y: number } }[] = [];

    for (let bp = 0; bp < totalLength; bp += interval) {
      const angle = bpToAngle(bp, totalLength);
      ticks.push({
        bp,
        angle,
        inner: polarToCartesian(angle, BASE_RADIUS - 2 / zoom),
        outer: polarToCartesian(angle, BASE_RADIUS - 8 / zoom),
        labelPos: polarToCartesian(angle, BASE_RADIUS - 16 / zoom),
      });
    }

    return ticks;
  }, [totalLength, zoom]);

  const annotationData = useMemo(() => {
    const hasDesignAdditions = annotations.some(
      (a) => (a as any).added_by_design === true,
    );
    // First pass: calculate basic data and lengths for all annotations
    const basicData = annotations.map((ann, idx) => ({
      annotation: ann,
      index: idx,
      isReverse: ann.direction === -1,
      length: getAnnotationLength(ann.start, ann.end, totalLength),
      midAngle: getAnnotationMidAngle(ann.start, ann.end, totalLength),
      color: getAnnotationColor(ann, idx, hasDesignAdditions),
      spansOrigin: ann.start > ann.end,
    }));

    // Separate by direction
    const forwardAnns = basicData.filter(d => !d.isReverse);
    const reverseAnns = basicData.filter(d => d.isReverse);

    // Sort each group by length (longest first - they get base layer)
    forwardAnns.sort((a, b) => b.length - a.length);
    reverseAnns.sort((a, b) => b.length - a.length);

    // Assign layer levels based on overlaps
    // Layer 0 = base layer (anchored to backbone)
    // Higher layers = further from backbone
    const assignLayers = (anns: typeof basicData) => {
      const layers: Map<number, number> = new Map();

      for (const ann of anns) {
        // Find the minimum layer that doesn't overlap with any annotation in that layer
        let layer = 0;
        let foundLayer = false;

        while (!foundLayer) {
          const annsInLayer = anns.filter(a => layers.get(a.index) === layer);
          const hasOverlap = annsInLayer.some(a =>
            annotationsOverlap(ann.annotation, a.annotation, totalLength)
          );

          if (!hasOverlap) {
            foundLayer = true;
          } else {
            layer++;
          }
        }

        layers.set(ann.index, layer);
      }

      return layers;
    };

    const forwardLayers = assignLayers(forwardAnns);
    const reverseLayers = assignLayers(reverseAnns);

    // Second pass: compute radii with layer offsets
    return basicData.map((data) => {
      const { annotation, index, isReverse, midAngle, color, spansOrigin } = data;

      const layer = isReverse ? (reverseLayers.get(index) || 0) : (forwardLayers.get(index) || 0);
      const layerOffset = layer * featureStrokeWidth;

      let innerR: number;
      let outerR: number;

      if (isReverse) {
        // Reverse: top of feature anchored below the reverse DNA strand
        // Leave room for DNA display between backbone and features
        // Higher layers go further inward (lower radius)
        outerR = BASE_RADIUS - BACKBONE_STROKE / 2 / zoom - DNA_STRAND_GAP / zoom - layerOffset;
        innerR = outerR - featureStrokeWidth;
      } else {
        // Forward or none: bottom of feature anchored above the forward DNA strand
        // Leave room for DNA display between backbone and features
        // Higher layers go further outward (higher radius)
        innerR = BASE_RADIUS + BACKBONE_STROKE / 2 / zoom + DNA_STRAND_GAP / zoom + layerOffset;
        outerR = innerR + featureStrokeWidth;
      }

      const featureRadius = (innerR + outerR) / 2;
      const arcLength = calculateArcLength(annotation.start, annotation.end, totalLength, featureRadius);

      const labelText = annotation.name + (spansOrigin ? " ⟳" : "");
      const estimatedLabelWidth = labelText.length * (labelFontSize * 0.6);
      const labelFitsInside = arcLength > estimatedLabelWidth * 1.3;

      return {
        annotation,
        index,
        featureRadius,
        innerR,
        outerR,
        midAngle,
        arcLength,
        color,
        spansOrigin,
        isReverse,
        labelFitsInside,
        labelText,
        layer,
      };
    });
  }, [annotations, totalLength, labelFontSize, featureStrokeWidth]);

  // Resolve a click at a given bp to the most-specific (shortest) feature that
  // contains it. Prevents double-clicking an overlap/sub-feature from surfacing
  // the larger underlying feature it's layered on top of.
  const resolveTopmostFeatureAt = useCallback((bp: number) => {
    const containing = annotationData.filter(d => {
      const s = d.annotation.start;
      const e = d.annotation.end;
      return s <= e ? (bp >= s && bp <= e) : (bp >= s || bp <= e);
    });
    if (containing.length === 0) return null;
    containing.sort((a, b) =>
      getAnnotationLength(a.annotation.start, a.annotation.end, totalLength) -
      getAnnotationLength(b.annotation.start, b.annotation.end, totalLength)
    );
    return { annotation: containing[0].annotation, index: containing[0].index };
  }, [annotationData, totalLength]);

  const labelPlacements = useMemo(() => {
    const placements: LabelPlacement[] = [];
    // Track all external labels with their positions for collision detection
    const externalLabels: { angle: number; yOffset: number; isReverse: boolean; labelWidth: number }[] = [];

    // Sort annotations by arc length (longer features get priority for inside labels)
    const sortedData = [...annotationData].sort((a, b) => b.arcLength - a.arcLength);

    for (const data of sortedData) {
      if (data.labelFitsInside) {
        placements.push({
          index: data.index,
          midAngle: data.midAngle,
          isInside: true,
          labelRadius: data.featureRadius,
          yOffset: 0,
          isReverse: data.isReverse,
        });
      } else {
        const baseOffset = data.isReverse ? EXTERNAL_LABEL_OFFSET_REVERSE : EXTERNAL_LABEL_OFFSET_FORWARD;
        const estimatedLabelWidth = data.labelText.length * (labelFontSize * 0.6);

        // More aggressive overlap detection
        // Convert label width to approximate angular span at the label radius
        const labelRadius = BASE_RADIUS + baseOffset;
        const angularWidth = estimatedLabelWidth / labelRadius;
        const offsetStep = labelFontSize + 6;

        let yOffset = 0;
        let foundSlot = false;
        const maxLayers = 8; // Maximum number of stacking layers

        for (let layer = 0; layer < maxLayers && !foundSlot; layer++) {
          const testOffset = layer * offsetStep * (data.isReverse ? -1 : 1);
          let hasCollision = false;

          for (const existing of externalLabels) {
            // Only check collision with labels on the same side
            if (existing.isReverse !== data.isReverse) continue;
            // Only check collision with labels at the same offset level
            if (Math.abs(existing.yOffset - testOffset) > 2) continue;

            // Calculate angular distance
            const angleDiff = Math.abs(normalizeAngle(data.midAngle) - normalizeAngle(existing.angle));
            const effectiveAngleDiff = Math.min(angleDiff, 2 * Math.PI - angleDiff);

            // Check if labels would overlap (considering both labels' widths)
            const minAngularDistance = (angularWidth + existing.labelWidth / labelRadius) / 2 + 0.05;

            if (effectiveAngleDiff < minAngularDistance) {
              hasCollision = true;
              break;
            }
          }

          if (!hasCollision) {
            yOffset = testOffset;
            foundSlot = true;
          }
        }

        externalLabels.push({
          angle: data.midAngle,
          yOffset,
          isReverse: data.isReverse,
          labelWidth: estimatedLabelWidth
        });

        placements.push({
          index: data.index,
          midAngle: data.midAngle,
          isInside: false,
          labelRadius: labelRadius,
          yOffset,
          isReverse: data.isReverse,
        });
      }
    }

    return placements;
  }, [annotationData, labelFontSize]);

  const selectionLength = useMemo(() => {
    if (selectionStart === null || selectionEnd === null) return 0;
    if (selectionStart <= selectionEnd) {
      return selectionEnd - selectionStart;
    } else {
      return (totalLength - selectionStart) + selectionEnd;
    }
  }, [selectionStart, selectionEnd, totalLength]);

  // Sequence font size - sized to fit between bases without overlap
  const sequenceFontSize = useMemo(() => {
    // Calculate the angular spacing between bases at the backbone
    const bpAngle = (2 * Math.PI) / totalLength;
    // Pixel spacing between bases in SVG coordinates
    const bpSpacing = bpAngle * BASE_RADIUS;

    // Font size should be slightly smaller than spacing to prevent overlap
    // but large enough to read. Use 85% of spacing.
    const fontSize = bpSpacing * 0.85;

    // Clamp to reasonable bounds
    // Min: readable at high zoom (small in SVG coords = larger on screen when zoomed)
    // Max: not too large at low zoom
    return Math.max(0.05, Math.min(fontSize, 6));
  }, [totalLength]);

  // Compute which features are visible in the current view
  const visibleFeatures = useMemo(() => {
    if (zoom < 2) return []; // Only show when zoomed in enough

    // Calculate the visible angle range based on zoom and rotation
    const visibleAngleRange = (2 * Math.PI) / zoom;
    const centerAngle = normalizeAngle(-rotation - Math.PI / 2); // Top of view
    const startAngle = normalizeAngle(centerAngle - visibleAngleRange / 2);
    const endAngle = normalizeAngle(centerAngle + visibleAngleRange / 2);

    return annotationData.filter(data => {
      const featureStartAngle = normalizeAngle(bpToAngle(data.annotation.start, totalLength));
      const featureEndAngle = normalizeAngle(bpToAngle(data.annotation.end, totalLength));
      const featureMidAngle = normalizeAngle(data.midAngle);

      // Check if feature midpoint is in visible range
      const isInRange = (angle: number) => {
        if (startAngle < endAngle) {
          return angle >= startAngle && angle <= endAngle;
        } else {
          // Range wraps around 0
          return angle >= startAngle || angle <= endAngle;
        }
      };

      return isInRange(featureMidAngle) || isInRange(featureStartAngle) || isInRange(featureEndAngle);
    });
  }, [annotationData, zoom, rotation, totalLength]);

  // Compute label positions with collision detection for visible features.
  // Labels are drawn counter-rotated so they stay axis-aligned on screen;
  // collision is therefore checked in screen-space axis-aligned boxes with
  // explicit padding to keep visible gaps between neighbouring labels.
  const visibleFeatureLabelOffsets = useMemo(() => {
    if (visibleFeatures.length === 0) return new Map<number, number>();

    const offsets = new Map<number, number>();
    const zoomedFontSize = Math.max(1.5, Math.min(12 / zoom, featureStrokeWidth * 0.6));
    const labelHeight = zoomedFontSize * 1.2;
    const hPad = zoomedFontSize * 0.8; // enforced horizontal gap between labels
    const vPad = zoomedFontSize * 0.4; // enforced vertical gap between stacked rows

    const cosR = Math.cos(rotation);
    const sinR = Math.sin(rotation);

    type PlacedLabel = { screenX: number; screenY: number; width: number; height: number; isReverse: boolean };
    const placedLabels: PlacedLabel[] = [];

    // Sort by arc length so larger features get priority for the closest slot
    const sortedFeatures = [...visibleFeatures].sort((a, b) => b.arcLength - a.arcLength);

    for (const data of sortedFeatures) {
      const { annotation, midAngle, isReverse, innerR, outerR, featureRadius } = data;

      const showOutside = shouldShowSequence || (!data.labelFitsInside && zoom >= 8);

      if (!showOutside) {
        offsets.set(data.index, 0);
        continue;
      }

      // Anchor collision grid to the feature midline (matches the render).
      const baseRadius = featureRadius;
      const estimatedWidth = annotation.name.length * zoomedFontSize * 0.6;

      let yOffset = 0;
      let foundSlot = false;
      const maxLayers = 10;
      const offsetStep = labelHeight + vPad;

      for (let layer = 0; layer < maxLayers && !foundSlot; layer++) {
        const testOffset = layer * offsetStep * (isReverse ? -1 : 1);
        const adjustedRadius = baseRadius + testOffset;
        const pos = polarToCartesian(midAngle, adjustedRadius);
        // Transform into the on-screen (unrotated) frame so axis-aligned
        // box overlap matches what the user actually sees.
        const screenX = pos.x * cosR - pos.y * sinR;
        const screenY = pos.x * sinR + pos.y * cosR;

        let hasCollision = false;
        for (const placed of placedLabels) {
          const dx = Math.abs(screenX - placed.screenX);
          const dy = Math.abs(screenY - placed.screenY);
          const minDx = (estimatedWidth + placed.width) / 2 + hPad;
          const minDy = (labelHeight + placed.height) / 2 + vPad;
          if (dx < minDx && dy < minDy) {
            hasCollision = true;
            break;
          }
        }

        if (!hasCollision) {
          yOffset = testOffset;
          foundSlot = true;
          placedLabels.push({
            screenX,
            screenY,
            width: estimatedWidth,
            height: labelHeight,
            isReverse,
          });
        }
      }

      // Even if every slot collides, still record the final offset so the
      // label renders at the furthest layer we tried.
      if (!foundSlot) {
        const testOffset = (maxLayers - 1) * offsetStep * (isReverse ? -1 : 1);
        const adjustedRadius = baseRadius + testOffset;
        const pos = polarToCartesian(midAngle, adjustedRadius);
        placedLabels.push({
          screenX: pos.x * cosR - pos.y * sinR,
          screenY: pos.x * sinR + pos.y * cosR,
          width: estimatedWidth,
          height: labelHeight,
          isReverse,
        });
        yOffset = testOffset;
      }

      offsets.set(data.index, yOffset);
    }

    return offsets;
  }, [visibleFeatures, zoom, rotation, featureStrokeWidth, shouldShowSequence]);

  return (
    <div
      ref={containerRef}
      style={{
        width: "100%",
        height,
        position: "relative",
        overflow: "hidden",
        backgroundColor: "#ffffff",
        borderRadius: "1rem",
        touchAction: "none",
        cursor: isDragging ? "crosshair" : "default",
      }}
    >
      {/* Title */}
      <div
        style={{
          position: "absolute",
          top: "1rem",
          left: "1rem",
          fontSize: "0.875rem",
          fontWeight: 500,
          color: COLORS.textDark,
          zIndex: 10,
        }}
      >
        {title}
      </div>

      {/* Zoom indicator */}
      <div
        style={{
          position: "absolute",
          top: "1rem",
          right: "1rem",
          fontSize: "0.75rem",
          color: COLORS.textDark,
          opacity: 0.7,
          zIndex: 10,
        }}
      >
        {zoom.toFixed(1)}x
      </div>

      {/* Selection indicator moved to the parent InteractiveSequenceViewer
          header (next to the Annotations count). The viewer still drives
          state via onSelectionChange. */}

      {/* Controls hint */}
      <div
        style={{
          position: "absolute",
          bottom: "0.5rem",
          left: "50%",
          transform: "translateX(-50%)",
          fontSize: "0.65rem",
          color: COLORS.textDark,
          opacity: 0.5,
          zIndex: 10,
          whiteSpace: "nowrap",
        }}
      >
        Drag/click: select • Shift+click: extend • Scroll up/down: rotate • {typeof navigator !== "undefined" && navigator.platform?.includes("Mac") ? "⌘" : "Ctrl"}+scroll: zoom
      </div>

      {/* Length indicator */}
      <div
        style={{
          position: "absolute",
          bottom: "0.5rem",
          right: "1rem",
          fontSize: "0.75rem",
          color: COLORS.textDark,
          opacity: 0.7,
          zIndex: 10,
        }}
      >
        {totalLength.toLocaleString()} bp
      </div>

      {/* SVG Plasmid */}
      <svg
        ref={svgRef}
        viewBox={`${-viewBoxSize / 2 + viewBoxOffset.x} ${-viewBoxSize / 2 + viewBoxOffset.y} ${viewBoxSize} ${viewBoxSize}`}
        style={{
          width: "100%",
          height: "100%",
        }}
        onMouseDown={handleMouseDown}
      >
        <g transform={`rotate(${(rotation * 180) / Math.PI})`}>
          {/* Selection arc - spans both forward and reverse feature areas.
              Offset each endpoint by -0.5 bp so the box edges sit on the
              boundary BETWEEN two nucleotides rather than through a letter's
              center. */}
          {selectionStart !== null && selectionEnd !== null && selectionLength > 0 && (
            <path
              d={createFilledArcPath(
                selectionStart - 0.5,
                selectionEnd - 0.5,
                totalLength,
                BASE_RADIUS - featureStrokeWidth - 5,
                BASE_RADIUS + featureStrokeWidth + 5
              )}
              fill={COLORS.selection}
              stroke={COLORS.selectionStroke}
              strokeWidth={Math.max(0.2, featureStrokeWidth * 0.035)}
              style={{ pointerEvents: "none" }}
            />
          )}

          {/* Backbone circle - thin black line at 30% opacity, always visible */}
          <path
            d={backbonePath}
            fill="none"
            stroke={COLORS.backbone}
            strokeWidth={BACKBONE_STROKE / zoom}
            opacity={BACKBONE_OPACITY}
          />

          {/* Tick marks — remain visible at all zooms; interval refines as zoom increases. */}
          {tickMarks.map((tick) => (
            <g key={`tick-${tick.bp}`}>
              <line
                x1={tick.inner.x}
                y1={tick.inner.y}
                x2={tick.outer.x}
                y2={tick.outer.y}
                stroke={COLORS.textDark}
                strokeWidth={1 / zoom}
                opacity={0.5}
              />
              <text
                x={tick.labelPos.x}
                y={tick.labelPos.y}
                fontSize={labelFontSize}
                fill={COLORS.textDark}
                fontWeight={500}
                textAnchor="middle"
                dominantBaseline="middle"
                opacity={0.5}
                transform={`rotate(${(-rotation * 180) / Math.PI}, ${tick.labelPos.x}, ${tick.labelPos.y})`}
              >
                {tick.bp === 0 ? "1" : tick.bp.toLocaleString()}
              </text>
            </g>
          ))}

          {/* Annotation arcs - rectangular for no direction, arrow for directional */}
          {annotationData.map((data) => {
            const { annotation, index, featureRadius, innerR, outerR, midAngle, color, isReverse } = data;

            // Determine if feature has direction
            const hasDirection = annotation.direction === 1 || annotation.direction === -1;

            // Create path based on direction (no hover expansion — arc radii stay constant)
            const pathD = hasDirection
              ? createArrowArcPath(
                  annotation.start,
                  annotation.end,
                  totalLength,
                  innerR,
                  outerR,
                  annotation.direction as number
                )
              : createFilledArcPath(
                  annotation.start,
                  annotation.end,
                  totalLength,
                  innerR,
                  outerR
                );

            return (
              <g key={`ann-${index}`} data-annotation>
                <path
                  d={pathD}
                  data-annotation-index={index}
                  fill={color}
                  fillOpacity={FEATURE_OPACITY}
                  stroke="#000000"
                  strokeWidth={Math.max(0.2, featureStrokeWidth * 0.035)}
                  strokeOpacity={1}
                  style={{
                    cursor: onAnnotationClick ? "pointer" : "default",
                    transition: "all 0.15s ease",
                  }}
                  onMouseEnter={() => setHoveredAnnotation(index)}
                  onMouseLeave={() => setHoveredAnnotation(null)}
                  onClick={(e) => {
                    e.stopPropagation();
                    handleAnnotationSelect(annotation, e);
                  }}
                  onDoubleClick={(e) => {
                    e.stopPropagation();
                    // Resolve from the set of feature paths actually under the
                    // cursor — not all features whose bp range covers this
                    // point — so overlapping features at different radii don't
                    // bleed into each other's click targets.
                    const hits = typeof document !== "undefined"
                      ? document.elementsFromPoint(e.clientX, e.clientY)
                      : [];
                    const hitIndices: number[] = [];
                    for (const el of hits) {
                      const raw = (el as Element).getAttribute?.("data-annotation-index");
                      if (raw != null) hitIndices.push(parseInt(raw, 10));
                    }
                    const candidates = hitIndices
                      .map((i) => annotationData[i])
                      .filter(Boolean);
                    if (candidates.length > 0) {
                      candidates.sort((a, b) =>
                        getAnnotationLength(a.annotation.start, a.annotation.end, totalLength) -
                        getAnnotationLength(b.annotation.start, b.annotation.end, totalLength)
                      );
                      onAnnotationClick?.(candidates[0].annotation, candidates[0].index);
                    } else {
                      onAnnotationClick?.(annotation, index);
                    }
                  }}
                />
              </g>
            );
          })}

          {/* Interaction chords — curved Bezier links connecting participants
              of each functional interaction (promoter→CDS→polyA, RBS→CDS,
              Pol3→sgRNA, recombination pairs, insulator boundaries). Rendered
              as a translucent curved path inside the backbone circle, color-
              coded by SBO term. A legend lives in the parent viewer. */}
          {showInteractions && interactions && interactions.length > 0 && (
            <g style={{ pointerEvents: "none" }}>
              {interactions.flatMap((ix, ixIdx) => {
                const INTERACTION_COLORS: Record<string, string> = {
                  // SBO term → hex
                  "http://identifiers.org/biomodels.sbo/SBO:0000589": "#2563eb", // genetic production — blue
                  "http://identifiers.org/biomodels.sbo/SBO:0000183": "#0ea5e9", // transcription — cyan
                  "http://identifiers.org/biomodels.sbo/SBO:0000184": "#16a34a", // translation — green
                  "http://identifiers.org/biomodels.sbo/SBO:0000170": "#10b981", // stimulation — emerald
                  "http://identifiers.org/biomodels.sbo/SBO:0000169": "#dc2626", // inhibition — red
                  "http://identifiers.org/biomodels.sbo/SBO:0000182": "#f59e0b", // recombination — amber
                  "http://identifiers.org/biomodels.sbo/SBO:0000178": "#f43f5e", // cleavage — rose
                };
                const color = INTERACTION_COLORS[ix.sbo_term || ""] || "#7c3aed";
                const positioned = (ix.participants || []).filter(
                  (p) => typeof p.start === "number" && typeof p.end === "number"
                );
                if (positioned.length < 2) return [];

                // Sort participants by role ordering for consistent arrowheads:
                // stimulator/inhibitor → template → modifier
                const roleOrder: Record<string, number> = {
                  stimulator: 0, inhibitor: 0, reactant: 0,
                  template: 1,
                  modifier: 2, product: 2,
                };
                const sorted = [...positioned].sort(
                  (a, b) => (roleOrder[(a.role || "").toLowerCase()] ?? 99)
                          - (roleOrder[(b.role || "").toLowerCase()] ?? 99)
                );

                // Radius for chord endpoints (just inside backbone)
                const chordR = BASE_RADIUS - 4;
                const centerBp = (p: InteractionParticipant) =>
                  (((p.start as number) + (p.end as number)) / 2);

                const elements: React.ReactElement[] = [];
                for (let i = 0; i < sorted.length - 1; i++) {
                  const a = sorted[i];
                  const b = sorted[i + 1];
                  const aBp = centerBp(a);
                  const bBp = centerBp(b);
                  const aAng = bpToAngle(aBp, totalLength);
                  const bAng = bpToAngle(bBp, totalLength);
                  const aPt = polarToCartesian(aAng, chordR);
                  const bPt = polarToCartesian(bAng, chordR);

                  // Quadratic Bezier through center biases chord inward —
                  // use a control point scaled toward origin for a gentle
                  // dished arc.
                  const midAng = (aAng + bAng) / 2;
                  // If angular distance > PI, rotate mid to short-arc side
                  const diff = Math.abs(aAng - bAng);
                  const adjMidAng = diff > Math.PI ? midAng + Math.PI : midAng;
                  const ctrlRadius = chordR * 0.25;
                  const ctrl = polarToCartesian(adjMidAng, ctrlRadius);

                  const isHL = highlightedInteractionId === ix.interaction_id;
                  const strokeOpacity = isHL ? 0.95 : 0.45;
                  const strokeW = isHL ? 2.25 : 1.25;

                  elements.push(
                    <path
                      key={`ix-${ixIdx}-${i}`}
                      d={`M ${aPt.x} ${aPt.y} Q ${ctrl.x} ${ctrl.y} ${bPt.x} ${bPt.y}`}
                      fill="none"
                      stroke={color}
                      strokeWidth={strokeW}
                      strokeOpacity={strokeOpacity}
                      strokeLinecap="round"
                    />
                  );

                  // small endpoint dot on b (the downstream participant) to
                  // communicate direction
                  elements.push(
                    <circle
                      key={`ix-${ixIdx}-${i}-dot`}
                      cx={bPt.x}
                      cy={bPt.y}
                      r={isHL ? 2.5 : 1.75}
                      fill={color}
                      fillOpacity={strokeOpacity}
                    />
                  );
                }
                return elements;
              })}
            </g>
          )}

          {/* Cut-profile glyphs for cloning features with a cut_profile.
              Draws an EcoRI-style staggered bracket at the DNA strand level so
              the cut visualization sits on the DNA it cuts, not out past the
              feature arc:
                - forward-strand tick at the forward DNA radius (cut_top)
                - reverse-strand tick at the reverse DNA radius (cut_bottom)
                - arc connecting the two ticks through the backbone mid-line
              Overhang direction is encoded by which strand cuts first (5' vs 3'). */}
          {annotationData.map((data) => {
            const ann: any = data.annotation;
            if (ann?.layer !== "cloning_feature") return null;
            const cp = ann?.cut_profile;
            if (!cp || typeof cp.cut_top !== "number" || typeof cp.cut_bottom !== "number") return null;
            if (ann?.feature_family === "primer_design_warning") return null; // warnings have no cut profile

            // Cuts happen on inter-base boundaries. Forward-strand cut_top
            // shifts CCW by half a base; reverse-strand cut_bottom shifts CW
            // by half a base (antiparallel convention) so each tick lands
            // between the two bases it cleaves.
            const topAngle = bpToAngle(cp.cut_top - 0.5, totalLength);
            const botAngle = bpToAngle(cp.cut_bottom + 0.5, totalLength);
            // Anchor cut ticks to the DNA strand radii (where the bases sit),
            // so the cut glyph visualizes the cut on the DNA it belongs to.
            const forwardDnaR = BASE_RADIUS + featureStrokeWidth / 2;
            const reverseDnaR = BASE_RADIUS - featureStrokeWidth / 2;
            const tickSpan = featureStrokeWidth * 0.35; // short radial tick, scales with feature height
            const topOut = polarToCartesian(topAngle, forwardDnaR + tickSpan);
            const topIn  = polarToCartesian(topAngle, forwardDnaR - tickSpan);
            const botOut = polarToCartesian(botAngle, reverseDnaR + tickSpan);
            const botIn  = polarToCartesian(botAngle, reverseDnaR - tickSpan);

            // Overhang connector along the backbone midline between the two cuts.
            const midR = BASE_RADIUS;
            const topMid = polarToCartesian(topAngle, midR);
            const botMid = polarToCartesian(botAngle, midR);
            const lowAngle = Math.min(topAngle, botAngle);
            const highAngle = Math.max(topAngle, botAngle);
            const sweep = highAngle - lowAngle;
            const largeArc = sweep > Math.PI ? 1 : 0;
            const arcPath = `M ${topMid.x} ${topMid.y} A ${midR} ${midR} 0 ${largeArc} 1 ${botMid.x} ${botMid.y}`;
            const cutColor = ann?.feature_family === "gateway_att" ? "#6A1B9A" : "#C2185B";
            // 4x thinner than before (1.75 → ~0.44).
            const cutStroke = 0.44;

            return (
              <g key={`cut-${data.index}`} style={{ pointerEvents: "none" }}>
                {/* Forward-strand cut tick at the forward DNA radius */}
                <line
                  x1={topOut.x} y1={topOut.y} x2={topIn.x} y2={topIn.y}
                  stroke={cutColor} strokeWidth={cutStroke} strokeLinecap="round"
                />
                {/* Reverse-strand cut tick at the reverse DNA radius */}
                <line
                  x1={botOut.x} y1={botOut.y} x2={botIn.x} y2={botIn.y}
                  stroke={cutColor} strokeWidth={cutStroke} strokeLinecap="round"
                />
                {/* Overhang connector between the two cut positions */}
                {cp.cut_top !== cp.cut_bottom && (
                  <path
                    d={arcPath}
                    fill="none"
                    stroke={cutColor}
                    strokeWidth={cutStroke}
                    strokeLinecap="round"
                  />
                )}
              </g>
            );
          })}

          {/* Labels at low zoom - show when not using visibleFeatures labels.
              Below 2x zoom we only render labels that fit fully inside their
              feature arc; labels that would otherwise render outside the
              feature are hidden until the user zooms in to ≥2x (where
              visibleFeatures takes over and renders them with collision
              avoidance). This keeps the low-zoom view uncluttered. */}
          {zoom < 2 && labelPlacements.map((placement) => {
            const data = annotationData[placement.index];
            if (!data) return null;

            const { annotation, midAngle, spansOrigin, isReverse, featureRadius } = data;
            const isHovered = hoveredAnnotation === placement.index;
            const labelText = annotation.name + (spansOrigin ? " ⟳" : "");

            // Below 8x, suppress outside leader-line labels here too —
            // the visibleFeatures.map path already enforces this for >=2x, and
            // this branch covers <2x where visibleFeatures is empty.
            if (!placement.isInside) {
              return null;
            }

            if (placement.isInside) {
              // Use curved text for inside labels
              const pathId = `label-path-${placement.index}`;

              return (
                <g key={`label-${placement.index}`}>
                  <defs>
                    <path
                      id={pathId}
                      d={createTextArcPath(
                        annotation.start,
                        annotation.end,
                        totalLength,
                        featureRadius
                      )}
                      fill="none"
                    />
                  </defs>
                  <text
                    fontSize={labelFontSize}
                    fill={COLORS.textDark}
                    fontWeight={isHovered ? 600 : 500}
                    style={{ pointerEvents: "none" }}
                  >
                    <textPath
                      href={`#${pathId}`}
                      startOffset="50%"
                      textAnchor="middle"
                    >
                      {labelText}
                    </textPath>
                  </text>
                </g>
              );
            } else {
              const featureEdgeRadius = featureRadius + (isReverse ? -featureStrokeWidth / 2 - 2 : featureStrokeWidth / 2 + 2);
              const featureEdge = polarToCartesian(midAngle, featureEdgeRadius);
              const labelPos = polarToCartesian(midAngle, placement.labelRadius);
              const adjustedLabelY = labelPos.y + placement.yOffset;

              return (
                <g key={`label-${placement.index}`}>
                  <line
                    x1={featureEdge.x}
                    y1={featureEdge.y}
                    x2={labelPos.x}
                    y2={adjustedLabelY}
                    stroke={COLORS.labelLine}
                    strokeWidth={0.5}
                    opacity={0.4}
                  />
                  <text
                    x={labelPos.x}
                    y={adjustedLabelY}
                    fontSize={labelFontSize}
                    fill={COLORS.textDark}
                    fontWeight={isHovered ? 600 : 500}
                    textAnchor="middle"
                    dominantBaseline="middle"
                    transform={`rotate(${(-rotation * 180) / Math.PI}, ${labelPos.x}, ${adjustedLabelY})`}
                    style={{ pointerEvents: "none" }}
                  >
                    {labelText}
                  </text>
                </g>
              );
            }
          })}

          {/* DNA Sequence - forward and complementary strands */}
          {/* Positioned within the base layer feature area (same level as layer 0 features) */}
          {shouldShowSequence && (
            <g>
              {visibleBases.map(({ bp, base, complement, angle }) => {
                // Position DNA within the base layer feature margins
                // Forward strand: above the backbone line, within forward feature area
                // Reverse strand: below the backbone line, within reverse feature area
                const forwardRadius = BASE_RADIUS + featureStrokeWidth / 2;
                const reverseRadius = BASE_RADIUS - featureStrokeWidth / 2;

                const forwardPos = polarToCartesian(angle, forwardRadius);
                const complementPos = polarToCartesian(angle, reverseRadius);

                return (
                  <g key={`base-${bp}`}>
                    {/* Forward strand (5' → 3') - above backbone, angled inward toward center */}
                    <text
                      x={forwardPos.x}
                      y={forwardPos.y}
                      fontSize={sequenceFontSize}
                      fill={BASE_COLORS[base] || COLORS.textDark}
                      fontFamily="'Droid Sans Mono', 'Courier New', Courier, monospace"
                      fontWeight={500}
                      textAnchor="middle"
                      dominantBaseline="middle"
                      transform={`rotate(${(angle * 180) / Math.PI + 90}, ${forwardPos.x}, ${forwardPos.y})`}
                    >
                      {base}
                    </text>
                    {/* Complementary strand (3' → 5') - below backbone, angled inward toward center */}
                    <text
                      x={complementPos.x}
                      y={complementPos.y}
                      fontSize={sequenceFontSize}
                      fill={BASE_COLORS[complement] || COLORS.textDark}
                      fontFamily="'Droid Sans Mono', 'Courier New', Courier, monospace"
                      fontWeight={500}
                      textAnchor="middle"
                      dominantBaseline="middle"
                      opacity={0.7}
                      transform={`rotate(${(angle * 180) / Math.PI + 90}, ${complementPos.x}, ${complementPos.y})`}
                    >
                      {complement}
                    </text>
                  </g>
                );
              })}
            </g>
          )}

          {/* Feature labels - curved text inside when no DNA, outside when DNA visible */}
          {visibleFeatures.map((data) => {
            const { annotation, midAngle, isReverse, color, layer, featureRadius, innerR, outerR, labelFitsInside, arcLength } = data;

            // Font size scales down as zoom increases to keep labels from becoming huge
            const zoomedFontSize = Math.max(1.5, Math.min(12 / zoom, featureStrokeWidth * 0.6));

            // When DNA is visible (>=16x), always show labels outside.
            // When no DNA, show inside if the label fits in the feature arc.
            // Outside labels with leader lines are visually noisy — suppress
            // them until the user is zoomed in far enough (>=8x) for the
            // remaining placements to breathe.
            const showOutside = shouldShowSequence || (!labelFitsInside && zoom >= 8);

            // Below 8x, hide labels that do not fit inside their feature
            // rather than rendering a chaotic outside leader line.
            if (!labelFitsInside && !showOutside) {
              return null;
            }

            // Get collision-resolved offset for this label
            const labelYOffset = visibleFeatureLabelOffsets.get(data.index) || 0;

            if (showOutside) {
              // Anchor the label to the feature's midline so the name sits
              // inside the feature arc rather than floating above/below it.
              // Collision-resolution still stacks via labelYOffset if multiple
              // short labels land on the same slot.
              const baseRadius = featureRadius;
              const adjustedRadius = baseRadius + labelYOffset;
              const rawLabelPos = polarToCartesian(midAngle, adjustedRadius);

              // Clamp the label into the visible viewBox so in-frame features
              // keep their full name on screen even when the label's natural
              // position would be cut off by the viewport edge. Clamping is
              // done in screen-space (after the outer <g> rotation) then
              // inverse-rotated back into the pre-rotation frame.
              const cosR = Math.cos(rotation);
              const sinR = Math.sin(rotation);
              const sxRaw = rawLabelPos.x * cosR - rawLabelPos.y * sinR;
              const syRaw = rawLabelPos.x * sinR + rawLabelPos.y * cosR;
              const vXmin = viewBoxOffset.x - viewBoxSize / 2;
              const vXmax = viewBoxOffset.x + viewBoxSize / 2;
              const vYmin = viewBoxOffset.y - viewBoxSize / 2;
              const vYmax = viewBoxOffset.y + viewBoxSize / 2;
              const estimatedLabelWidth = annotation.name.length * zoomedFontSize * 0.6;
              const padX = estimatedLabelWidth / 2 + zoomedFontSize * 0.5;
              const padY = zoomedFontSize * 0.8;
              const sxC = Math.max(vXmin + padX, Math.min(vXmax - padX, sxRaw));
              const syC = Math.max(vYmin + padY, Math.min(vYmax - padY, syRaw));
              const labelPos = {
                x: sxC * Math.cos(-rotation) - syC * Math.sin(-rotation),
                y: sxC * Math.sin(-rotation) + syC * Math.cos(-rotation),
              };

              // Leader line from feature edge to the resolved label position.
              // When the label has been pushed away by collision-resolution this
              // gives the user a visible connection back to its feature.
              const featureEdgeR = featureRadius + (isReverse ? -featureStrokeWidth / 2 - 1 : featureStrokeWidth / 2 + 1);
              const featureEdge = polarToCartesian(midAngle, featureEdgeR);

              return (
                <g key={`feature-label-${data.index}`}>
                  <line
                    x1={featureEdge.x}
                    y1={featureEdge.y}
                    x2={labelPos.x}
                    y2={labelPos.y}
                    stroke={COLORS.labelLine}
                    strokeWidth={Math.max(0.3, 0.6 / zoom)}
                    opacity={0.55}
                  />
                  {/* Background for readability */}
                  <text
                    x={labelPos.x}
                    y={labelPos.y}
                    fontSize={zoomedFontSize}
                    fill="white"
                    stroke="white"
                    strokeWidth={2 / zoom}
                    fontWeight={600}
                    textAnchor="middle"
                    dominantBaseline="middle"
                    transform={`rotate(${(-rotation * 180) / Math.PI}, ${labelPos.x}, ${labelPos.y})`}
                    style={{ pointerEvents: "none" }}
                  >
                    {annotation.name}
                  </text>
                  {/* Label text */}
                  <text
                    x={labelPos.x}
                    y={labelPos.y}
                    fontSize={zoomedFontSize}
                    fill={color}
                    fontWeight={600}
                    textAnchor="middle"
                    dominantBaseline="middle"
                    transform={`rotate(${(-rotation * 180) / Math.PI}, ${labelPos.x}, ${labelPos.y})`}
                    style={{ pointerEvents: "none" }}
                  >
                    {annotation.name}
                  </text>
                </g>
              );
            } else {
              // Inside label using textPath for curved text along the arc
              const pathId = `text-path-${data.index}`;
              const textPathRadius = featureRadius;

              return (
                <g key={`feature-label-${data.index}`}>
                  <defs>
                    <path
                      id={pathId}
                      d={createTextArcPath(
                        annotation.start,
                        annotation.end,
                        totalLength,
                        textPathRadius
                      )}
                      fill="none"
                    />
                  </defs>
                  {/* Background stroke for readability */}
                  <text
                    fontSize={zoomedFontSize}
                    fill="white"
                    stroke="white"
                    strokeWidth={2 / zoom}
                    fontWeight={600}
                    style={{ pointerEvents: "none" }}
                  >
                    <textPath
                      href={`#${pathId}`}
                      startOffset="50%"
                      textAnchor="middle"
                    >
                      {annotation.name}
                    </textPath>
                  </text>
                  {/* Main text */}
                  <text
                    fontSize={zoomedFontSize}
                    fill={COLORS.textDark}
                    fontWeight={600}
                    style={{ pointerEvents: "none" }}
                  >
                    <textPath
                      href={`#${pathId}`}
                      startOffset="50%"
                      textAnchor="middle"
                    >
                      {annotation.name}
                    </textPath>
                  </text>
                </g>
              );
            }
          })}
        </g>
      </svg>
    </div>
  );
}
