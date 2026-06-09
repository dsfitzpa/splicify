import { NextResponse } from "next/server";

function deriveBackendBase(rawUrl: string): string {
  const parsed = new URL(rawUrl);
  const basePath = parsed.pathname.replace(/\/chat\/?$/, "").replace(/\/$/, "");
  return `${parsed.origin}${basePath}`;
}

function originOnly(rawUrl: string): string {
  return new URL(rawUrl).origin;
}

export async function POST(req: Request) {
  const BACKEND_URL = process.env.BACKEND_API_URL || process.env.BACKEND_URL;

  if (!BACKEND_URL) {
    return NextResponse.json(
      { ok: false, error: "Missing env var BACKEND_API_URL or BACKEND_URL" },
      { status: 500 }
    );
  }

  const backendBase = deriveBackendBase(BACKEND_URL);
  const backendOrigin = originOnly(BACKEND_URL);

  let body: unknown;
  try {
    body = await req.json();
  } catch {
    return NextResponse.json({ ok: false, error: "Invalid JSON body" }, { status: 400 });
  }

  // Try the new analyze_plasmid endpoint first, fall back to analyze_intent
  const candidates = [
    `${backendBase}/plannotate/analyze_plasmid`,
    `${backendOrigin}/plannotate/analyze_plasmid`,
    `${backendBase}/plannotate/analyze_intent`,
    `${backendOrigin}/plannotate/analyze_intent`,
  ];

  let lastRaw = "";
  let lastStatus = 502;
  for (const analyzeUrl of candidates) {
    let upstreamRes: Response;
    try {
      upstreamRes = await fetch(analyzeUrl, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(body),
      });
    } catch {
      continue;
    }

    const raw = await upstreamRes.text();
    try {
      const data = raw ? JSON.parse(raw) : null;
      if (upstreamRes.status === 404 || data?.detail === "Not Found") {
        lastStatus = upstreamRes.status;
        continue;
      }
      return NextResponse.json(data, { status: upstreamRes.status });
    } catch {
      lastRaw = raw.slice(0, 500);
      lastStatus = upstreamRes.status || 502;
    }
  }

  return NextResponse.json(
    { ok: false, error: "Backend returned non-JSON", raw: lastRaw || undefined },
    { status: lastStatus }
  );
}
