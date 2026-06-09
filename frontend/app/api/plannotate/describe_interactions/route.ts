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

  const urls = [
    `${backendBase}/plannotate/describe_interactions`,
    `${backendOrigin}/plannotate/describe_interactions`,
  ];

  let lastErr = "No backend reachable";
  let lastStatus = 502;
  for (const url of urls) {
    try {
      const upstream = await fetch(url, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(body),
      });
      const raw = await upstream.text();
      let data: any = null;
      try {
        data = raw ? JSON.parse(raw) : null;
      } catch {
        lastErr = "Non-JSON response";
        lastStatus = upstream.status || 502;
        continue;
      }
      if (upstream.ok) {
        return NextResponse.json(data, { status: upstream.status });
      }
      lastErr = (data && (data.error || data.detail)) || "Upstream error";
      lastStatus = upstream.status;
    } catch (err) {
      lastErr = `fetch failed: ${err}`;
      continue;
    }
  }
  return NextResponse.json(
    { ok: false, error: lastErr, summary: "", bullets: [], markdown: "" },
    { status: lastStatus }
  );
}
