import { NextResponse } from "next/server";

function deriveBackendBase(rawUrl: string): string {
  const parsed = new URL(rawUrl);
  const basePath = parsed.pathname.replace(/\/chat\/?$/, "").replace(/\/$/, "");
  return `${parsed.origin}${basePath}`;
}

function originOnly(rawUrl: string): string {
  return new URL(rawUrl).origin;
}

function extractErrorMessage(data: any): string {
  if (!data) return "Unknown backend error";
  if (typeof data.error === "string" && data.error) return data.error;
  if (typeof data.detail === "string" && data.detail) return data.detail;
  if (typeof data.message === "string" && data.message) return data.message;
  if (data.detail && typeof data.detail === "object") {
    if (typeof data.detail.error === "string" && data.detail.error) return data.detail.error;
    if (typeof data.detail.details === "string" && data.detail.details) return data.detail.details;
  }
  return "Guide-design request failed";
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

  async function postTo(url: string): Promise<Response> {
    return fetch(url, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(body),
    });
  }

  const candidates = [
    {
      url: `${backendBase}/plannotate/design_guides`,
      proxyTarget: "design_guides",
      backendRef: backendBase,
    },
    {
      url: `${backendOrigin}/plannotate/design_guides`,
      proxyTarget: "design_guides",
      backendRef: backendOrigin,
    },
  ];

  let lastStatus = 502;
  let lastError = "Backend returned non-JSON";
  let lastRaw = "";

  for (const candidate of candidates) {
    let upstreamRes: Response;
    try {
      upstreamRes = await postTo(candidate.url);
    } catch (err) {
      lastError = `Could not reach backend: ${err}`;
      continue;
    }

    const raw = await upstreamRes.text();
    let data: any = null;
    try {
      data = raw ? JSON.parse(raw) : null;
    } catch {
      lastStatus = upstreamRes.status || 502;
      lastError = "Backend returned non-JSON";
      lastRaw = raw.slice(0, 500);
      continue;
    }

    const isNotFoundLike =
      upstreamRes.status === 404 ||
      upstreamRes.status === 405 ||
      data?.detail === "Not Found";

    if (isNotFoundLike) {
      lastStatus = upstreamRes.status || 404;
      lastError = extractErrorMessage(data);
      continue;
    }

    if (!upstreamRes.ok) {
      return NextResponse.json(
        {
          ok: false,
          error: extractErrorMessage(data),
          details: data,
          proxy_target: candidate.proxyTarget,
          backend_base: candidate.backendRef,
        },
        { status: upstreamRes.status || 500 }
      );
    }

    return NextResponse.json(
      {
        ...data,
        proxy_target: candidate.proxyTarget,
        backend_base: candidate.backendRef,
      },
      { status: upstreamRes.status }
    );
  }

  return NextResponse.json(
    {
      ok: false,
      error: lastError,
      backend_base: backendBase,
      backend_origin: backendOrigin,
      raw: lastRaw || undefined,
    },
    { status: lastStatus || 502 }
  );
}
