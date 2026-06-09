import { NextResponse } from "next/server";

// Vercel proxy for the AI Agent (v2) — forwards multipart/JSON to the VPS
// /agent_v2/chat (Claude Sonnet 4.6 tool-use loop with AIPlasmidDesign tools).
// Up to 5 min: full-pipeline runs (triage + 3 Explore + Plan + Main + Summarizer)
// take 3-5 min on Sonnet 4.6 across ~10-20 tool calls.
export const maxDuration = 300;

export async function POST(req: Request) {
  // Resolve backend URL. Prefer explicit BACKEND_AGENT_V2_URL, fall back to
  // deriving from BACKEND_API_URL by replacing /api/chat -> /agent_v2/chat.
  let AGENT_URL = process.env.BACKEND_AGENT_V2_URL || "";
  if (!AGENT_URL) {
    const base = process.env.BACKEND_API_URL || process.env.BACKEND_URL || "";
    if (base) AGENT_URL = base.replace(/\/api\/chat\/?$/, "/agent_v2/chat");
  }
  if (!AGENT_URL) {
    console.error("[agent_v2/chat] No backend URL. Set BACKEND_AGENT_V2_URL or BACKEND_API_URL.");
    return NextResponse.json(
      {
        ok: false,
        reply: "Server configuration error: BACKEND_AGENT_V2_URL is not set.",
        error: "Missing env var BACKEND_AGENT_V2_URL",
      },
      { status: 500 }
    );
  }

  const outForm = new FormData();
  const contentTypeIn = req.headers.get("content-type") || "";

  if (contentTypeIn.includes("multipart/form-data")) {
    let incomingForm: FormData;
    try {
      incomingForm = await req.formData();
    } catch {
      return NextResponse.json(
        { ok: false, reply: "Invalid multipart/form-data body.", error: "Bad request body" },
        { status: 400 }
      );
    }
    for (const [k, v] of incomingForm.entries()) outForm.append(k, v);
  } else {
    let payload: Record<string, unknown>;
    try {
      payload = await req.json();
    } catch {
      return NextResponse.json(
        { ok: false, reply: "Invalid JSON body.", error: "Bad request body" },
        { status: 400 }
      );
    }
    for (const [k, v] of Object.entries(payload ?? {})) {
      if (v !== null && v !== undefined) outForm.append(k, String(v));
    }
  }

  let upstreamRes: Response;
  try {
    upstreamRes = await fetch(AGENT_URL, {
      method: "POST",
      body: outForm,
      // 290 s upstream cap — leaves a small margin under Vercel's 300 s limit.
      signal: AbortSignal.timeout(290_000),
    });
  } catch (err: unknown) {
    const msg = err instanceof Error ? err.message : String(err);
    console.error(`[agent_v2/chat] Failed to reach backend at ${AGENT_URL}:`, msg);
    return NextResponse.json(
      {
        ok: false,
        reply: `Could not reach the agent backend. Please try again in a moment. (${msg})`,
        error: `Upstream fetch failed: ${msg}`,
      },
      { status: 502 }
    );
  }

  const contentTypeOut = upstreamRes.headers.get("content-type") || "";
  if (contentTypeOut.includes("application/json")) {
    try {
      const data = await upstreamRes.json();
      return NextResponse.json(data, { status: upstreamRes.status });
    } catch {
      const raw = await upstreamRes.text().catch(() => "");
      console.error("[agent_v2/chat] Backend returned non-JSON despite content-type:", raw.slice(0, 200));
      return NextResponse.json(
        { ok: false, reply: "Backend returned malformed JSON.", error: "Upstream non-JSON" },
        { status: 502 }
      );
    }
  }

  const raw = await upstreamRes.text().catch(() => "");
  return NextResponse.json(
    {
      ok: false,
      reply: raw.slice(0, 1000) || "Empty backend response.",
      error: "Upstream returned non-JSON",
    },
    { status: upstreamRes.status }
  );
}
