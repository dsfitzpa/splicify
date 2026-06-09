import { NextResponse } from "next/server";

// Eager annotate-on-upload proxy. Forwards a multipart .gb upload to the VPS
// /agent_v2/annotate-on-upload, which runs the v1 annotation pipeline
// (annotate_cached depth=full) and returns a `viz` payload identical to
// what the chat envelope emits. No agent involvement, no LLM cost.
// Lets the frontend paint the circular viewer immediately, before the chat call.
export const maxDuration = 60;

export async function POST(req: Request) {
  let AGENT_URL = process.env.BACKEND_AGENT_V2_URL || "";
  if (!AGENT_URL) {
    const base = process.env.BACKEND_API_URL || process.env.BACKEND_URL || "";
    if (base) AGENT_URL = base.replace(/\/api\/chat\/?$/, "/agent_v2/chat");
  }
  if (!AGENT_URL) {
    return NextResponse.json(
      { ok: false, error: "Missing env var BACKEND_AGENT_V2_URL" },
      { status: 500 }
    );
  }
  const ANNOTATE_URL = AGENT_URL.replace(/\/chat\/?$/, "/annotate-on-upload");

  const outForm = new FormData();
  let incoming: FormData;
  try {
    incoming = await req.formData();
  } catch {
    return NextResponse.json(
      { ok: false, error: "Invalid multipart/form-data body." },
      { status: 400 }
    );
  }
  for (const [k, v] of incoming.entries()) outForm.append(k, v);

  let upstreamRes: Response;
  try {
    upstreamRes = await fetch(ANNOTATE_URL, {
      method: "POST",
      body: outForm,
      signal: AbortSignal.timeout(55_000),
    });
  } catch (err: unknown) {
    const msg = err instanceof Error ? err.message : String(err);
    console.error(`[agent_v2/annotate-on-upload] upstream fetch failed:`, msg);
    return NextResponse.json(
      { ok: false, error: `Upstream fetch failed: ${msg}` },
      { status: 502 }
    );
  }

  try {
    const data = await upstreamRes.json();
    return NextResponse.json(data, { status: upstreamRes.status });
  } catch {
    const raw = await upstreamRes.text().catch(() => "");
    return NextResponse.json(
      { ok: false, error: `Upstream non-JSON: ${raw.slice(0, 200)}` },
      { status: 502 }
    );
  }
}
