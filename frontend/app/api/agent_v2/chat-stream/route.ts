import { NextResponse } from "next/server";

// SSE proxy for the AI Agent (v2). Forwards JSON body to the VPS
// /agent_v2/chat-stream and streams the upstream text/event-stream body
// back to the client unchanged. The first SSE event is `event: shorthand`
// (the triage classifier's prompt summary), followed by `event: envelope`
// with the full chat envelope.
// Edge Runtime: streaming responses are not subject to function-duration
// caps on any Vercel tier. The agent v2 pipeline can take 2-5 minutes
// on Sonnet 4.6; the SSE stream needs to stay alive that long. Node
// runtime functions (even with maxDuration=300) get killed early on
// Hobby tier and may buffer SSE responses on serverless wrappers.
export const runtime = "edge";

export async function POST(req: Request) {
  let AGENT_URL = process.env.BACKEND_AGENT_V2_URL || "";
  if (!AGENT_URL) {
    const base = process.env.BACKEND_API_URL || process.env.BACKEND_URL || "";
    if (base) AGENT_URL = base.replace(/\/api\/chat\/?$/, "/agent_v2/chat");
  }
  if (!AGENT_URL) {
    return NextResponse.json(
      {
        ok: false,
        reply: "Server configuration error: BACKEND_AGENT_V2_URL is not set.",
        error: "Missing env var BACKEND_AGENT_V2_URL",
      },
      { status: 500 }
    );
  }
  // Swap /chat suffix for /chat-stream
  const STREAM_URL = AGENT_URL.replace(/\/chat\/?$/, "/chat-stream");

  // Forward the body as raw text (the upstream accepts JSON or multipart).
  const contentType = req.headers.get("content-type") || "application/json";
  const bodyBuf = await req.arrayBuffer();

  let upstreamRes: Response;
  try {
    upstreamRes = await fetch(STREAM_URL, {
      method: "POST",
      headers: { "content-type": contentType },
      body: bodyBuf,
      signal: AbortSignal.timeout(290_000),
    });
  } catch (err: unknown) {
    const msg = err instanceof Error ? err.message : String(err);
    console.error(`[agent_v2/chat-stream] upstream fetch failed:`, msg);
    return NextResponse.json(
      {
        ok: false,
        reply: `Could not reach the agent backend. (${msg})`,
        error: `Upstream fetch failed: ${msg}`,
      },
      { status: 502 }
    );
  }

  // Stream upstream body verbatim back to the client.
  return new Response(upstreamRes.body, {
    status: upstreamRes.status,
    headers: {
      "content-type": "text/event-stream",
      "cache-control": "no-cache",
      "x-accel-buffering": "no",
    },
  });
}
