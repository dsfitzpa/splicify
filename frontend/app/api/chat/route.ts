import { NextResponse } from "next/server";

// Allow up to 60 s — Gibson design on large fragments can take 8-12 s on the VPS.
// Without this Vercel's default (10 s on Hobby, 15 s on Pro) silently kills the request.
export const maxDuration = 60;

export async function POST(req: Request) {
  // Support both new BACKEND_API_URL and legacy BACKEND_URL / N8N_WEBHOOK_URL
  const BACKEND_URL =
    process.env.BACKEND_API_URL ||
    process.env.BACKEND_URL ||
    process.env.N8N_WEBHOOK_URL;

  if (!BACKEND_URL) {
    console.error("[chat/route] No backend URL configured. Set BACKEND_API_URL in Vercel environment variables.");
    return NextResponse.json(
      {
        ok: false,
        reply: "Server configuration error: BACKEND_API_URL is not set. Please contact the administrator.",
        error: "Missing env var BACKEND_API_URL",
      },
      { status: 500 }
    );
  }

  // Build a FormData body to forward.
  // The Python backend (FastAPI) only accepts multipart/form-data for Form() fields.
  // Whether the client sent JSON or multipart, we always forward as FormData.
  const outForm = new FormData();
  const contentTypeIn = req.headers.get("content-type") || "";

  if (contentTypeIn.includes("multipart/form-data")) {
    // CASE 1: client sent multipart (file uploads) — forward fields as-is
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
    // CASE 2: client sent JSON — convert each field to a FormData string entry
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

  // Forward to the Python backend with a generous timeout
  let upstreamRes: Response;
  try {
    upstreamRes = await fetch(BACKEND_URL, {
      method: "POST",
      body: outForm,
      // Node 18+ fetch supports AbortController; give the VPS 55 s to respond
      signal: AbortSignal.timeout(55_000),
    });
  } catch (err: unknown) {
    const msg = err instanceof Error ? err.message : String(err);
    console.error(`[chat/route] Failed to reach backend at ${BACKEND_URL}:`, msg);
    return NextResponse.json(
      {
        ok: false,
        reply: `Could not reach the design backend. Please try again in a moment. (${msg})`,
        error: `Upstream fetch failed: ${msg}`,
      },
      { status: 502 }
    );
  }

  // Parse and forward the JSON response
  const contentTypeOut = upstreamRes.headers.get("content-type") || "";
  if (contentTypeOut.includes("application/json")) {
    try {
      const data = await upstreamRes.json();
      return NextResponse.json(data, { status: upstreamRes.status });
    } catch {
      const raw = await upstreamRes.text().catch(() => "");
      console.error("[chat/route] Backend returned non-JSON despite content-type:", raw.slice(0, 200));
      return NextResponse.json(
        { ok: false, reply: "Backend returned an unreadable response.", raw },
        { status: 502 }
      );
    }
  }

  // Fallback: text body that might still be JSON
  const raw = await upstreamRes.text().catch(() => "");
  try {
    const parsed = JSON.parse(raw);
    return NextResponse.json(parsed, { status: upstreamRes.status });
  } catch {
    console.error("[chat/route] Backend returned non-JSON text:", raw.slice(0, 200));
    return NextResponse.json(
      { ok: false, reply: "Backend returned an unexpected response format.", raw },
      { status: 502 }
    );
  }
}
