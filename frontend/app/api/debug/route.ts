import { NextResponse } from "next/server";

/**
 * GET /api/debug
 * Returns backend connectivity info. Use this to verify BACKEND_API_URL is wired up on Vercel.
 * Does NOT expose the actual URL value for security — just whether it is set and reachable.
 */
export async function GET() {
  const BACKEND_URL =
    process.env.BACKEND_API_URL ||
    process.env.BACKEND_URL ||
    process.env.N8N_WEBHOOK_URL;

  if (!BACKEND_URL) {
    return NextResponse.json({
      backend_configured: false,
      error: "BACKEND_API_URL env var is not set on this deployment.",
    });
  }

  // Probe the backend with a minimal ping (empty message → help text is fine, we just want 200)
  let reachable = false;
  let status: number | null = null;
  let probe_error: string | null = null;

  try {
    const form = new FormData();
    form.append("message", "ping");
    form.append("include_ai_explanation", "false");

    const res = await fetch(BACKEND_URL, {
      method: "POST",
      body: form,
      signal: AbortSignal.timeout(10_000),
    });
    status = res.status;
    reachable = res.ok || res.status < 500;
  } catch (err: unknown) {
    probe_error = err instanceof Error ? err.message : String(err);
  }

  return NextResponse.json({
    backend_configured: true,
    // Show only the hostname so the user can confirm which VPS is being used
    backend_host: (() => { try { return new URL(BACKEND_URL).host; } catch { return "invalid URL"; } })(),
    reachable,
    http_status: status,
    probe_error,
  });
}
