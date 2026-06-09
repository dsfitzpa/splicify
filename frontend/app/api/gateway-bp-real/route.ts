import { NextResponse } from "next/server";

export const dynamic = "force-dynamic";

export async function GET() {
  // Extract base URL from BACKEND_API_URL (remove /api/chat if present)
  const backendUrl = process.env.BACKEND_API_URL || "http://localhost:8000";
  const baseUrl = backendUrl.replace(/\/api\/chat$/, '');

  const candidateUrls = [
    `${baseUrl}/cloning/demo/gateway_bp_real`,
    "http://localhost:8000/cloning/demo/gateway_bp_real",
    "http://127.0.0.1:8000/cloning/demo/gateway_bp_real",
  ];

  let lastError: any = null;

  for (const url of candidateUrls) {
    try {
      const controller = new AbortController();
      const timeoutId = setTimeout(() => controller.abort(), 15000);

      const res = await fetch(url, {
        method: "GET",
        signal: controller.signal,
        headers: {
          "Content-Type": "application/json",
        },
      });

      clearTimeout(timeoutId);

      const data = await res.json();

      if (!res.ok) {
        lastError = data;
        continue;
      }

      return NextResponse.json(data);
    } catch (err: any) {
      lastError = err;
      continue;
    }
  }

  return NextResponse.json(
    {
      error: "Failed to reach backend Gateway BP Real demo endpoint",
      details: lastError?.message || String(lastError),
    },
    { status: 500 }
  );
}
