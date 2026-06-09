import { NextResponse } from "next/server";
import type { NextRequest } from "next/server";

// Hosts that should serve the agent-first /splicify experience instead of
// the engine UI at the project root.
const SPLICIFY_HOSTS = new Set([
  "splicify.ai",
  "www.splicify.ai",
]);


export function middleware(request: NextRequest) {
  const host = (request.headers.get("host") || "").toLowerCase().split(":")[0];
  if (!SPLICIFY_HOSTS.has(host)) {
    return NextResponse.next();
  }

  const url = request.nextUrl;
  // Do not rewrite API routes, Next.js internals, or files that already live
  // under /splicify (otherwise the rewrite loops).
  if (
    url.pathname.startsWith("/api") ||
    url.pathname.startsWith("/_next") ||
    url.pathname.startsWith("/splicify") ||
    url.pathname === "/favicon.ico"
  ) {
    return NextResponse.next();
  }

  // Every other path under splicify.ai serves the focused agent page.
  const rewritten = url.clone();
  rewritten.pathname = "/splicify";
  return NextResponse.rewrite(rewritten);
}


export const config = {
  // Skip static assets but otherwise inspect every request.
  matcher: ["/((?!_next/static|_next/image|favicon.ico|.*\\.(?:png|jpg|jpeg|svg|gif|webp|ico|css|js|woff2?)).*)"],
};
