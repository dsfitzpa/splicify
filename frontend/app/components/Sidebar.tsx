"use client";

import React from "react";
import Link from "next/link";
import { usePathname } from "next/navigation";

export default function Sidebar() {
  const pathname = usePathname();

  const navItems = [
    { href: "/", label: "Home" },
    { href: "/about", label: "About" },
    { href: "/tutorial", label: "Tutorial" },
  ];

  return (
    <header className="z-40">
      <div
        className="mx-auto flex w-full max-w-7xl items-center justify-center px-4 py-3 sm:px-6"
        style={{
          color: "#dbefe7",
        }}
      >
        <nav className="flex items-center justify-center gap-1 sm:gap-2">
          {navItems.map((item) => {
            const isActive = pathname === item.href;
            return (
              <Link
                key={item.href}
                href={item.href}
                aria-current={isActive ? "page" : undefined}
                className="px-4 py-2 text-sm font-medium transition-all hover:opacity-100"
                style={{
                  backgroundColor: isActive ? "rgba(219, 239, 231, 0.14)" : "transparent",
                  color: "#dbefe7",
                  opacity: isActive ? 1 : 0.78,
                }}
              >
                {item.label}
              </Link>
            );
          })}
        </nav>
      </div>
    </header>
  );
}
