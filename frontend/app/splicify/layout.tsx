import type { Metadata } from "next";

export const metadata: Metadata = {
  title: "Splicify — AI Molecular Biologist",
  description:
    "Design plasmids and CRISPR experiments with Claude Sonnet 4.6. Drop a .gb file, describe what you want, get assembled.gb + parts + protocol.",
};

export default function SplicifyLayout({ children }: { children: React.ReactNode }) {
  return <div data-splicify>{children}</div>;
}
