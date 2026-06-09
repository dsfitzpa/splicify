import { TopNav } from "../components/splicify/Landing";
import Chat from "../components/Chat";

export const metadata = {
  title: "Splicify · Clone v1.0",
  description: "Splicify dashboard — plasmid design, solved.",
};

export default function DashboardPage() {
  return (
    <div
      style={{
        background: "var(--forest-800)",
        color: "var(--mint-200)",
        minHeight: "100vh",
        fontFamily: "var(--font-body)",
      }}
    >
      <TopNav variant="dark" active="engine" />
      <main
        style={{
          maxWidth: 1200,
          width: "100%",
          margin: "0 auto",
          padding: "40px 40px 56px",
        }}
      >
        <Chat />
      </main>
    </div>
  );
}
