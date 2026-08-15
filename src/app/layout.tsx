import type { Metadata } from "next";
import { Geist, Geist_Mono } from "next/font/google";
import "./globals.css";
import { Toaster } from "@/components/ui/toaster";

const geistSans = Geist({
  variable: "--font-geist-sans",
  subsets: ["latin"],
});

const geistMono = Geist_Mono({
  variable: "--font-geist-mono",
  subsets: ["latin"],
});

export const metadata: Metadata = {
  title: "Helix Codec — Synthetic DNA Data Storage",
  description: "Encode files into synthetic DNA sequences, simulate mutations and sequencing errors, and recover the original data via Reed-Solomon error correction.",
  keywords: ["DNA storage", "synthetic biology", "Reed-Solomon", "error correction", "archival storage", "codec"],
  authors: [{ name: "Helix Codec" }],
  icons: {
    icon: "https://z-cdn.chatglm.cn/z-ai/static/logo.svg",
  },
  openGraph: {
    title: "Helix Codec — Synthetic DNA Data Storage",
    description: "Encode files into synthetic DNA. Simulate mutations. Recover via Reed-Solomon.",
    url: "https://chat.z.ai",
    siteName: "Helix Codec",
    type: "website",
  },
  twitter: {
    card: "summary_large_image",
    title: "Helix Codec — Synthetic DNA Data Storage",
    description: "Encode files into synthetic DNA. Simulate mutations. Recover via Reed-Solomon.",
  },
};

export default function RootLayout({
  children,
}: Readonly<{
  children: React.ReactNode;
}>) {
  return (
    <html lang="en" suppressHydrationWarning>
      <body
        className={`${geistSans.variable} ${geistMono.variable} antialiased bg-background text-foreground`}
      >
        {children}
        <Toaster />
      </body>
    </html>
  );
}
