import js from "@eslint/js";
import tsPlugin from "@typescript-eslint/eslint-plugin";
import tsParser from "@typescript-eslint/parser";
import reactPlugin from "eslint-plugin-react";
import globals from "globals";

export default [
  {
    ignores: [
      "node_modules/**",
      ".next/**",
      "dist/**",
      "rust/**",
      "wasm-src/**",
      "scripts/**",
      "test-data/**",
      "src/lib/dna/wasm-pkg/**",
      "src/lib/dna/pkg/**",
      "*.lock",
      "bun.lock",
      "package-lock.json",
    ],
  },
  {
    files: ["**/*.{ts,tsx,js,jsx}"],
    plugins: {
      "@typescript-eslint": tsPlugin,
      "react": reactPlugin,
    },
    languageOptions: {
      parser: tsParser,
      parserOptions: {
        ecmaVersion: "latest",
        sourceType: "module",
        ecmaFeatures: { jsx: true },
      },
      globals: {
        ...globals.browser,
        ...globals.node,
        ...globals.es2021,
        React: "readonly",
      },
    },
    rules: {
      ...js.configs.recommended.rules,
      // Warnings only (not errors) for code quality issues
      "no-unused-vars": "warn",
      "no-console": "off",
      "no-empty": "warn",
      "no-self-assign": "warn",
      // TypeScript-specific
      "@typescript-eslint/no-explicit-any": "off",
      "@typescript-eslint/no-unused-vars": ["warn", { argsIgnorePattern: "^_" }],
      "@typescript-eslint/no-require-imports": "off",
      "@typescript-eslint/ban-ts-comment": "off",
      // React
      "react/react-in-jsx-scope": "off",
    },
  },
];
