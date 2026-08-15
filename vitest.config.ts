import { defineConfig } from 'vitest/config';
import path from 'path';

export default defineConfig({
  test: {
    globals: true,
    environment: 'node',
    include: ['src/**/*.test.ts'],
    coverage: {
      provider: 'v8',
      reporter: ['text', 'json', 'html'],
      include: ['src/lib/dna/**/*.ts'],
      exclude: ['src/lib/dna/pkg/**', 'src/lib/dna/wasm-pkg/**', 'src/lib/dna/**/*.test.ts', 'src/lib/dna/types.ts'],
    },
    timeout: 60000,
  },
  resolve: {
    alias: {
      '@': path.resolve(__dirname, './src'),
    },
  },
});
