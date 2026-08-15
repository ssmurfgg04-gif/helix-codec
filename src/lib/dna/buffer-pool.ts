/**
 * Buffer Pool — Reuse TypedArrays to eliminate GC pressure
 *
 * Allocating new Uint8Array for each read is expensive (GC pressure).
 * This pool reuses buffers of fixed sizes.
 */

export class BufferPool {
  private pools: Map<number, Uint8Array[]> = new Map();
  private maxSize: number;

  constructor(maxSize: number = 10000) {
    this.maxSize = maxSize;
  }

  /**
   * Acquire a buffer of the given size.
   * Returns a reused buffer if available, otherwise allocates new.
   */
  acquire(size: number): Uint8Array {
    const pool = this.pools.get(size);
    if (pool && pool.length > 0) {
      return pool.pop()!;
    }
    return new Uint8Array(size);
  }

  /**
   * Release a buffer back to the pool for reuse.
   */
  release(buffer: Uint8Array) {
    const size = buffer.length;
    let pool = this.pools.get(size);
    if (!pool) {
      pool = [];
      this.pools.set(size, pool);
    }
    if (pool.length < this.maxSize) {
      pool.push(buffer);
    }
  }

  /**
   * Clear all pooled buffers.
   */
  clear() {
    this.pools.clear();
  }

  /**
   * Get pool statistics.
   */
  stats(): { sizes: number[]; totalBuffers: number; totalBytes: number } {
    let totalBuffers = 0;
    let totalBytes = 0;
    const sizes: number[] = [];
    for (const [size, pool] of this.pools) {
      sizes.push(size);
      totalBuffers += pool.length;
      totalBytes += pool.length * size;
    }
    return { sizes, totalBuffers, totalBytes };
  }
}

// Global buffer pool
export const globalBufferPool = new BufferPool();
