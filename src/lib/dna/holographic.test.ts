import { describe, it, expect } from 'vitest';
import { holographicShuffle, holographicUnshuffle } from './holographic';

describe('Holographic shuffle', () => {
  it('should be symmetric: unshuffle(shuffle(x)) === x', () => {
    const data = new Uint8Array(64);
    for (let i = 0; i < 64; i++) data[i] = i;
    const shuffled = holographicShuffle(data, 8, 8);
    const unshuffled = holographicUnshuffle(shuffled, 8, 8);
    expect(Array.from(unshuffled)).toEqual(Array.from(data));
  });
});
