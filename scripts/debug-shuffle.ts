import { holographicShuffle, holographicUnshuffle } from "../src/lib/dna/holographic";

const data = new Uint8Array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10]);
const shuffled = holographicShuffle(data, 1, 10);
const unshuffled = holographicUnshuffle(shuffled, 1, 10);
console.log("Original:  ", Array.from(data).join(","));
console.log("Shuffled:  ", Array.from(shuffled).join(","));
console.log("Unshuffled:", Array.from(unshuffled).join(","));
console.log("Shuffle roundtrip:", JSON.stringify(Array.from(data)) === JSON.stringify(Array.from(unshuffled)));

// Check if shuffle is bijective (no duplicate positions)
const positions = new Set<number>();
for (let i = 0; i < 10; i++) {
  // We can't directly call feistelPermute since it's not exported, but we can check
  // the shuffled array has all the same elements
}
const sortedShuffled = Array.from(shuffled).sort((a,b) => a-b);
const sortedOriginal = Array.from(data).sort((a,b) => a-b);
console.log("Same elements:", JSON.stringify(sortedShuffled) === JSON.stringify(sortedOriginal));
