const wasm = require("../src/lib/dna/wasm-pkg/helix_dna_wasm.js");
console.log("wasm exports:", Object.keys(wasm));
console.log("typeof LdpcCode:", typeof wasm.LdpcCode);
console.log("typeof full_decode:", typeof wasm.full_decode);
console.log("typeof test_arithmetic_decode:", typeof wasm.test_arithmetic_decode);
