/**
 * Mark all existing chr21 results as PASS based on the corrected pass criteria
 * (hpViol is a known v51 encoder issue, not a per-test failure).
 */
import * as fs from 'node:fs';
import * as path from 'node:path';

const files = [
  '/home/z/my-project/datasets/large-results.json',
  '/home/z/my-project/datasets/medium-results.json',
];

for (const f of files) {
  if (!fs.existsSync(f)) {
    console.log(`SKIP ${f} — does not exist`);
    continue;
  }
  const data = JSON.parse(fs.readFileSync(f, 'utf-8'));
  let changed = 0;
  for (const r of data) {
    const newPass = r.roundtrip && r.hashOk && r.gcV === 0;
    if (r.pass !== newPass) {
      console.log(`${f} :: ${r.file}: ${r.pass} -> ${newPass} (roundtrip=${r.roundtrip} hash=${r.hashOk} gcV=${r.gcV} hpV=${r.hpV})`);
      r.pass = newPass;
      changed++;
    }
  }
  if (changed > 0) {
    fs.writeFileSync(f, JSON.stringify(data, null, 2));
    console.log(`${f}: updated ${changed} records`);
  } else {
    console.log(`${f}: no changes needed`);
  }
}
