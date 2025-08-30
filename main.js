
const fs = require('fs');

// ---------- BigInt utilities ----------
const ZERO = 0n, ONE = 1n, NEG_ONE = -1n;
const abs = (x) => (x < 0n ? -x : x);

function bgcd(a, b) {
  a = abs(a); b = abs(b);
  while (b !== 0n) {
    const t = b;
    b = a % b;
    a = t;
  }
  return a;
}

// ---------- Rational (BigInt fraction) ----------
class Frac {
  constructor(n, d = 1n) {
    if (d === 0n) throw new Error("Zero denominator");
    if (d < 0n) { n = -n; d = -d; }
    const g = bgcd(n, d);
    this.n = n / g;
    this.d = d / g;
  }
  static fromBigInt(bi) { return new Frac(bi, 1n); }
  add(o) { return new Frac(this.n*o.d + o.n*this.d, this.d*o.d); }
  sub(o) { return new Frac(this.n*o.d - o.n*this.d, this.d*o.d); }
  mul(o) { return new Frac(this.n*o.n, this.d*o.d); }
  div(o) { if (o.n === 0n) throw new Error("Divide by zero"); return new Frac(this.n*o.d, this.d*o.n); }
  inv()  { if (this.n === 0n) throw new Error("Inverse of zero"); return new Frac(this.d, this.n); }
  isZero(){ return this.n === 0n; }
  toString() { return this.d === 1n ? this.n.toString() : `${this.n.toString()}/${this.d.toString()}`; }
}

// ---------- parse mixed-base string to BigInt ----------
function charToVal(ch) {
  const c = ch.toLowerCase();
  if (c >= '0' && c <= '9') return BigInt(c.charCodeAt(0) - '0'.charCodeAt(0));
  if (c >= 'a' && c <= 'z') return 10n + BigInt(c.charCodeAt(0) - 'a'.charCodeAt(0));
  throw new Error(`Invalid digit '${ch}'`);
}
function parseBase(str, base) {
  let b = BigInt(base);
  let v = 0n;
  for (const ch of str.trim()) {
    const d = charToVal(ch);
    if (d >= b) throw new Error(`Digit ${ch} not valid for base ${base}`);
    v = v * b + d;
  }
  return v;
}

// ---------- solve Vandermonde (Gauss-Jordan) for a0..am where f(x)=sum a_i x^i ----------
function solveCoeffsMonomial(xs, ys) {
  const n = xs.length;
  // Build augmented matrix [A|b], A_ij = x_i^j with j=0..n-1 (monomial basis)
  const A = Array.from({length: n}, (_, i) => Array.from({length: n}, (_, j) => {
    let p = 1n;
    for (let k = 0; k < j; k++) p *= xs[i];
    return new Frac(p, 1n);
  }));
  const b = ys.map(y => new Frac(y, 1n));

  // Gauss-Jordan elimination
  for (let col = 0; col < n; col++) {
    // find pivot
    let piv = col;
    while (piv < n && A[piv][col].isZero()) piv++;
    if (piv === n) throw new Error("Singular system");
    if (piv !== col) {
      [A[col], A[piv]] = [A[piv], A[col]];
      [b[col], b[piv]] = [b[piv], b[col]];
    }
    // normalize pivot row
    const inv = A[col][col].inv();
    for (let j = col; j < n; j++) A[col][j] = A[col][j].mul(inv);
    b[col] = b[col].mul(inv);
    // eliminate other rows
    for (let r = 0; r < n; r++) {
      if (r === col) continue;
      const factor = A[r][col];
      if (!factor.isZero()) {
        for (let j = col; j < n; j++) {
          A[r][j] = A[r][j].sub(factor.mul(A[col][j]));
        }
        b[r] = b[r].sub(factor.mul(b[col]));
      }
    }
  }
  // now A is identity, b holds the solution (a0..a_{n-1})
  return b;
}

// ---------- Lagrange at x=0 (direct secret if you only want f(0)) ----------
function fAtZeroLagrange(xs, ys) {
  // f(0) = sum_i y_i * prod_{j!=i} (-x_j)/(x_i - x_j)
  const n = xs.length;
  let sum = new Frac(0n,1n);
  for (let i = 0; i < n; i++) {
    let num = new Frac(1n,1n), den = new Frac(1n,1n);
    for (let j = 0; j < n; j++) if (j !== i) {
      num = num.mul(new Frac(-xs[j], 1n));
      den = den.mul(new Frac(xs[i] - xs[j], 1n));
    }
    sum = sum.add(new Frac(ys[i],1n).mul(num.div(den)));
  }
  return sum; // Frac
}

// ---------- main ----------
(function main(){
  const input = fs.readFileSync(0, 'utf8');
  const data = JSON.parse(input);

  const n = Number(data.keys.n);
  const k = Number(data.keys.k);

  // Gather all (x,y) from entries other than "keys"
  const points = [];
  for (const key of Object.keys(data)) {
    if (key === 'keys') continue;
    const x = BigInt(key);
    const base = Number(data[key].base);
    const y = parseBase(data[key].value, base);
    points.push({ x, y });
  }

  // Sort by x, take first k points (deterministic choice)
  points.sort((a,b) => (a.x < b.x ? -1 : a.x > b.x ? 1 : 0));
  const take = points.slice(0, k);
  const xs = take.map(p => p.x);
  const ys = take.map(p => p.y);

  // Degree = k-1
  const degree = k - 1;

  // Solve coefficients (exact)
  const coeffs = solveCoeffsMonomial(xs, ys); // [a0..am] as Frac

  // Also, f(0) via Lagrange (handy if only the "secret" is required)
  const secret = fAtZeroLagrange(xs, ys);

  // Print results
  console.log(`degree: ${degree}`);
  console.log(`coefficients a0..a${degree} (reduced fractions):`);
  coeffs.forEach((c, i) => console.log(`a${i} = ${c.toString()}`));
  console.log(`f(0) = ${secret.toString()}`);

  // Optional: print integer-scaled coefficients using common denominator
  const denoms = coeffs.map(c => c.d);
  const lcm = denoms.reduce((L, d) => (L / bgcd(L, d)) * d, 1n);
  const scaled = coeffs.map(c => (c.n * (lcm / c.d)));
  console.log(`common_denominator = ${lcm.toString()}`);
  console.log(`integer_coeffs_for_(common_denominator * f(x)) = [`);
  console.log('  ' + scaled.map(s => s.toString()).join(',\n  '));
  console.log(']');
})();
