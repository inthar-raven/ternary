/**
 * Ternary Scale Analysis Library
 *
 * Provides scale analysis, guide frame detection, and tuning computations
 * for ternary scales.
 */

// ============================================================================
// CONSTANTS
// ============================================================================

const EDO_BOUND = 111;
const S_LOWER_BOUND = 20.0;
const S_UPPER_BOUND = 200.0;

const STEP_LETTERS = [
  "", // 0
  "X", // 1
  "Ls", // 2
  "Lms", // 3
  "Lmns", // 4
  "HLmns", // 5
  "HLmnst", // 6
  "BHLmnst", // 7
  "BHLmnstw", // 8
  "BCHLmnstw", // 9
  "BCHLmnpstw", // 10
  "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz", // >= 11
];

// ============================================================================
// UTILITY FUNCTIONS
// ============================================================================

/**
 * Greatest Common Divisor using Euclidean algorithm
 */
function gcd(a, b) {
  a = Math.abs(a);
  b = Math.abs(b);
  while (b !== 0) {
    [a, b] = [b, a % b];
  }
  return a;
}

/**
 * Modular inverse: find x such that (x*a) mod b = 1
 */
function modinv(a, b) {
  const [g, x, _] = extendedGcd(a, b);
  if (g !== 1) {
    throw new Error("Non-coprime generator error");
  }
  return ((x % b) + b) % b;
}

/**
 * Extended Euclidean algorithm: returns [gcd, x, y] such that ax + by = gcd
 */
function extendedGcd(a, b) {
  let [oldR, r] = [a, b];
  let [oldS, s] = [1, 0];
  let [oldT, t] = [0, 1];

  while (r !== 0) {
    const quotient = Math.floor(oldR / r);
    [oldR, r] = [r, oldR - quotient * r];
    [oldS, s] = [s, oldS - quotient * s];
    [oldT, t] = [t, oldT - quotient * t];
  }
  return [oldR, oldS, oldT];
}

/**
 * Python-style modulo (always positive result)
 */
function mod(n, m) {
  return ((n % m) + m) % m;
}

/**
 * Rotate an array by the given amount
 */
function rotate(arr, amount) {
  if (arr.length === 0) return [...arr];
  const n = mod(amount, arr.length);
  if (n === 0) return [...arr];
  return [...arr.slice(n), ...arr.slice(0, n)];
}

/**
 * Check if two arrays are equal
 */
function arraysEqual(a, b) {
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) {
    if (a[i] !== b[i]) return false;
  }
  return true;
}

/**
 * Compare two arrays lexicographically
 * Returns -1 if a < b, 0 if a == b, 1 if a > b
 */
function compareArrays(a, b) {
  const minLen = Math.min(a.length, b.length);
  for (let i = 0; i < minLen; i++) {
    if (a[i] < b[i]) return -1;
    if (a[i] > b[i]) return 1;
  }
  return a.length - b.length;
}

/**
 * Get all rotations of a word
 */
function rotations(word) {
  return Array.from({ length: word.length }, (_, i) => rotate(word, i));
}

/**
 * Booth's algorithm for lexicographically least rotation
 * Returns the rotation index
 */
function booth(scale) {
  const n = scale.length;
  const f = new Array(2 * n).fill(-1);
  let k = 0;

  for (let j = 1; j < 2 * n; j++) {
    let i = f[j - k - 1];
    while (i !== -1 && scale[j % n] !== scale[(k + i + 1) % n]) {
      if (scale[j % n] < scale[(k + i + 1) % n]) {
        k = j - i - 1;
      }
      i = f[i];
    }
    if (i === -1 && scale[j % n] !== scale[(k + i + 1) % n]) {
      if (scale[j % n] < scale[(k + i + 1) % n]) {
        k = j;
      }
      f[j - k] = -1;
    } else {
      f[j - k] = i + 1;
    }
  }
  return k;
}

/**
 * Get the lexicographically least mode of a word
 */
function leastMode(scale) {
  return rotate(scale, booth(scale));
}

// ============================================================================
// COUNT VECTOR (Multiset Operations)
// ============================================================================

/**
 * A CountVector is a multiset represented as a Map from elements to counts.
 */
class CountVector {
  constructor(map = new Map()) {
    this.map = new Map(map);
  }

  static ZERO = new CountVector();

  static fromSlice(arr) {
    const map = new Map();
    for (const elem of arr) {
      map.set(elem, (map.get(elem) || 0) + 1);
    }
    return new CountVector(map);
  }

  static fromObject(obj) {
    return new CountVector(
      new Map(Object.entries(obj).map(([k, v]) => [parseInt(k), v])),
    );
  }

  add(other) {
    const result = new Map(this.map);
    for (const [key, value] of other.map) {
      const newVal = (result.get(key) || 0) + value;
      if (newVal === 0) {
        result.delete(key);
      } else {
        result.set(key, newVal);
      }
    }
    return new CountVector(result);
  }

  neg() {
    const result = new Map();
    for (const [key, value] of this.map) {
      result.set(key, -value);
    }
    return new CountVector(result);
  }

  scalarMul(lambda) {
    const result = new Map();
    for (const [key, value] of this.map) {
      result.set(key, lambda * value);
    }
    return new CountVector(result);
  }

  get(key) {
    return this.map.get(key) || 0;
  }

  len() {
    let sum = 0;
    for (const value of this.map.values()) {
      sum += Math.abs(value);
    }
    return sum;
  }

  isEmpty() {
    return this.len() === 0;
  }

  keys() {
    return [...this.map.keys()].sort((a, b) => a - b);
  }

  toObject() {
    const obj = {};
    for (const [key, value] of this.map) {
      obj[key] = value;
    }
    return obj;
  }

  toArray() {
    const maxKey = Math.max(...this.keys(), 0);
    const arr = new Array(maxKey + 1).fill(0);
    for (const [key, value] of this.map) {
      arr[key] = value;
    }
    return arr;
  }

  equals(other) {
    if (this.map.size !== other.map.size) return false;
    for (const [key, value] of this.map) {
      if (other.get(key) !== value) return false;
    }
    return true;
  }

  toString() {
    return JSON.stringify(this.toObject());
  }
}

/**
 * Get the interval (as CountVector) spanning a subword
 */
function wordOnDegree(scale, degree, subwordLength) {
  const rotated = rotate(scale, degree);
  if (subwordLength <= scale.length) {
    return rotated.slice(0, subwordLength);
  } else {
    const fullCopies = Math.floor(subwordLength / scale.length);
    const remainder = subwordLength % scale.length;
    let result = [];
    for (let i = 0; i < fullCopies; i++) {
      result = result.concat(rotated);
    }
    return result.concat(rotated.slice(0, remainder));
  }
}

/**
 * Get the dyad on a degree as a CountVector
 */
function dyadOnDegree(scale, degree, intervalClass) {
  return CountVector.fromSlice(wordOnDegree(scale, degree, intervalClass));
}

/**
 * Get the distinct spectrum of subwords of a given length
 */
function distinctSpectrum(scale, subwordLength) {
  const seen = new Set();
  const result = [];
  for (let degree = 0; degree < scale.length; degree++) {
    const cv = dyadOnDegree(scale, degree, subwordLength);
    const key = cv.toString();
    if (!seen.has(key)) {
      seen.add(key);
      result.push(cv);
    }
  }
  return result;
}

/**
 * Compute the maximum variety of a scale
 */
function maximumVariety(scale) {
  if (scale.length === 0) return 0;
  let maxVar = 0;
  const floorHalf = Math.floor(scale.length / 2);
  for (let subwordLength = 1; subwordLength <= floorHalf; subwordLength++) {
    const spectrum = distinctSpectrum(scale, subwordLength);
    maxVar = Math.max(maxVar, spectrum.length);
  }
  return maxVar;
}

// ============================================================================
// MOS (Moment of Symmetry) SCALE GENERATION
// ============================================================================

/**
 * Generate the darkest mode of a MOS scale and its dark generator
 * Using Bresenham's line algorithm
 */
function darkestMosModeAndGenBresenham(a, b) {
  const d = gcd(a, b);
  if (d === 1) {
    const countGenerSteps = modinv(a, a + b);
    const resultScale = [];
    let [currentX, currentY] = [0, 0];

    while (currentX !== a || currentY !== b) {
      if (a * currentY >= b * (currentX + 1)) {
        currentX += 1;
        resultScale.push(0);
      } else {
        currentY += 1;
        resultScale.push(1);
      }
    }

    const resultGener = CountVector.fromSlice(
      resultScale.slice(0, countGenerSteps),
    );
    return [resultScale, resultGener];
  } else {
    const [primMos, gener] = darkestMosModeAndGenBresenham(a / d, b / d);
    const repeated = [];
    for (let i = 0; i < d; i++) {
      repeated.push(...primMos);
    }
    return [repeated, gener];
  }
}

// ============================================================================
// WORD/SCALE OPERATIONS
// ============================================================================

/**
 * Get the step variety (number of distinct steps) in a scale
 */
function stepVariety(scale) {
  return new Set(scale).size;
}

/**
 * Replace all instances of one letter with another
 */
function replace(scale, from, to) {
  return scale.map((x) => (x === from ? to : x));
}

/**
 * Delete all instances of a letter
 */
function deleteStep(scale, letter) {
  return scale.filter((x) => x !== letter);
}

/**
 * Check if L=m results in a MOS (monotone property)
 */
function monotoneLm(scale) {
  return maximumVariety(replace(scale, 1, 0)) === 2;
}

/**
 * Check if m=s results in a MOS (monotone property)
 */
function monotoneMs(scale) {
  return maximumVariety(replace(scale, 2, 1)) === 2;
}

/**
 * Check if s=0 results in a MOS (monotone property)
 */
function monotoneS0(scale) {
  return maximumVariety(deleteStep(scale, 2)) === 2;
}

/**
 * Letterwise substitution for scale words
 */
function subst(template, x, filler) {
  if (filler.length === 0) {
    return template.filter((letter) => letter !== x);
  }

  const result = [];
  let i = 0;
  for (const letter of template) {
    if (letter === x) {
      result.push(filler[i % filler.length]);
      i++;
    } else {
      result.push(letter);
    }
  }
  return result;
}

/**
 * Generate MOS substitution scales for one permutation
 */
function mosSubstitutionScalesOnePerm(n0, n1, n2) {
  const [template, _] = darkestMosModeAndGenBresenham(n0, n1 + n2);
  const [filler, gener] = darkestMosModeAndGenBresenham(n1, n2);
  const fillerShifted = filler.map((x) => x + 1);
  const generSize = gener.len();

  const results = [];
  for (let i = 0; i < n1 + n2; i++) {
    const rotatedFiller = rotate(
      fillerShifted,
      (i * generSize) % fillerShifted.length,
    );
    results.push(subst(template, 1, rotatedFiller));
  }
  return results;
}

/**
 * Generate all MOS substitution scales for a given step signature
 */
function mosSubstitutionScales(sig) {
  const [n0, n1, n2] = sig;

  // n0L (n1m n2s)
  const scales1 = mosSubstitutionScalesOnePerm(n0, n1, n2);

  // n1m (n0L n2s) -> map letters
  const scales2 = mosSubstitutionScalesOnePerm(n1, n2, n0).map((scale) =>
    scale.map((x) => (x + 1) % 3),
  );

  // n2s (n0L n1m) -> map letters
  const scales3 = mosSubstitutionScalesOnePerm(n2, n0, n1).map((scale) =>
    scale.map((x) => (x === 0 ? 2 : (x - 1) % 3)),
  );

  const allScales = [...scales1, ...scales2, ...scales3];

  // Canonicalize and deduplicate
  const seen = new Set();
  const unique = [];
  for (const scale of allScales) {
    const canonical = leastMode(scale);
    const key = canonical.join(",");
    if (!seen.has(key)) {
      seen.add(key);
      unique.push(canonical);
    }
  }

  return unique.sort((a, b) => compareArrays(a, b));
}

/**
 * Check if scale is a MOS substitution scale
 */
function isMosSubst(scale, t, f1, f2) {
  return (
    stepVariety(scale) === 3 &&
    maximumVariety(deleteStep(scale, t)) === 2 &&
    maximumVariety(replace(scale, f1, f2)) === 2
  );
}

/**
 * Get the chirality of a scale word
 */
function chirality(word) {
  const leastModeWord = leastMode(word);
  const wordRev = [...word].reverse();
  const leastModeWordRev = leastMode(wordRev);

  const cmp = compareArrays(leastModeWord, leastModeWordRev);
  if (cmp < 0) return "Right";
  if (cmp > 0) return "Left";
  return "Achiral";
}

// ============================================================================
// GUIDE FRAME ANALYSIS
// ============================================================================

/**
 * Get stacked k-steps for a scale
 */
function stackedKSteps(k, neck) {
  return Array.from({ length: neck.length }, (_, i) => {
    const subword = wordOnDegree(neck, k * i, k);
    return CountVector.fromSlice(subword);
  });
}

/**
 * Get guided GS chains from a chain of intervals
 */
function guidedGsChains(chain) {
  const len = chain.length;
  const results = [];

  for (const rotation of rotations(chain)) {
    // Check if the last element is different from all previous elements
    const last = rotation[len - 1];
    const prefix = rotation.slice(0, len - 1);

    let foundDuplicate = false;
    for (const cv of prefix) {
      if (cv.equals(last)) {
        foundDuplicate = true;
        break;
      }
    }

    if (!foundDuplicate) {
      results.push(weakPeriodCountVector(prefix));
    }
  }

  return results;
}

/**
 * Get weak period for CountVector arrays
 */
function weakPeriodCountVector(arr) {
  for (let l = 1; l <= arr.length; l++) {
    const prefix = arr.slice(0, l);
    let matches = true;
    for (let i = 0; i < arr.length; i++) {
      if (!arr[i].equals(prefix[i % l])) {
        matches = false;
        break;
      }
    }
    if (matches) {
      return prefix;
    }
  }
  return arr;
}

/**
 * Get all k-step guided GS list
 */
function kStepGuidedGsList(k, neck) {
  return guidedGsChains(stackedKSteps(k, neck));
}

/**
 * Get all guided generator sequences for a scale
 */
function guidedGsList(neck) {
  const len = neck.length;
  const results = [];

  for (let k = 2; k <= len - 2; k++) {
    if (gcd(k, len) === 1) {
      results.push(...kStepGuidedGsList(k, neck));
    }
  }

  return results;
}

/**
 * Get offset between two CountVector sequences if they are rotationally equivalent
 */
function offsetVec(vec1, vec2) {
  if (vec1.length !== vec2.length) return null;

  for (let i = 0; i < vec1.length; i++) {
    const rotated = rotate(vec2, i);
    let matches = true;
    for (let j = 0; j < vec1.length; j++) {
      if (!vec1[j].equals(rotated[j])) {
        matches = false;
        break;
      }
    }
    if (matches) {
      return i;
    }
  }
  return null;
}

/**
 * GuideFrame class representing a guide frame structure
 */
class GuideFrame {
  constructor(gs, polyoffset = [CountVector.ZERO]) {
    this.gs = gs;
    this.polyoffset = polyoffset;
  }

  complexity() {
    return this.gs.length * this.polyoffset.length;
  }

  multiplicity() {
    return this.polyoffset.length;
  }

  static trySimple(scale, k) {
    if (scale.length === 0 || gcd(scale.length, k) !== 1) {
      return [];
    }

    return kStepGuidedGsList(k, scale).map(
      (gs) => new GuideFrame(gs, [CountVector.ZERO]),
    );
  }

  static tryMultiple(scale, multiplicity, k) {
    if (
      multiplicity === 1 ||
      scale.length === 0 ||
      scale.length % multiplicity !== 0
    ) {
      return [];
    }

    const d = gcd(k, scale.length);
    const coD = scale.length / d;

    if (coD % multiplicity !== 0) {
      if (d === multiplicity) {
        // Interleaved scale
        const subscales = [];
        for (let degree = 0; degree < d; degree++) {
          const rotation = rotate(scale, degree);
          const stacked = stackedKSteps(d, rotation);
          subscales.push(stacked.slice(0, scale.length / d));
        }

        const subscaleOnRoot = subscales[0];
        const offsets = [];

        for (let i = 0; i < subscales.length; i++) {
          const offset = offsetVec(subscaleOnRoot, subscales[i]);
          if (offset === null) return [];
          offsets.push(
            CountVector.fromSlice(wordOnDegree(scale, 0, offset * d + i)),
          );
        }

        offsets.sort((a, b) => a.len() - b.len());

        if (offsets.length === 1 && offsets[0].equals(CountVector.ZERO)) {
          const gsList = guidedGsListForSubscale(subscaleOnRoot);
          return gsList.map((gs) => new GuideFrame(gs, [CountVector.ZERO]));
        }

        const gsList = guidedGsListForSubscale(subscaleOnRoot);
        return gsList.map((gs) => new GuideFrame(gs, offsets));
      }
    }

    return [];
  }
}

/**
 * Get guided GS list for a subscale (represented as CountVector array)
 */
function guidedGsListForSubscale(subscale) {
  if (subscale.length === 2) {
    return [[subscale[0]]];
  }

  const len = subscale.length;
  const results = [];

  for (let k = 2; k <= len / 2; k++) {
    if (gcd(k, len) === 1) {
      const stacked = [];
      for (let i = 0; i < len; i++) {
        // Stack k intervals
        let sum = CountVector.ZERO;
        for (let j = 0; j < k; j++) {
          sum = sum.add(subscale[(i * k + j) % len]);
        }
        stacked.push(sum);
      }
      results.push(...guidedGsChains(stacked));
    }
  }

  return results;
}

/**
 * Get all guide frames for a scale
 */
function guideFrames(scale) {
  const results = [];
  const len = scale.length;

  // Simple guide frames
  for (let k = 2; k <= len - 2; k++) {
    if (gcd(k, len) === 1) {
      results.push(...GuideFrame.trySimple(scale, k));
    }
  }

  // Multiple/interleaved guide frames
  const factors = [];
  for (let m = 2; m <= Math.sqrt(len); m++) {
    if (len % m === 0) {
      factors.push(m);
      if (m !== len / m) {
        factors.push(len / m);
      }
    }
  }

  for (const m of factors) {
    for (let k = 2; k <= len - 2; k++) {
      results.push(...GuideFrame.tryMultiple(scale, m, k));
    }
  }

  // Sort by complexity
  results.sort((a, b) => a.complexity() - b.complexity());

  // Remove duplicates
  const seen = new Set();
  const unique = [];
  for (const frame of results) {
    const key = JSON.stringify({
      gs: frame.gs.map((cv) => cv.toString()),
      polyoffset: frame.polyoffset.map((cv) => cv.toString()),
    });
    if (!seen.has(key)) {
      seen.add(key);
      unique.push(frame);
    }
  }

  return unique;
}

/**
 * Convert a GuideFrame to a result object
 * Returns objects with string keys "0", "1", "2" to match the expected format
 */
function guideFrameToResult(structure) {
  const { gs, polyoffset } = structure;

  // Convert CountVector to object format with string keys
  const cvToObj = (cv) => ({
    0: cv.get(0) || 0,
    1: cv.get(1) || 0,
    2: cv.get(2) || 0,
  });

  const gsObjects = gs.map(cvToObj);

  const aggregate = gs.reduce((acc, cv) => acc.add(cv), CountVector.ZERO);
  const aggregateObj = cvToObj(aggregate);

  const polyoffsetObjects = polyoffset.map(cvToObj);

  return {
    gs: gsObjects,
    aggregate: aggregateObj,
    polyoffset: polyoffsetObjects,
    multiplicity: structure.multiplicity(),
    complexity: structure.complexity(),
  };
}

/**
 * Compute 3x3 determinant for lattice basis check
 * Works with both arrays and objects with string keys "0", "1", "2"
 */
function det3(v0, v1, v2) {
  // Handle both array format and object format
  const get = (v, i) => (Array.isArray(v) ? v[i] : (v[String(i)] ?? v[i] ?? 0));

  return (
    get(v0, 0) * get(v1, 1) * get(v2, 2) +
    get(v0, 1) * get(v1, 2) * get(v2, 0) +
    get(v0, 2) * get(v1, 0) * get(v2, 1) -
    get(v0, 2) * get(v1, 1) * get(v2, 0) -
    get(v0, 1) * get(v1, 0) * get(v2, 2) -
    get(v0, 0) * get(v1, 2) * get(v2, 1)
  );
}

/**
 * Convert step vector to object format
 */
function toStepVectorObj(v) {
  if (Array.isArray(v)) {
    return { 0: v[0] || 0, 1: v[1] || 0, 2: v[2] || 0 };
  }
  return { 0: v["0"] || 0, 1: v["1"] || 0, 2: v["2"] || 0 };
}

/**
 * Get unimodular basis if one exists
 */
function getUnimodularBasis(structures, stepSig) {
  for (const structure of structures) {
    if (structure.multiplicity() === 1) {
      const result = guideFrameToResult(structure);
      const gs = result.gs;

      // Check pairs from gs
      for (let i = 0; i < gs.length; i++) {
        for (let j = i; j < gs.length; j++) {
          if (Math.abs(det3(stepSig, gs[i], gs[j])) === 1) {
            return [toStepVectorObj(gs[i]), toStepVectorObj(gs[j]), result];
          }
        }
      }

      // Check polyoffset with gs
      const polyoffset = result.polyoffset;
      for (const v of polyoffset) {
        for (const w of gs) {
          if (Math.abs(det3(stepSig, v, w)) === 1) {
            return [toStepVectorObj(v), toStepVectorObj(w), result];
          }
        }
      }
    } else {
      const result = guideFrameToResult(structure);
      const vecForGsElement = result.gs[0];
      const polyoffset = result.polyoffset;
      if (polyoffset.length > 0) {
        const vecForOffset = polyoffset[polyoffset.length - 1];
        if (Math.abs(det3(stepSig, vecForGsElement, vecForOffset)) === 1) {
          return [
            toStepVectorObj(vecForGsElement),
            toStepVectorObj(vecForOffset),
            result,
          ];
        }
      }
    }
  }
  return null;
}

// ============================================================================
// NECKLACE ENUMERATION (Sawada's Algorithm)
// ============================================================================

/**
 * Generate all necklaces with a given content (step signature)
 * Using Sawada's algorithm
 */
function necklacesFixedContent(content) {
  if (content.length === 0) return [];

  // Filter out zeros and track the permutation
  const nonZeroContent = [];
  const letterMap = [];
  for (let i = 0; i < content.length; i++) {
    if (content[i] > 0) {
      letterMap.push(i);
      nonZeroContent.push(content[i]);
    }
  }

  if (nonZeroContent.length === 0) return [];

  const arity = nonZeroContent.length;
  const scaleLen = nonZeroContent.reduce((a, b) => a + b, 0);

  // Initialize
  const remContent = [...nonZeroContent];
  remContent[0] -= 1;

  const word = [0];
  for (let i = 1; i < scaleLen; i++) {
    word.push(arity - 1);
  }

  const availLetters = [];
  for (let i = arity - 1; i >= 0; i--) {
    if (i === 0 && remContent[0] === 0) continue;
    availLetters.push(i);
  }

  const runs = new Array(scaleLen).fill(0);
  const collection = [];

  sawadaMut(
    remContent,
    runs,
    availLetters,
    word,
    1,
    1,
    1,
    collection,
    arity,
    scaleLen,
  );

  // Map back to original letters
  return collection.map((necklace) =>
    necklace.map((letter) => letterMap[letter]),
  );
}

/**
 * Recursive Sawada algorithm implementation
 */
function sawadaMut(
  remContent,
  runs,
  availLetters,
  a,
  t,
  p,
  s,
  currColl,
  arity,
  scaleLen,
) {
  if (remContent[arity - 1] === scaleLen - t) {
    if (
      (remContent[arity - 1] === runs[t - p] && scaleLen % p === 0) ||
      remContent[arity - 1] > runs[t - p]
    ) {
      currColl.push([...a]);
    }
  } else if (remContent[0] !== scaleLen - t) {
    if (availLetters.length > 0) {
      let jIdx = 0;
      let j = availLetters[jIdx];

      while (j >= a[t - p]) {
        runs[s] = t - s;

        const wasLastOfKind = remContent[j] === 1;
        if (wasLastOfKind) {
          const idx = availLetters.indexOf(j);
          if (idx !== -1) availLetters.splice(idx, 1);
        }
        remContent[j] -= 1;
        a[t] = j;

        sawadaMut(
          remContent,
          runs,
          availLetters,
          a,
          t + 1,
          j === a[t - p] ? p : t + 1,
          j === arity - 1 ? s : t + 1,
          currColl,
          arity,
          scaleLen,
        );

        if (wasLastOfKind) {
          // Find where to insert j back (maintaining descending order)
          let insertIdx = 0;
          while (
            insertIdx < availLetters.length &&
            availLetters[insertIdx] > j
          ) {
            insertIdx++;
          }
          availLetters.splice(insertIdx, 0, j);
        }
        remContent[j] += 1;

        // Find next smaller available letter
        jIdx++;
        while (jIdx < availLetters.length && availLetters[jIdx] >= j) {
          jIdx++;
        }
        if (jIdx >= availLetters.length) break;
        j = availLetters[jIdx];
      }
    }
    a[t] = arity - 1;
  }
}

// ============================================================================
// SCALE CONVERSION UTILITIES
// ============================================================================

/**
 * Convert a string scale word to array of numbers
 */
function stringToNumbers(word) {
  const result = [];
  const distinctChars = new Set(word);
  const arity = distinctChars.size;
  const letters = STEP_LETTERS[Math.min(arity, 11)];

  for (const c of word) {
    const idx = letters.indexOf(c);
    if (idx !== -1) {
      result.push(idx);
    }
  }
  return result;
}

/**
 * Convert array of numbers to string scale word
 */
function numbersToString(word) {
  const distinctNums = new Set(word);
  const arity = distinctNums.size;
  const letters = STEP_LETTERS[Math.min(arity, 11)];

  let result = "";
  for (const i of word) {
    if (i < letters.length) {
      result += letters[i];
    }
  }
  return result;
}

/**
 * Get step signature from a word
 */
function wordToSig(input) {
  const result = [0, 0, 0];
  for (const i of input) {
    if (i < 3) {
      result[i]++;
    }
  }
  return result;
}

// ============================================================================
// EQUAL TUNINGS
// ============================================================================

/**
 * Convert steps in ed(equave) to cents
 */
function stepsAsCents(steps, ed, equaveCents = 1200) {
  return (steps / ed) * equaveCents;
}

/**
 * Find all ED tunings for a ternary step signature
 */
function edTuningsForTernary(
  stepSig,
  edBound = EDO_BOUND,
  aberLower = S_LOWER_BOUND,
  aberUpper = S_UPPER_BOUND,
  equaveCents = 1200,
) {
  const results = [];

  for (let l = 3; l < edBound; l++) {
    for (let m = 2; m < l; m++) {
      for (let s = 1; s < m; s++) {
        const ed = l * stepSig[0] + m * stepSig[1] + s * stepSig[2];
        if (ed <= edBound) {
          const aberSize = stepsAsCents(s, ed, equaveCents);
          if (aberLower <= aberSize && aberSize <= aberUpper) {
            results.push([l, m, s]);
          }
        }
      }
    }
  }

  return results;
}

// ============================================================================
// MAIN API FUNCTIONS
// ============================================================================

/**
 * Generate a scale profile for a word
 */
function wordToProfile(query) {
  const brightest = numbersToString(leastMode(query));
  const lm = monotoneLm(query);
  const ms = monotoneMs(query);
  const s0 = monotoneS0(query);
  const chi = chirality(query);
  const reversed = numbersToString(leastMode([...query].reverse()));
  const mv = maximumVariety(query);
  const stepSig = wordToSig(query);
  const substLMs = isMosSubst(query, 0, 1, 2);
  const substMLs = isMosSubst(query, 1, 0, 2);
  const substSLm = isMosSubst(query, 2, 0, 1);

  const frames = guideFrames(query);
  const basisResult = getUnimodularBasis(frames, stepSig);

  if (basisResult) {
    const [basis1, basis2, structure] = basisResult;
    return {
      word: brightest,
      lattice_basis: [basis1, basis2],
      chirality: chi,
      reversed: reversed,
      structure: structure,
      lm: lm,
      ms: ms,
      s0: s0,
      subst_l_ms: substLMs,
      subst_m_ls: substMLs,
      subst_s_lm: substSLm,
      mv: mv,
    };
  } else {
    return {
      word: brightest,
      lattice_basis: null,
      chirality: chi,
      reversed: reversed,
      structure: null,
      lm: lm,
      ms: ms,
      s0: s0,
      subst_l_ms: substLMs,
      subst_m_ls: substMLs,
      subst_s_lm: substSLm,
      mv: mv,
    };
  }
}

/**
 * Get ED tunings formatted as strings
 */
function sigToEdTunings(stepSig) {
  const edTunings = edTuningsForTernary(stepSig);
  return edTunings.map((v) => {
    const edo = v[0] * stepSig[0] + v[1] * stepSig[1] + v[2] * stepSig[2];
    return v.map((i) => `${i}\\${edo}`);
  });
}

/**
 * Process a step signature query (main API function)
 */
function sigResult(
  query,
  {
    lm = false,
    ms = false,
    s0 = false,
    ggsLen = 0,
    ggsLenConstraint = "at most",
    complexity = 0,
    complexityConstraint = "at most",
    mv = 0,
    mvConstraint = "at most",
    mosSubst = "off",
  } = {},
) {
  const stepSig = query;

  // Get all scales with this signature
  let scales;
  if (mosSubst === "on") {
    scales = mosSubstitutionScales(stepSig);
  } else {
    scales = necklacesFixedContent(stepSig);
  }

  // Apply filters
  const filteringCond = (scale) => {
    if (lm && !monotoneLm(scale)) return false;
    if (ms && !monotoneMs(scale)) return false;
    if (s0 && !monotoneS0(scale)) return false;

    if (ggsLen > 0) {
      const frames = guideFrames(scale);
      if (frames.length === 0) return false;
      if (ggsLenConstraint === "exactly") {
        if (frames[0].gs.length !== ggsLen) return false;
      } else {
        if (frames[0].gs.length > ggsLen) return false;
      }
    }

    if (mv > 0) {
      const scaleMaxVar = maximumVariety(scale);
      if (mvConstraint === "exactly") {
        if (scaleMaxVar !== mv) return false;
      } else {
        if (scaleMaxVar > mv) return false;
      }
    }

    if (complexity > 0) {
      const frames = guideFrames(scale);
      if (frames.length === 0) return false;
      if (complexityConstraint === "exactly") {
        if (frames[0].complexity() !== complexity) return false;
      } else {
        if (frames[0].complexity() > complexity) return false;
      }
    }

    return true;
  };

  scales = scales.filter(filteringCond);

  // Generate profiles
  const profiles = scales.map((scale) => wordToProfile(scale));

  // Sort by complexity
  profiles.sort((a, b) => {
    const complexityA = a.structure ? a.structure.complexity : 255;
    const complexityB = b.structure ? b.structure.complexity : 255;
    return complexityA - complexityB;
  });

  return {
    profiles: profiles,
    ji_tunings: [], // JI tunings are disabled for now (would require more complex computation)
    ed_tunings: sigToEdTunings(stepSig),
  };
}

/**
 * Process a word query (main API function)
 */
function wordResult(query) {
  const wordAsNumbers = stringToNumbers(query);
  const stepSig = wordToSig(wordAsNumbers);

  return {
    profile: wordToProfile(wordAsNumbers),
    ji_tunings: [],
    ed_tunings: sigToEdTunings(stepSig),
  };
}

// ============================================================================
// EXPORTS
// ============================================================================

// Export for use in browser (window) or module systems
const TernaryLib = {
  // Constants
  EDO_BOUND,
  S_LOWER_BOUND,
  S_UPPER_BOUND,
  STEP_LETTERS,

  // Utility functions
  gcd,
  modinv,
  extendedGcd,
  mod,
  rotate,
  arraysEqual,
  compareArrays,
  rotations,
  booth,
  leastMode,

  // CountVector class
  CountVector,

  // Word operations
  wordOnDegree,
  dyadOnDegree,
  distinctSpectrum,
  maximumVariety,

  // MOS generation
  darkestMosModeAndGenBresenham,

  // Scale properties
  stepVariety,
  replace,
  deleteStep,
  monotoneLm,
  monotoneMs,
  monotoneS0,
  subst,
  mosSubstitutionScales,
  isMosSubst,
  chirality,

  // Guide frames
  stackedKSteps,
  guidedGsChains,
  kStepGuidedGsList,
  guidedGsList,
  GuideFrame,
  guideFrames,
  guideFrameToResult,
  det3,
  toStepVectorObj,
  getUnimodularBasis,

  // Necklace enumeration
  necklacesFixedContent,

  // Conversion utilities
  stringToNumbers,
  numbersToString,
  wordToSig,

  // ED tunings
  stepsAsCents,
  edTuningsForTernary,

  // Main API
  wordToProfile,
  sigToEdTunings,
  sigResult,
  wordResult,
};

// Export for different module systems
if (typeof module !== "undefined" && module.exports) {
  module.exports = TernaryLib;
}
if (typeof window !== "undefined") {
  window.TernaryLib = TernaryLib;
}

// ES module export
export default TernaryLib;
