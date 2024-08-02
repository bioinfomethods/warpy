import { assert } from "@vueuse/core";

const suffixes: string[] = ["", "k", "m", "g", "t", "p", "e"];

/**
 * Convert a number to a corresponding string with an appropriate magnitude suffix.
 *
 * @param x a number to turn into a string
 * @returns a human readable string corresponding to the number with an appropriate magnitude suffix.
 */
export function humanize(x: number): string {
  const neg = x < 0;
  if (neg) {
    x = -x;
  }
  let n = 0;
  while (x >= 1000) {
    n += 1;
    x /= 1000;
  }
  x = Math.round(x);
  const prefix = neg ? "-" : "";
  const suffix = suffixes[n];
  return `${prefix}${x}${suffix}bp`;
}

/**
 * Compute the SHA-1 digest of a given string.
 *
 * @param text
 * @returns SHA1 digest of the given string.
 */
export async function computeSha1(text: string): Promise<string> {
  const encoder = new TextEncoder();
  const data = encoder.encode(text);
  const hash = await window.crypto.subtle.digest("SHA-1", data);
  const hashArray = Array.from(new Uint8Array(hash)); // convert buffer to byte array
  const hashHex = hashArray.map((b) => b.toString(16).padStart(2, "0")).join(""); // convert bytes to hex string
  return hashHex;
}

export function stableSort<T>(xs: T[], cmp: (a: T, b: T) => number): T[] {
  return xs
    .map((item, index) => ({ item, index }))
    .sort((a, b) => cmp(a.item, b.item) || a.index - b.index)
    .map(({ item }) => item);
}

export function depermute<T>(xs: T[], perm: number[]) {
  assert(xs.length == perm.length);
  for (let i = 0; i < xs.length; ++i) {
    let x = xs[i];
    let j = i;
    while (true) {
      const k = perm[j];
      perm[j] = j;
      if (k == i) {
        break;
      }
      xs[j] = xs[k];
      j = k;
    }
    xs[j] = x;
  }
}
