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
