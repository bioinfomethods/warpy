const suffixes: string[] = ["", "k", "m", "g", "t", "p", "e"];

/**
 * Convert a number to a corresponding string with an appropriate magnitude suffix.
 * @param x a number to turn into a string
 * @returns a human readable string corresponding to the number with an appropriate magnitude suffix.
 */
export function humanize(x: number): string {
  const neg = (x < 0);
  if (neg) {
    x = -x;
  }
  let n = 0;
  while (x >= 1000) {
    n += 1;
    x /= 1000;
  }
  x = Math.round(x);
  const prefix = (neg? '-' : '');
  const suffix = suffixes[n];
  return `${prefix}${x}${suffix}bp`;
}