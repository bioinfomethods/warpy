import { depermute } from "./utils";

export type Comparator<T> = (a: T, b: T) => number;

export interface ColumnSpecificationStr<T> {
  kind: "str";
  title: string;
  field: keyof T;
  value: (item: T) => string;
  sort?: Comparator<T>;
}

export interface ColumnSpecificationNum<T> {
  kind: "num";
  title: string;
  field: keyof T;
  value: (item: T) => number;
  sort?: Comparator<T>;
}

export type ColumnSpecification<T> = ColumnSpecificationStr<T> | ColumnSpecificationNum<T>;

function cmpStr(a: string, b: string): number {
  if (a < b) {
    return -1;
  }
  if (a > b) {
    return +1;
  }
  return 0;
}

function cmpNum(a: number, b: number): number {
  return a - b;
}

export function makeComparator<T>(spec: ColumnSpecification<T>): Comparator<T> {
  if (spec.sort) {
    return spec.sort;
  }
  if (spec.kind == "str") {
    return (a: T, b: T) => {
      const ak = spec.value(a);
      const bk = spec.value(b);
      return cmpStr(ak, bk);
    };
  } else {
    return (a: T, b: T) => {
      const ak = spec.value(a);
      const bk = spec.value(b);
      return cmpNum(ak, bk);
    };
  }
}

export function maybeFlipComparator<T>(cmp: Comparator<T>, ascending: boolean): Comparator<T> {
  if (ascending) {
    return cmp;
  } else {
    return (a: T, b: T) => cmp(b, a);
  }
}

export function sortInPlace<T>(xs: T[], spec: ColumnSpecification<T>, ascending: boolean) {
  const cmp = maybeFlipComparator(makeComparator(spec), ascending);

  function permCmp(i: number, j: number): number {
    return cmp(xs[i], xs[j]) || (ascending ? i - j : j - i);
  }

  const perm: number[] = new Array(xs.length);
  for (let i = 0; i < xs.length; ++i) {
    perm[i] = i;
  }
  perm.sort(permCmp);
  depermute(xs, perm);
}
