export type Segment = {
  readid: string;
  chrom: string;
  pos: number;
  strand: string;
  qual: number;
  offset: number;
  rlen: number;
  qlen: number;
};

export function slug(seg: Segment): string {
  return `${seg.readid}-${seg.chrom}-${seg.pos}-${seg.strand}-${seg.qual}-${seg.offset}-${seg.rlen}-${seg.qlen}`;
}

export type SegmentGroupInfo = {
  grp: string;
  grpMin: number;
  grpMax: number;
  grpWidth: number;
  ratio: number;
  winMin: number;
  winMax: number;
  winWidth: number;
  ticks: number;
  segments: Segment[];
};

export type ReadInfo = {
  start: number;
  length: number;
  begin: number;
  flip: boolean;
};

export type Locus = {
  chrom: string;
  start: number;
  end: number;
};

export type ReadItem = {
  selected: boolean;
  colour: string;
  readid: string;
  mapped: Locus[];
};

export function validLocus(txt: string): string | boolean {
  const m = (txt || "").match(/^(chr([0-9]+|X|Y|MY)):([0-9]+)[-]([0-9]+)$/);
  if (!m) {
    return false;
  }
  //const chrom = m[1];
  const begin = m[3];
  const end = m[4];
  if (end != undefined) {
    if (parseInt(begin) > parseInt(end)) {
      return "end must be greater than start";
    }
  }
  return true;
}

export function parseLocus(txt: string): Locus | null {
  const m = (txt || "").match(/^(chr([0-9]+|X|Y|MY)):([0-9]+)[-]([0-9]+)$/);
  if (!m) {
    return null;
  }
  const chrom = m[1];
  const start = parseInt(m[3]);
  const end = parseInt(m[4]);

  return { chrom: chrom, start: start, end: end };
}
