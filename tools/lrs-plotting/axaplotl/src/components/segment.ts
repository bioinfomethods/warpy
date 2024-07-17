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

export type ReadItem = {
  selected: boolean;
  colour: string;
  readid: string;
}