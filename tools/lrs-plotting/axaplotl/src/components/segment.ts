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
  return `${seg.readid}-${seg.chrom}-${seg.pos}`;
}
