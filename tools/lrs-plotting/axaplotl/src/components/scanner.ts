import { BamFile } from "@gmod/bam";
import { Locus, makeSegment, RawSegment } from "./segment";

/**
 * Options for scanning aligned segments.
 */
export type ScannerOptions = {
  /**
   * Maximum number of reads to examine. Setting this to 0 allows
   * the scanner to scan an unlimited number, which could lead to
   * performance issues. (Default: 200)
   */
  maxReads: number;

  /**
   * Usually, it is undesirable to include reads that have only a
   * single aligned segment, so we drop them. Setting this to true
   * will cause reads with a single segment to be included.
   * (Default: false)
   */
  includeSingletons: boolean;
};

export const defaultScannerOptions: ScannerOptions = {
  maxReads: 200,
  includeSingletons: false,
};

function lengths(parts: [string, number][]): [number, number] {
  let p = 0;
  let q = 0;
  for (const [o, n] of parts) {
    switch (o) {
      case "M": {
        p += n;
        q += n;
        break;
      }
      case "I": {
        q += n;
        break;
      }
      case "D": {
        p += n;
        break;
      }
      case "S":
      case "H": {
        q += n;
        break;
      }
    }
  }
  return [p, q];
}

function compress(parts: [string, number][]): [string, number][] {
  let match = 0;
  let delta = 0;
  for (const [o, n] of parts) {
    switch (o) {
      case "M": {
        match += n;
        break;
      }
      case "I": {
        delta += n;
        break;
      }
      case "D": {
        delta -= n;
        break;
      }
    }
  }
  if (delta == 0) {
    return [["M", match]];
  } else if (delta > 0) {
    return [
      ["M", match],
      ["I", delta],
    ];
  } else {
    return [
      ["M", match],
      ["D", -delta],
    ];
  }
}

function splitCigar(cig: string, strand: string): [number, number, string, number, number][] {
  const parts: [string, number][] = Array.from(cig.matchAll(/(\d+)(\D)/g)).map((m) => {
    return [m[2].toUpperCase(), Number.parseInt(m[1])];
  });
  if (strand == "-") {
    parts.reverse();
  }
  const leftClip: [string, number][] = [];
  if ((parts.length > 0 && parts[0][0] === "S") || parts[0][0] === "H") {
    const item = parts.shift();
    if (item !== undefined) {
      leftClip.push(item);
    }
  }

  if (parts[parts.length - 1][0] === "S" || parts[parts.length - 1][0] === "H") {
    parts.pop();
  }

  let [p, q] = lengths(leftClip);
  let [p0, q0] = lengths(parts);
  return [
    [
      p,
      q,
      compress(parts)
        .map(([o, n]) => `${n}${o}`)
        .join(""),
      p0,
      q0,
    ],
  ];
}

export async function scanSegments(bam: BamFile, loci: Locus[], opts: Partial<ScannerOptions> = {}): Promise<RawSegment[]> {
  const options: ScannerOptions = { ...defaultScannerOptions, ...opts };

  await bam.getHeader();
  const bamOpts = {};
  const resItems: Set<[string, string, number, string, number, number, number, number]> = new Set();
  for (const locus of loci) {
    const chrom = locus.chrom;
    const start = locus.start - 1;
    const end = locus.end;
    let recCount: number = 0;
    for await (const recs of bam.streamRecordsForRange(chrom, start, end, bamOpts)) {
      for (const rec of recs) {
        recCount += 1;
        if (options.maxReads > 0 && recCount > options.maxReads) {
          console.log(`too many BAM records for locus ${locus.chrom}:${locus.start}-${locus.end}`);
          break;
        }
        if (rec.isSegmentUnmapped()) {
          continue;
        }
        if (rec.isSecondary()) {
          continue;
        }
        const qname = rec.name();
        const pos = (rec.get("start") as number) + 1;
        const items: Set<[string, number, string, string, number]> = new Set();
        if (true) {
          const strand = rec.isReverseComplemented() ? "-" : "+";
          const cig = rec.cigar() as string;
          const qual = parseInt(rec.qual() || "0");
          items.add([chrom, pos, strand, cig, qual]);
        }
        const wholeSA = rec.get("SA");
        if (wholeSA == undefined) {
          continue;
        }
        for (const sa of wholeSA.split(";")) {
          if (sa === "") {
            continue;
          }
          const parts = sa.split(",");
          const chrom = parts[0];
          const pos = parseInt(parts[1]);
          const strand = parts[2];
          const cig = parts[3];
          const qual = parseInt(parts[4]);
          items.add([chrom, pos, strand, cig, qual]);
        }
        for (const seg of items) {
          const [chrom, pos, strand, cig, qual] = seg;
          for (const [p, off, _subSig, rlen, qlen] of splitCigar(cig, strand)) {
            resItems.add([qname, chrom, pos + p, strand, qual, off, rlen, qlen]);
          }
        }
      }
      if (options.maxReads > 0 && recCount > options.maxReads) {
        break;
      }
    }
  }
  const res: RawSegment[] = [];
  const seen: Set<string> = new Set();
  for (const item of resItems) {
    const [qname, chrom, pos, strand, qual, off, rlen, qlen] = item;
    if (rlen > 0 && qlen > 0) {
      const raw: RawSegment = { readid: qname, chrom: chrom, pos: pos, strand: strand, qual: qual, offset: off, rlen: rlen, qlen: qlen };
      const seg = makeSegment(raw);
      if (seen.has(seg.id)) {
        continue;
      }
      seen.add(seg.id);
      res.push(raw);
    }
  }
  return res;
}
