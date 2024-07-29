"""Scan BAM files for supplementary alignments.

Usage:
    prepare-alignments.py [options] <bam> <sv-vcf> <zip-name>
    prepare-alignments.py -L [options] <bam> <locus>...

Options:
    -L                      Take one or more loci from the command line and write to stdout.
    --split-indels          Split alignments at insertions/deletions
    --min-indel N           Threshold for splitting indels [default: 50].
"""

import re
import sys
import json
import docopt
import pysam
import zipfile
import hashlib

def parse_locus(locus):
    m = re.match("([^:]+):([0-9,]+)-([0-9,]+)$", locus)
    if m is None:
        return None
    g = m.groups()
    return (g[0], int(g[1].replace(",", "")), int(g[2].replace(",", "")))

def mkpos(s):
    if s == "*":
        return None
    return int(s)

def splitCigar(cig, rc):
    parts = re.findall("([0-9]+)([MIDSH])", cig)
    if rc:
        parts = parts[::-1]
    return parts

def mappedLen(cig):
    t = 0
    for (l, o) in cig:
        if o == "M" or o == "D":
            t += int(l)
    return t

def queryLen(cig):
    t = 0
    for (l, o) in cig:
        if o == "M" or o == "I":
            t += int(l)
    return t

def partition(parts, d):
    breaks = []
    for i in range(len(parts)):
        (n, o) = parts[i]
        if o in 'ID' and n >= d:
            breaks.append(i)
    partitions = []
    j = 0
    for i in breaks:
        partitions.append(parts[j:i])
        partitions.append([parts[i]])
        j = i + 1
    partitions.append(parts[j:])
    return partitions

def lengths(cig):
    p = 0
    q = 0
    for (n, o) in cig:
        if o == 'M':
            p += n
            q += n
        if o == 'I':
            q += n
        if o == 'D':
            p += n
        if o == 'S' or o == 'H':
            q += n
    return (p, q)

def compress(cig):
    match = 0
    delta = 0
    for (n, o) in cig:
        if o == 'M':
            match += n
        if o == 'I':
            delta += n
        if o == 'D':
            delta -= n
    if match == 0:
        if delta == 0:
            return []
        elif delta > 0:
            return [(delta, 'I')]
        else:
            return [(-delta, 'D')]
    if delta == 0:
        return [(match, 'M')]
    elif delta > 0:
        return [(match, 'M'), (delta, 'I')]
    else:
        return [(match, 'M'), (-delta, 'D')]

def split_cigar(cig, d, strand):
    parts = [(int(n), o) for (n,o) in re.findall("([0-9]+)([MIDSH])", cig)]
    if strand == "-":
        parts = parts[::-1]

    (_, readLength) = lengths(parts)

    leftClip = []
    if parts[0][1] == 'S' or parts[0][1] == 'H':
        leftClip = [parts[0]]
        parts = parts[1:]
    
    if parts[-1][1] in 'SH':
        parts = parts[:-1]
    
    if d is not None:
        parts = partition(parts, d)
    else:
        parts = [parts]

    (p, q) = lengths(leftClip)
    res = []
    for piece in parts:
        (p0, q0) = lengths(piece)
        res.append((p, q, ''.join([f'{n}{o}' for (n, o) in compress(piece)]), p0, q0, readLength - (q + q0)))
        p += p0
        q += q0
    return res

def scan_reads_simple(bam, chrom, start, end, seen):
    for a in bam.fetch(chrom, start, end):
        flg = a.flag
        if flg & 4 == 4: # unmapped
            continue
        if flg & 256 == 256: # secondary
            continue
        if not a.has_tag("SA"):
            continue

        qname = a.query_name

        chrom = a.reference_name
        pos = a.reference_start + 1

        slug = (qname, chrom, pos)
        if slug in seen:
            continue
        seen.add(slug)

        strand = "+"
        if flg & 16 > 0:
            strand = "-"
        cig = a.cigarstring
        qual = a.mapping_quality

        items = [(chrom, pos, strand, cig, qual)]

        for sa in a.get_tag("SA").split(";"):
            if sa == "":
                continue
            parts = sa.split(",")
            chrom = parts[0]
            pos = int(parts[1])

            slug = (qname, chrom, pos)
            if slug in seen:
                continue
            seen.add(slug)

            strand = parts[2]
            cig = parts[3]
            qual = int(parts[4])
            items.append((chrom, pos, strand, cig, qual))

        yield (qname, items)

def scan_reads_deep(bam, chrom, start, end, seen):
    res = {}
    todo = {}

    for a in bam.fetch(chrom, start, end):
        flg = a.flag
        if flg & 4 == 4: # unmapped
            continue
        if flg & 256 == 256: # secondary
            continue
        if not a.has_tag("SA"):
            continue

        qname = a.query_name

        chrom = a.reference_name
        pos = a.reference_start + 1
        strand = "+"
        if flg & 16 > 0:
            strand = "-"
        cig = a.cigarstring
        qual = a.mapping_quality

        slug = (qname, chrom, pos)
        if slug in seen:
            continue
        seen.add(slug)

        if qname not in res:
                res[qname] = set()
        item = (chrom, pos, strand, cig, qual)
        res[qname].add(item)

        for sa in a.get_tag("SA").split(";"):
            if sa == "":
                continue
            parts = sa.split(",")
            chrom = parts[0]
            pos = int(parts[1])
            if chrom not in todo:
                todo[chrom] = {}
            if pos not in todo[chrom]:
                todo[chrom][pos] = set()
            todo[chrom][pos].add(qname)

    for chrom in todo:
        for pos in todo[chrom]:
            qnames = todo[chrom][pos]
            for a in bam.fetch(chrom, pos - 1, pos + 1):
                if a.reference_start + 1 != pos:
                    continue
                if a.query_name not in qnames:
                    continue
                flg = a.flag
                if flg & 4 == 4: # unmapped
                    continue
                if flg & 256 == 256: # secondary
                    continue

                slug = (qname, chrom, pos)
                if slug in seen:
                    continue
                seen.add(slug)

                strand = "+"
                if flg & 16 > 0:
                    strand = "-"
                qname = a.query_name
                cig = a.cigarstring
                qual = a.mapping_quality
                item = (chrom, pos, strand, cig, qual)
                res[qname].add(item)

    for qname in res:
        yield (qname, res[qname])

def sig(alleles, svlen):
    h = hashlib.sha256()
    for a in alleles:
        h.update(a.encode())
    h.update(str(svlen).encode())
    return h.hexdigest()[:8]

def parseLocus(txt):
    m = re.match(r'(chr[^:]*):(\d+)-(\d+)', txt)
    if m is None:
        return None
    g = m.groups()
    return (g[0], int(g[1]), int(g[2]))

def mainFromLoci(args):
    bamName = args["<bam>"]
    lociText = args["<locus>"]

    BIG = 1000
    border = 50

    loci = []
    for locusText in lociText:
        locus = parseLocus(locusText)
        if locus is None:
            print(f'could not parse locus "{locusText}"', file=sys.stderr)
            return
        loci.append(locus)

    out = {}
    bam = pysam.AlignmentFile(bamName, "rb")
    for (chrom, start, stop) in loci:
        seen = set()
        items = set()
        if start + border >= stop - border:
            intervals = [(max(1, start - border), stop + border)]
        else:
            intervals = [(max(1, start - border), start + border), (max(1, stop  - border), stop + border)]
        for (ivlStart, ivlEnd) in intervals:
            for item in scan_reads_simple(bam, chrom, ivlStart, ivlEnd, seen):
                (nm, segs) = item
                for seg in sorted(segs):
                    (chrom, pos, strand, cig, qual) = seg
                    for (p, off, subCig, rlen, qlen, r) in split_cigar(cig, None, strand):
                        items.add((nm, chrom, pos + p, strand, qual, off, rlen, qlen))
                if len(items) > BIG:
                    break

        if len(items) == 0:
            # no supplementary mappings!
            # This is a todo case where
            # we might want to split
            # the CIGAR mappings.
            continue

        if len(items) > BIG:
            continue
        items = list(sorted(items))
        for i in range(len(items)):
            (nm, chrom, pos, strand, qual, off, rlen, qlen) = items[i]
            items[i] = {"readid": nm, "chrom": chrom, "pos": pos, "strand": strand, "qual": qual, "offset": off, "rlen": rlen, "qlen": qlen}
        out[f'{chrom}:{start}-{stop}'] = items
    json.dump(out, sys.stdout, indent=2)

def mainShallow(args):
    bamName = args["<bam>"]
    vcfName = args["<sv-vcf>"]
    zipName = args["<zip-name>"]

    BIG = 1000
    border = 50

    nn = 0
    out = zipfile.ZipFile(zipName, mode="w", compression=zipfile.ZIP_DEFLATED)
    bam = pysam.AlignmentFile(bamName, "rb")
    for rec in pysam.VariantFile(vcfName):
        chrom = rec.chrom
        start = rec.pos
        stop = rec.stop
        kind = rec.info["SVTYPE"]
        svlen = 0
        if "SVLEN" in rec.info:
            svlen = rec.info["SVLEN"]
        h = sig(rec.alleles, svlen)
        key = f'{chrom}:{start}-{stop}-{kind}-{h}.json'

        nn += 1
        if nn & 255 == 0:
            print(key)

        seen = set()
        items = set()
        if start + border >= stop - border:
            intervals = [(max(1, start - border), stop + border)]
        else:
            intervals = [(max(1, start - border), start + border), (max(1, stop  - border), stop + border)]
        for (ivlStart, ivlEnd) in intervals:
            for item in scan_reads_simple(bam, chrom, ivlStart, ivlEnd, seen):
                (nm, segs) = item
                for seg in sorted(segs):
                    (chrom, pos, strand, cig, qual) = seg
                    for (p, off, subCig, rlen, qlen, r) in split_cigar(cig, None, strand):
                        items.add((nm, chrom, pos + p, strand, qual, off, rlen, qlen))
                if len(items) > BIG:
                    break

        if len(items) == 0:
            # no supplementary mappings!
            # This is a todo case where
            # we might want to split
            # the CIGAR mappings.
            continue

        if len(items) > BIG:
            print(f'dropping {key}')
            continue
        items = list(sorted(items))
        for i in range(len(items)):
            (nm, chrom, pos, strand, qual, off, rlen, qlen) = items[i]
            items[i] = {"readid": nm, "chrom": chrom, "pos": pos, "strand": strand, "qual": qual, "offset": off, "rlen": rlen, "qlen": qlen}
        out.writestr(key, json.dumps(items))
    out.close()

def mainDeep(args):
    bamName = args["<bam>"]
    vcfName = args["<sv-vcf>"]
    zipName = args["<zip-name>"]

    splitThreshold = int(args["--min-indel"])

    BIG = 1000
    border = 50

    nn = 0
    out = zipfile.ZipFile(zipName, mode="w", compression=zipfile.ZIP_DEFLATED)
    bam = pysam.AlignmentFile(bamName, "rb")
    for rec in pysam.VariantFile(vcfName):
        chrom = rec.chrom
        start = rec.pos
        stop = rec.stop
        kind = rec.info["SVTYPE"]
        svlen = 0
        if "SVLEN" in rec.info:
            svlen = rec.info["SVLEN"]
        h = sig(rec.alleles, svlen)
        key = f'{chrom}:{start}-{stop}-{kind}-{h}.json'
        print(f'{key} {svlen}')
        nn += 1
        if nn & 255 == 0:
            print(key)
            break

        seen = set()
        items = set()
        if start + border >= stop - border:
            intervals = [(max(1, start - border), stop + border)]
        else:
            intervals = [(max(1, start - border), start + border), (max(1, stop  - border), stop + border)]
        for (ivlStart, ivlEnd) in intervals:
            for item in scan_reads_deep(bam, chrom, ivlStart, ivlEnd, seen):
                (nm, segs) = item
                for seg in sorted(segs):
                    (chrom, pos, strand, cig, qual) = seg
                    for (p, off, subCig, rlen, qlen, r) in split_cigar(cig, splitThreshold, strand):
                        items.add((nm, chrom, pos + p, strand, qual, off, rlen, qlen))
                if len(items) > BIG:
                    break

        if len(items) == 0:
            # no supplementary mappings!
            # This is a todo case where
            # we might want to split
            # the CIGAR mappings.
            continue

        if len(items) > BIG:
            print(f'dropping {key}')
            continue
        items = list(sorted(items))
        for i in range(len(items)):
            (nm, chrom, pos, strand, qual, off, rlen, qlen) = items[i]
            items[i] = {"readid": nm, "chrom": chrom, "pos": pos, "strand": strand, "qual": qual, "offset": off, "rlen": rlen, "qlen": qlen}
        out.writestr(key, json.dumps(items))
    out.close()

def main(args):
    if args['-L']:
        return mainFromLoci(args)
    if args['--split-indels']:
        return mainDeep(args)
    mainShallow(args)

if __name__ == "__main__":
    args = docopt.docopt(__doc__)
    main(args)
