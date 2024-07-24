<script setup lang="ts">
import { computed, ComputedRef, ref } from "vue";

import * as d3 from "d3";
import { optionDefaults, Options } from "./options";
import ChromPlot from "./ChromPlot.vue";
import ReadTable from "./ReadTable.vue";
import { Locus, ReadItem, Segment } from "./segment";

const props = defineProps<{
  locus: string;
  segments: Segment[];
  options: Partial<Options>;
}>();

const options = computed(() => {
  return { ...optionDefaults, ...props.options };
});

const dataByChrom: ComputedRef<d3.InternMap<string, Segment[]>> = computed(() => d3.group(props.segments, (d) => d.chrom));

const readColours = computed(() => {
  const groups = d3.group(props.segments, (seg) => seg.readid);
  let nextColour = 0;
  const res: Map<string, string> = new Map();
  groups.forEach((_segs, readid) => {
    const colour = nextColour++ / groups.size;
    res.set(readid, d3.interpolateTurbo(colour));
  });
  return res;
});

function constructReadItems(segs: Segment[]): ReadItem[] {
  const lociByReadid: Map<string, Locus[]> = new Map();
  const strands: Map<string, string> = new Map();
  const starts: Map<string, [number, string]> = new Map();
  const ends: Map<string, [number, string]> = new Map();
  const lengths: Map<string, number> = new Map();
  segs.forEach((seg) => {
    const readid: string = seg.readid;
    const locus: Locus = { chrom: seg.chrom, start: seg.pos, end: seg.pos + seg.rlen };
    const loci = lociByReadid.get(readid);
    if (loci) {
      loci.push(locus);
    } else {
      lociByReadid.set(readid, [locus]);
    }

    strands.set(seg.readid, seg.strand);

    const oldStart = starts.get(seg.readid) || [1e20, ""];
    const segStart = seg.offset;
    if (segStart < oldStart[0]) {
      starts.set(seg.readid, [segStart, seg.id]);
    }

    const oldEnd = ends.get(seg.readid) || [-1, ""];
    const segEnd = seg.offset + seg.qlen;
    if (segEnd > oldEnd[0]) {
      ends.set(seg.readid, [segEnd, seg.id]);
    }

    const oldLength = lengths.get(seg.readid) || 0;
    if (segEnd > oldLength) {
      lengths.set(seg.readid, segEnd);
    }
  });

  function getId(ids: Map<string, [number, string]>, readid: string): string {
    const item = ids.get(readid);
    if (item) {
      return item[1];
    }
    return "";
  }

  const res: ReadItem[] = [];
  readColours.value.forEach((colour, readid) => {
    res.push({
      selected: true,
      colour: colour,
      readid: readid,
      flipped: false,
      strand: strands.get(readid) || "+",
      start: getId(starts, readid),
      end: getId(ends, readid),
      length: lengths.get(readid) || 0,
      mapped: lociByReadid.get(readid) || [],
    });
  });
  return res;
}

const readItems = ref<ReadItem[]>(constructReadItems(props.segments));

function collectReadIds(segs: Segment[]): string[] {
  const seen: Set<string> = new Set();
  segs.forEach((seg) => {
    seen.add(seg.readid);
  });
  const res = Array.from(seen);
  res.sort();
  return res;
}

const selectedReadIds = ref<string[]>(collectReadIds(props.segments));

const selectedDataByChrom: ComputedRef<Map<string, Segment[]>> = computed(() => {
  const wantedReads: Set<string> = new Set();
  selectedReadIds.value.forEach((readid) => {
    wantedReads.add(readid);
  });

  const res: Map<string, Segment[]> = new Map();
  dataByChrom.value.forEach((segs, locus) => {
    const selectedSegments: Segment[] = d3.filter(segs, (seg) => wantedReads.has(seg.readid));
    if (selectedSegments.length > 0) {
      res.set(locus, selectedSegments);
    }
  });
  return res;
});

const chromtabs = defineModel("chromtabs");
</script>

<template>
  <v-sheet rounded border>
    <v-tabs v-model="chromtabs">
      <v-tab v-for="chrom in selectedDataByChrom.keys()" :key="chrom" :value="chrom">{{ chrom }}</v-tab>
    </v-tabs>
    <v-tabs-window v-model="chromtabs">
      <v-tabs-window-item v-for="chrom in selectedDataByChrom.keys()" :key="chrom" :value="chrom">
        <ChromPlot
          :locus="locus"
          :chrom="chrom"
          :segments="selectedDataByChrom.get(chrom) || []"
          :reads="readItems"
          :options="options"
          :colours="readColours"
        ></ChromPlot>
      </v-tabs-window-item>
    </v-tabs-window>
    <ReadTable v-model:selected="selectedReadIds" v-model:items="readItems"></ReadTable>
  </v-sheet>
</template>
