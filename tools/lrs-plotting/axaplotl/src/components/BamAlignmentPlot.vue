<script setup lang="ts">
import { computed, ComputedRef } from "vue";
import { computedAsync } from "@vueuse/core";
import { BamFile } from "@gmod/bam";
import * as d3 from "d3";
import { optionDefaults, Options } from "./options";
import ChromPlot from "./ChromPlot.vue";
import { Locus, makeSegment, RawSegment, Segment } from "./segment";
import { scanSegments } from "./scanner";
import { computeSha1 } from "./utils";

const props = defineProps<{
  bam: BamFile;
  loci: Locus[];
  options: Partial<Options>;
}>();

const locusString = computed(() => {
  return d3.map(props.loci, (locus) => `${locus.chrom}:${locus.start}-${locus.end}`).join(",");
});

const segments = computedAsync(async () => {
  console.log("scanning...");
  const segs = await scanSegments(props.bam, props.loci);
  console.log(segs);
  return segs;
});

const options = computed(() => {
  return { ...optionDefaults, ...props.options };
});

const chroms = computed(() => {
  if (!segments.value) {
    return [];
  }
  const seen: Set<string> = new Set();
  segments.value.forEach((seg) => {
    seen.add(seg.chrom);
  });
  const res = Array.from(seen);
  res.sort();
  return res;
});

const dataByChrom: ComputedRef<d3.InternMap<string, Segment[]>> = computed(() =>
  d3.group(
    d3.map(segments.value, (d) => makeSegment(d)),
    (d) => d.chrom
  )
);

const readColours = computed(() => {
  const groups = d3.group(segments.value, (seg) => seg.readid);
  let nextColour = 0;
  const res: Map<string, string> = new Map();
  groups.forEach((_segs, readid) => {
    const colour = nextColour++ / groups.size;
    res.set(readid, d3.interpolateTurbo(colour));
  });
  return res;
});

const chromtabs = defineModel("chromtabs");

const jsonBlob = computed<string>(() => {
  const loc = props.loci[0];
  const data = segments.value || null;
  const res: {[locus:string]: (RawSegment[] | null)} = {};
  res[`${loc.chrom}:${loc.start}-${loc.end}`] = data;
  return JSON.stringify(res);
});

async function saveJsonBlob() {
  const sha = await computeSha1(locusString.value);
  const dataStr = "data:text/json;charset=utf-8," + encodeURIComponent(jsonBlob.value);
  const dlAnchorElem = document.getElementById("downloadAnchor");
  dlAnchorElem?.setAttribute("href", dataStr);
  dlAnchorElem?.setAttribute("download", `${sha.slice(-12)}.json`);
  dlAnchorElem?.click();
}
</script>

<template>
  <div>
    <v-card v-if="segments && segments.length > 0">
      <v-card-title>
        Scanning the <span v-if="loci.length > 1">loci</span><span v-else>locus</span>{{ " " }}
        <span v-for="(locus, i) in loci" :key="`${locus.chrom}:${locus.start}-${locus.end}`">
          <span v-if="i > 0">, </span>{{ `${locus.chrom}:${locus.start}-${locus.end}` }} </span
        >{{ " " }} yields {{ segments?.length || 0 }} alignments.
        <v-btn icon="mdi-download" size="x-small" @click="saveJsonBlob()"></v-btn>
        <a id="downloadAnchor" style="display: none"></a>
      </v-card-title>
      <v-tabs v-model="chromtabs">
        <v-tab v-for="chrom in chroms" :key="chrom" :value="chrom">{{ chrom }}</v-tab>
      </v-tabs>
      <v-tabs-window v-model="chromtabs">
        <v-tabs-window-item v-for="chrom in chroms" :key="chrom" :value="chrom">
          <ChromPlot
            :chrom="chrom"
            :locus="locusString"
            :segments="dataByChrom.get(chrom) || []"
            :reads="[]"
            :options="options"
            :colours="readColours"
          ></ChromPlot>
        </v-tabs-window-item>
      </v-tabs-window>
    </v-card>
    <v-card v-else> </v-card>
  </div>
</template>
