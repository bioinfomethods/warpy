<script setup lang="ts">
import StaticPlot from "./components/StaticPlot.vue";
import DynamicPlot from "./components/DynamicPlot.vue";
import { Options } from "./components/options";

const options: Partial<Options> = {};

/*
const manifest = computedAsync(async () => {
  const res: Map<string, string> = new Map();
  const resp = await fetch("/manifest.json");
  const stuff = await resp.json();
  if (typeof stuff == "object") {
    const keys = Object.keys(stuff);
    for (const key of keys) {
      const value = stuff[key];
      if (typeof value == "string") {
        res.set(key, value);
      } else {
        console.log(`warning: dropping ${key} from manifest.json`);
      }
    }
  }
  return res;
});

const samples = computed(() => {
  if (manifest.value) {
    const res: string[] = [...manifest.value.keys()];
    return res;
  } else {
    const res: string[] = [];
    return res;
  }
});

const selectedSample = defineModel<string>();

const selectedBam = computed(() => {
  if (selectedSample.value && manifest.value) {
    return manifest.value.get(selectedSample.value);
  } else {
    return undefined;
  }
});

function validLocus(txt: string): string | boolean {
  const m = (txt || "").match(/^(chr([0-9]+|X|Y|MY)):([0-9]+)[-]([0-9]+)$/);
  if (!m) {
    return false;
  }
  //const chrom = m[1];
  const begin = m[3];
  const end = m[4];
  if (end != undefined) {
    if (parseInt(begin) > parseInt(end)) return "end must be greater than start";
  }
  return true;
}

function validLoci(txt: string): string | boolean {
  const parts = (txt || "").split(/[ ]+/);
  for (const part of parts) {
    const r = validLocus(part.trim());
    if (r != true) {
      return r;
    }
  }
  return true;
}
*/

/*
const rawLoci = defineModel<string>("rawLoci", { default: "chr21:44416000-44423001" });

const loci = computed(() => {
  const res: [string, number, number][] = [];
  (rawLoci.value || "").split(/[ ]+/).forEach((s) => {
    const m = s.match(/^(chr([0-9]+|X|Y|MY)):([0-9]+)[-]([0-9]+)$/);
    if (m) {
      const chrom = m[1];
      const begin = parseInt(m[3]);
      const end = parseInt(m[4]);
      res.push([chrom, begin, end]);
    }
  });
  return res;
});
*/

const staticOrDynamic = defineModel<string>("staticOrDynamic");
</script>

<template>
  <v-sheet>
    <v-tabs v-model="staticOrDynamic">
      <v-tab value="static">Precomputed Alignment Summaries</v-tab>
      <v-tab value="dynamic">Alignment Summaries from BAM</v-tab>
    </v-tabs>
    <v-tabs-window v-model="staticOrDynamic">
      <v-tabs-window-item value="static">
        <StaticPlot :options="options"></StaticPlot>
      </v-tabs-window-item>
      <v-tabs-window-item value="dynamic">
        <DynamicPlot :options="options"></DynamicPlot>
      </v-tabs-window-item>
    </v-tabs-window>
  </v-sheet>
</template>
