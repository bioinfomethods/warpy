<script setup lang="ts">
import { computedAsync } from "@vueuse/core";
//import BamAlignmentPlot from "./components/BamAlignmentPlot.vue";
import AlignmentPlot from "./components/AlignmentPlot.vue";
import { Options } from "./components/options";
import { makeSegment, RawSegment, Segment } from "./components/segment";
import { computed, Ref } from "vue";
import { computeSha1 } from "./components/utils";
import * as d3 from "d3";

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

async function fetchFileData(file: File): Promise<string> {
  const blob = await file.text();
  return blob;
}

const selectedFile = defineModel<File | null>("selectedFile");

const rawText: Ref<string | undefined> = computedAsync(() => {
  const file = selectedFile.value;
  if (file != undefined) {
    return fetchFileData(file);
  }
});

const signature = computedAsync(() => {
  return computeSha1(rawText.value || "");
});

const rawData: Ref<Map<string,Segment[]> | undefined> = computed(() => {
  if (rawText.value) {
    const raw: { [locus: string]: RawSegment[] } = JSON.parse(rawText.value);
    
    const res: Map<string,Segment[]> = new Map();
    for (const locus in raw) {
      const rawSegs = raw[locus] || [];
      res.set(locus, d3.map(rawSegs, makeSegment))
    }
    return res;
  }
});

const currentLocus = defineModel<string>("currentLocus");
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
</script>

<template>
  <div>
    <v-file-input
      v-model="selectedFile"
      label="source data"
      hint="select a json file with alignment segments"
      accept="text/json"
    ></v-file-input>
    <v-card v-if="rawData">
      <v-tabs v-model="currentLocus">
        <v-tab v-for="[locus, _segments] in rawData" :key="signature + locus" :value="locus">{{ locus }}</v-tab>
      </v-tabs>
      <v-tabs-window v-model="currentLocus">
        <v-tabs-window-item v-for="[locus, segments] in rawData" :key="signature + locus" :value="locus">
          <AlignmentPlot :locus="locus as string" :segments="segments" :options="options" />
        </v-tabs-window-item>
      </v-tabs-window>
    </v-card>
    <!--
      <v-select label="Sample" :items="samples" v-model="selectedSample"></v-select>
      <v-text-field
        clearable
        label="Locus/Loci"
        hint="space separated for multiple"
        placeholder="chr7:6007774-6098582"
        :rules="[validLoci]"
        v-model="rawLoci"
      ></v-text-field>
      <BamAlignmentPlot v-if="selectedBam" :bam-url="selectedBam" :bam-headers="bamHeaders" :loci="loci" :options="options" />
    -->
    <v-card> </v-card>
  </div>
</template>
