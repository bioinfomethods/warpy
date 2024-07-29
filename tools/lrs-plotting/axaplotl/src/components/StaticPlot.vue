<script setup lang="ts">
import { computedAsync } from "@vueuse/core";
import AlignmentPlot from "./AlignmentPlot.vue";
import { Options } from "./options";
import { makeSegment, RawSegment, Segment } from "./segment";
import { computed, Ref } from "vue";
import { computeSha1 } from "./utils";
import * as d3 from "d3";

defineProps<{
  options: Partial<Options>;
}>();

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

const rawData: Ref<Map<string, Segment[]> | undefined> = computed(() => {
  if (rawText.value) {
    const raw: { [locus: string]: RawSegment[] } = JSON.parse(rawText.value);

    const res: Map<string, Segment[]> = new Map();
    for (const locus in raw) {
      const rawSegs = raw[locus] || [];
      res.set(locus, d3.map(rawSegs, makeSegment));
    }
    return res;
  }
});

const currentLocus = defineModel<string>("currentLocus");
</script>

<template>
  <v-sheet>
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
  </v-sheet>
</template>
