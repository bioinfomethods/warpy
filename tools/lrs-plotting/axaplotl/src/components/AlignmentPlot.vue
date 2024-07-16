<script setup lang="ts">
import { computed, ComputedRef } from "vue";

import * as d3 from "d3";
import { optionDefaults, Options } from "./options";
import ChromPlot from "./ChromPlot.vue";
import { Segment } from "./segment";

const props = defineProps<{
  segments: Segment[];
  options: Partial<Options>;
}>();


const options = computed(() => {
  return { ...optionDefaults, ...props.options };
});

const chroms = computed(() => {
  if (!props.segments) {
    return [];
  }
  const seen: Set<string> = new Set();
  props.segments.forEach((seg) => {
    seen.add(seg.chrom);
  });
  const res = Array.from(seen);
  res.sort();
  return res;
});

const dataByChrom: ComputedRef<d3.InternMap<string, Segment[]>> = computed(() => d3.group(props.segments, (d) => d.chrom));

const readColours = computed(() => {
  const groups = d3.group(props.segments, (seg) => seg.readid);
  let nextColour = 0;
  const res: Map<string, number> = new Map();
  groups.forEach((_segs, readid) => {
    const colour = nextColour++ / groups.size;
    res.set(readid, colour);
  });
  return res;
});

const chromtabs = defineModel();
</script>

<template>
  <div>
    <v-card>
      <v-tabs v-model="chromtabs">
        <v-tab v-for="chrom in chroms" :key="chrom" :value="chrom">{{ chrom }}</v-tab>
      </v-tabs>
      <v-tabs-window v-model="chromtabs">
        <v-tabs-window-item v-for="chrom in chroms" :key="chrom" :value="chrom">
          <ChromPlot
          :chrom="chrom"
            :segments="dataByChrom.get(chrom) || []"
            :options="options"
            :colours="readColours"
          ></ChromPlot>
        </v-tabs-window-item>
      </v-tabs-window>
    </v-card>
  </div>
</template>
