<script setup lang="ts">
import { computed, ComputedRef } from "vue";
import { computedAsync } from "@vueuse/core";
import { BamFile } from "@gmod/bam";
import { RemoteFile, FilehandleOptions } from "generic-filehandle";
import * as d3 from "d3";
import { optionDefaults, Options } from "./options";
import ChromPlot from "./ChromPlot.vue";
import { Segment } from "./segment";
import { scanSegments } from "./scanner";

const props = defineProps<{
  bamUrl: string;
  bamHeaders?: any;
  loci: [string, number, number][];
  options: Partial<Options>;
}>();



const bam = computed(() => {
  const fileOpts: FilehandleOptions = {};
  if (props.bamHeaders) {
    fileOpts.headers = props.bamHeaders;
  }
  const bamFilehandle = new RemoteFile(props.bamUrl, fileOpts);
  const baiFilehandle = new RemoteFile(props.bamUrl + ".bai", fileOpts);
  return new BamFile({ bamFilehandle: bamFilehandle, baiFilehandle: baiFilehandle });
});

const segments = computedAsync(async () => {
  const segs = await scanSegments(bam.value, props.loci);
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

const dataByChrom: ComputedRef<d3.InternMap<string, Segment[]>> = computed(() => d3.group(segments.value, (d) => d.chrom));

const readColours = computed(() => {
  const groups = d3.group(segments.value, (seg) => seg.readid);
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
