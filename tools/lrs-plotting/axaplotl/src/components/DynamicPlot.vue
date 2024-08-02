<script setup lang="ts">
/// <reference lib="es2021" />
import { computed, onMounted, ref, watchEffect } from "vue";
import { Options } from "./options";
import { BamFile } from "@gmod/bam";
import { BlobFile } from "generic-filehandle";
import { Locus, makeSegment, parseLocus, Segment } from "./segment";
import AlignmentPlot from "./AlignmentPlot.vue";
import * as d3 from "d3";
import { scanSegments } from "./scanner";

defineProps<{
  options: Partial<Options>;
}>();

function checkFileNames(bamName: string, baiName: string): boolean | string {
  if (!bamName.endsWith(".bam")) {
    return "BAM file must end with .bam";
  }
  if (!baiName.endsWith(".bai")) {
    return "BAI file must end with .bai";
  }
  if (baiName.startsWith(bamName)) {
    // foo.bam, foo.bam.bai case
    return true;
  }
  const bamBaseName = bamName.slice(0, -4);
  const baiBaseName = baiName.slice(0, -4);
  if (bamBaseName != baiBaseName) {
    return `BAM and BAI file names were not corresponding`;
  }
  return true;
}

function validFileSelection(value: File[]): boolean | string {
  if (value.length != 2) {
    return "BAM and BAI file must be selected";
  }
  const file1 = value[0];
  const file2 = value[1];
  if (file1.name.endsWith(".bam")) {
    return checkFileNames(file1.name, file2.name);
  }
  if (file2.name.endsWith(".bam")) {
    return checkFileNames(file2.name, file1.name);
  }
  return "exactly 1 BAM and 1 BAI file must be selected.";
}

const bamAndBaiFile = defineModel<File | File[] | null>("bamAndBaiFile");

const inputFiles = computed<[File, File] | undefined>(() => {
  if (bamAndBaiFile.value && Array.isArray(bamAndBaiFile.value) && validFileSelection(bamAndBaiFile.value) == true) {
    const file1 = bamAndBaiFile.value[0];
    const file2 = bamAndBaiFile.value[1];
    if (file1.name.endsWith(".bam")) {
      const res: [File, File] = [file1, file2];
      return res;
    }
    if (file2.name.endsWith(".bam")) {
      const res: [File, File] = [file2, file1];
      return res;
    }
  }
});

const bamBaseName = computed<string>(() => {
  if (inputFiles.value) {
    const bamFile = inputFiles.value[0];
    return bamFile.name.slice(0, -4);
  }
  return "unknown";
});

const bam = computed<BamFile | undefined>(() => {
  if (inputFiles.value) {
    const bamHandle = new BlobFile(inputFiles.value[0]);
    const baiHandle = new BlobFile(inputFiles.value[1]);
    const res = new BamFile({ bamFilehandle: bamHandle, baiFilehandle: baiHandle });
    return res;
  }
});

const rawLocusString = defineModel<string>("rawLociString", { default: "" });

function validLocus(txt: string): string | boolean {
  const m = (txt || "").match(/^(chr([0-9]+|X|Y|MY)):([0-9,]+)[-]([0-9,]+)$/);
  if (!m) {
    return false;
  }
  //const chrom = m[1];
  const begin = m[3].replaceAll(",", "");
  const end = m[4].replaceAll(",", "");
  if (end != undefined) {
    if (parseInt(begin) > parseInt(end)) return "end must be greater than start";
  }
  return true;
}

const locus = computed<Locus | undefined>(() => {
  if (rawLocusString.value && validLocus(rawLocusString.value) == true) {
    const res = parseLocus(rawLocusString.value.trim().replaceAll(",", ""));
    if (res) {
      return res;
    }
  }
});

const loci = computed<Locus[]>(() => {
  const res: Locus[] = [];
  if (locus.value) {
    const l = locus.value;
    const width = l.end - l.start;
    if (width <= 1000) {
      res.push(l);
    } else {
      const lhs: Locus = { chrom: l.chrom, start: d3.max([1, l.start - 10]) || 1, end: l.start + 10 };
      const rhs: Locus = { chrom: l.chrom, start: d3.max([1, l.end - 10]) || 1, end: l.end + 10 };
      res.push(lhs);
      res.push(rhs);
    }
  }
  return res;
});

const locusString = computed<string | undefined>(() => {
  if (locus.value) {
    return `${locus.value.chrom}:${locus.value.start}-${locus.value.end}`;
  }
});

type LocusAndSegments = {
  locus: Locus;
  loci: Locus[];
  segments: Segment[];
  json: string;
};

const locusAndSegmentsTable = ref<Map<string, LocusAndSegments>>(new Map());

const scanningNow = ref<boolean>(false);

async function doScan() {
  if (bam.value && locus.value && loci.value) {
    const label = `${locus.value.chrom}:${locus.value.start}-${locus.value.end}`;
    scanningNow.value = true;
    console.log("scanning...");
    const segs = await scanSegments(bam.value, loci.value);
    const seen: Set<string> = new Set();
    const segments: Segment[] = [];
    segs.forEach((raw) => {
      const seg = makeSegment(raw);
      if (!seen.has(seg.id)) {
        seen.add(seg.id);
        segments.push(seg);
      }
    });
    const stuff: { [locus: string]: Segment[] } = {};
    stuff[label] = segments;
    const json = JSON.stringify(stuff);
    const item: LocusAndSegments = {
      locus: locus.value,
      loci: loci.value,
      segments: segments,
      json: json,
    };
    locusAndSegmentsTable.value.set(label, item);
    scanningNow.value = false;
  }
}

function removeLocus(locus: string) {
  if (locusAndSegmentsTable.value) {
    locusAndSegmentsTable.value.delete(locus);
  }
}

async function saveJsonBlob(item: LocusAndSegments) {
  if (locusString.value) {
    const loc = `${item.locus.chrom}_${item.locus.start}_${item.locus.end}`;
    const dataStr = "data:text/json;charset=utf-8," + encodeURIComponent(item.json);
    const dlAnchorElem = document.getElementById("downloadAnchor");
    dlAnchorElem?.setAttribute("href", dataStr);
    dlAnchorElem?.setAttribute("download", `${loc}.json`);
    dlAnchorElem?.click();
  }
}

const currentLocus = defineModel<string>("currentLocus");

onMounted(() => {
  watchEffect(() => {
    if (locusAndSegmentsTable.value && locusAndSegmentsTable.value.size > 0) {
      if (!currentLocus.value || !locusAndSegmentsTable.value.has(currentLocus.value)) {
        const locus = Array.from(locusAndSegmentsTable.value.keys())[0];
        currentLocus.value = locus;
      }
    }
  });
});
</script>

<template>
  <v-sheet>
    <v-card>
      <v-card-text>
        <v-file-input
          v-model="bamAndBaiFile"
          label="BAM & BAI"
          hint="select a BAM and corresponding BAI file with long read alignments"
          accept=".bam,.bai"
          multiple
          :rules="[validFileSelection]"
        ></v-file-input>
        <v-text-field v-model="rawLocusString" :rules="[validLocus]" label="Breakpoint locus" hint="Locus to scan (if large, select breakpoints)" clearable></v-text-field>
        <v-btn :disabled="!locus || !bam" @click="doScan()">Scan Alignments</v-btn>
        <span :style="{ visibility: scanningNow ? 'visible' : 'hidden', 'margin-left': '1rem' }">
          <v-progress-circular indeterminate></v-progress-circular>
        </span>
      </v-card-text>
      <v-card-text>
        <v-tabs v-model="currentLocus">
          <v-tab v-for="[locus, item] in locusAndSegmentsTable" :key="bamBaseName + locus" :value="locus">
            {{ locus }}
            <v-btn icon="mdi-download" size="x-small" variant="plain" @click="saveJsonBlob(item)"></v-btn>
            <v-btn icon="mdi-close" size="x-small" variant="plain" @click="removeLocus(locus)"></v-btn>
          </v-tab>
        </v-tabs>
        <v-tabs-window v-model="currentLocus">
          <v-tabs-window-item v-for="[locus, item] in locusAndSegmentsTable" :key="bamBaseName + locus" :value="locus">
            <AlignmentPlot :locus="locus" :segments="item.segments" :options="options"></AlignmentPlot>
          </v-tabs-window-item>
        </v-tabs-window>
      </v-card-text>
    </v-card>
    <a id="downloadAnchor" style="display: none"></a>
  </v-sheet>
</template>
